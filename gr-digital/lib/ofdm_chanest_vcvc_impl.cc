/* -*- c++ -*- */
/*
 * Copyright 2012,2013 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ofdm_chanest_vcvc_impl.h"
#include <gnuradio/io_signature.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace gr {
namespace digital {

ofdm_chanest_vcvc::sptr
ofdm_chanest_vcvc::make(const std::vector<gr_complex>& sync_symbol1,
                        const std::vector<gr_complex>& sync_symbol2,
                        const std::vector<gr_complex>& sync_symbol3,
                        const std::vector<gr_complex>& sync_symbol4,
                        const std::vector<gr_complex>& sync_symbol5,
                        int n_data_symbols,
                        int eq_noise_red_len,
                        int max_carr_offset,
                        bool force_one_sync_symbol)
{
    return gnuradio::get_initial_sptr(new ofdm_chanest_vcvc_impl(sync_symbol1,
                                                                 sync_symbol2,
                                                                 sync_symbol3,
                                                                 sync_symbol4,
                                                                 sync_symbol5,
                                                                 n_data_symbols,
                                                                 eq_noise_red_len,
                                                                 max_carr_offset,
                                                                 force_one_sync_symbol));
}

ofdm_chanest_vcvc_impl::ofdm_chanest_vcvc_impl(
    const std::vector<gr_complex>& sync_symbol1,
    const std::vector<gr_complex>& sync_symbol2,
    const std::vector<gr_complex>& sync_symbol3,
    const std::vector<gr_complex>& sync_symbol4,
    const std::vector<gr_complex>& sync_symbol5,
    int n_data_symbols,
    int eq_noise_red_len,
    int max_carr_offset,
    bool force_one_sync_symbol)
    : block("ofdm_chanest_vcvc",
            io_signature::make2(2, 2, sizeof(gr_complex) * sync_symbol1.size(), sizeof(gr_complex) * sync_symbol1.size()),
            io_signature::make2(2, 2, sizeof(gr_complex) * sync_symbol1.size(), sizeof(gr_complex) * sync_symbol1.size())),
      d_fft_len(sync_symbol1.size()),
      d_n_data_syms(n_data_symbols),
      d_n_sync_syms(1),
      d_eq_noise_red_len(eq_noise_red_len),
      d_ref_sym1_1(sync_symbol2),
      d_ref_sym1_2(sync_symbol3),
      d_ref_sym2_1(sync_symbol4),
      d_ref_sym2_2(sync_symbol5),
      d_corr_v1(sync_symbol2),
      d_corr_v2(sync_symbol4),
      d_known_symbol_diffs(0, 0),
      d_new_symbol_diffs(0, 0),
      d_first_active_carrier(0),
      d_last_active_carrier(sync_symbol2.size() - 1),
      d_interpolate(false)
{
    // Set index of first and last active carrier
    for (int i = 0; i < d_fft_len; i++) {
        if (d_ref_sym1_1[i] != gr_complex(0, 0)) {
            d_first_active_carrier = i;
            break;
        }
    }
    for (int i = d_fft_len - 1; i >= 0; i--) {
        if (d_ref_sym1_1[i] != gr_complex(0, 0)) {
            d_last_active_carrier = i;
            break;
        }
    }

    // Sanity checks
    if (!sync_symbol2.empty()) {
        if (sync_symbol1.size() != sync_symbol2.size()) {
            throw std::invalid_argument("sync symbols must have equal length.");
        }
        if (!force_one_sync_symbol) {
            d_n_sync_syms = 2;
        }
    } else {
        if (sync_symbol1[d_first_active_carrier + 1] == gr_complex(0, 0)) {
            d_last_active_carrier++;
            d_interpolate = true;
        }
    }

    // Set up coarse freq estimation info
    // Allow all possible values:
    d_max_neg_carr_offset1 = -d_first_active_carrier;
    d_max_neg_carr_offset2 = -d_first_active_carrier;
    d_max_pos_carr_offset1 = d_fft_len - d_last_active_carrier - 1;
    d_max_pos_carr_offset2 = d_fft_len - d_last_active_carrier - 1;
    if (max_carr_offset != -1) {
        d_max_neg_carr_offset1 = std::max(-max_carr_offset, d_max_neg_carr_offset1);
        d_max_neg_carr_offset2 = std::max(-max_carr_offset, d_max_neg_carr_offset2);
        d_max_pos_carr_offset1 = std::min(max_carr_offset, d_max_pos_carr_offset1);
        d_max_pos_carr_offset2 = std::min(max_carr_offset, d_max_pos_carr_offset2);
    }
    // Carrier offsets must be even
    if (d_max_neg_carr_offset1 % 2)
        d_max_neg_carr_offset1++;
    if (d_max_pos_carr_offset1 % 2)
        d_max_pos_carr_offset1--;
    if (d_max_neg_carr_offset2 % 2)
        d_max_neg_carr_offset2++;
    if (d_max_pos_carr_offset2 % 2)
        d_max_pos_carr_offset2--;

    if (d_n_sync_syms == 2) {
        for (int i = 0; i < d_fft_len; i++) {
            if (sync_symbol1[i] == gr_complex(0, 0)) {
                d_corr_v1[i] = gr_complex(0, 0);
                d_corr_v2[i] = gr_complex(0, 0);
            } else {
                d_corr_v1[i] /= sync_symbol1[i];
                d_corr_v2[i] /= sync_symbol1[i];
            }
        }
    } else {
        d_corr_v1.resize(0, 0);
        d_corr_v2.resize(0, 0);
        d_known_symbol_diffs.resize(d_fft_len, 0);
        d_new_symbol_diffs.resize(d_fft_len, 0);
        for (int i = d_first_active_carrier;
             i < d_last_active_carrier - 2 && i < d_fft_len - 2;
             i += 2) {
            d_known_symbol_diffs[i] = std::norm(sync_symbol1[i] - sync_symbol1[i + 2]);
        }
    }

    set_output_multiple(d_n_data_syms);
    set_relative_rate((uint64_t)d_n_data_syms, (uint64_t)(d_n_data_syms + 3));
    set_tag_propagation_policy(TPP_DONT);
}

ofdm_chanest_vcvc_impl::~ofdm_chanest_vcvc_impl() {}

void ofdm_chanest_vcvc_impl::forecast(int noutput_items,
                                      gr_vector_int& ninput_items_required)
{
    ninput_items_required[0] =
        (noutput_items / d_n_data_syms) * (d_n_data_syms + 3);
    ninput_items_required[1] =
        (noutput_items / d_n_data_syms) * (d_n_data_syms + 3);
}

int ofdm_chanest_vcvc_impl::get_carr_offset(const gr_complex* sync_sym1,
                                            const gr_complex* sync_sym2,
                                            std::vector<gr_complex> d_corr_v,
                                            int d_max_pos_carr_offset,
                                            int d_max_neg_carr_offset)
{
    int carr_offset = 0;
    if (!d_corr_v.empty()) {
        // Use Schmidl & Cox method
        float Bg_max = 0;
        // g here is 2g in the paper
        for (int g = d_max_neg_carr_offset; g <= d_max_pos_carr_offset; g += 2) {
            gr_complex tmp = gr_complex(0, 0);
            for (int k = 0; k < d_fft_len; k++) {
                if (d_corr_v[k] != gr_complex(0, 0)) {
                    tmp += std::conj(sync_sym1[k + g]) * std::conj(d_corr_v[k]) *
                           sync_sym2[k + g];
                }
            }
            if (std::abs(tmp) > Bg_max) {
                Bg_max = std::abs(tmp);
                carr_offset = g;
            }
        }
    } else {
        // Correlate
        std::fill(d_new_symbol_diffs.begin(), d_new_symbol_diffs.end(), 0);
        for (int i = 0; i < d_fft_len - 2; i++) {
            d_new_symbol_diffs[i] = std::norm(sync_sym1[i] - sync_sym1[i + 2]);
        }

        float sum;
        float max = 0;
        for (int g = d_max_neg_carr_offset; g <= d_max_pos_carr_offset; g += 2) {
            sum = 0;
            for (int j = 0; j < d_fft_len; j++) {
                if (d_known_symbol_diffs[j]) {
                    sum += (d_known_symbol_diffs[j] * d_new_symbol_diffs[j + g]);
                }
                if (sum > max) {
                    max = sum;
                    carr_offset = g;
                }
            }
        }
    }
    return carr_offset;
}


void ofdm_chanest_vcvc_impl::get_chan_taps(const gr_complex* sync_sym1_1,
                                           const gr_complex* sync_sym1_2,
                                           const gr_complex* sync_sym2_1,
                                           const gr_complex* sync_sym2_2,
                                           int carr_offset1,
                                           int carr_offset2,
                                           std::vector<std::vector<std::vector<gr_complex>>>& taps1,
                                           std::vector<std::vector<std::vector<gr_complex>>>& taps2)
{
    const gr_complex* sym1_1 = sync_sym1_1;
    const gr_complex* sym1_2 = sync_sym1_2;
    const gr_complex* sym2_1 = sync_sym2_1;
    const gr_complex* sym2_2 = sync_sym2_2;
    int loop_start1 = 0;
    int loop_start2 = 0;
    int loop_end1 = d_fft_len;
    int loop_end2 = d_fft_len;
    std::vector<std::vector<gr_complex>> S({{gr_complex(0, 0), gr_complex(0, 0)}, {gr_complex(0, 0), gr_complex(0, 0)}});
    std::vector<std::vector<gr_complex>> Y({{gr_complex(0, 0), gr_complex(0, 0)}, {gr_complex(0, 0), gr_complex(0, 0)}});
    if (carr_offset1 > 0) {
        loop_start1 = carr_offset1;
    } else if (carr_offset1 < 0) {
        loop_end1 = d_fft_len + carr_offset1;
    }
    if (carr_offset2 > 0) {
        loop_start2 = carr_offset2;
    } else if (carr_offset2 < 0) {
        loop_end2 = d_fft_len + carr_offset2;
    }

    for (int i = loop_start1; (i < loop_end1) && ((i - carr_offset2) < d_fft_len) && ((i - carr_offset2) > -1) ; i++)
    {
        if ((sync_sym1_1[i - carr_offset1] != gr_complex(0, 0)) && (sync_sym1_2[i - carr_offset1] != gr_complex(0, 0)) && (sync_sym2_1[i - carr_offset2] != gr_complex(0, 0)) && (sync_sym2_2[i - carr_offset2] != gr_complex(0, 0)))
        {
            S[0][0] = d_ref_sym1_1[i - carr_offset1];
            S[0][1] = d_ref_sym1_2[i - carr_offset1];
            S[1][0] = d_ref_sym2_1[i - carr_offset2];
            S[1][1] = d_ref_sym2_2[i - carr_offset2];
            Y[0][0] = sym1_1[i];
            Y[0][1] = sym1_2[i];
            Y[1][0] = sym2_1[i];
            Y[1][1] = sym2_2[i];
            gr_complex** S_ptr = calc_matrix_matrix(S, 2, 2);
            gr_complex** Y_ptr = calc_matrix_matrix(Y, 2, 2);
            gr_complex** S_H_ptr = calc_matrix_H(S_ptr, 2, 2);
            gr_complex** channel_temp = calc_matrix_mult(S_ptr, S_H_ptr, 2, 2, 2);
            channel_temp = calc_matrix_inv(channel_temp, 2);
            channel_temp = calc_matrix_mult(S_H_ptr, channel_temp, 2, 2, 2);
            gr_complex** channel_ptr = calc_matrix_mult(Y_ptr, channel_temp, 2, 2, 2);
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    taps1[i - carr_offset1][j][k] = channel_ptr[j][k];
                }
            }
        }
    }

    for (int i = loop_start2; (i < loop_end2) && ((i - carr_offset1) < d_fft_len) && ((i - carr_offset1) > -1) ; i++)
    {
        if ((sync_sym1_1[i - carr_offset1] != gr_complex(0, 0)) && (sync_sym1_2[i - carr_offset1] != gr_complex(0, 0)) && (sync_sym2_1[i - carr_offset2] != gr_complex(0, 0)) && (sync_sym2_2[i - carr_offset2] != gr_complex(0, 0)))
        {
            S[0][0] = d_ref_sym1_1[i - carr_offset1];
            S[0][1] = d_ref_sym1_2[i - carr_offset1];
            S[1][0] = d_ref_sym2_1[i - carr_offset2];
            S[1][1] = d_ref_sym2_2[i - carr_offset2];
            Y[0][0] = sync_sym1_1[i];
            Y[0][1] = sync_sym1_2[i];
            Y[1][0] = sync_sym2_1[i];
            Y[1][1] = sync_sym2_2[i];
            gr_complex** S_ptr = calc_matrix_matrix(S, 2, 2);
            gr_complex** Y_ptr = calc_matrix_matrix(Y, 2, 2);
            gr_complex** S_H_ptr = calc_matrix_H(S_ptr, 2, 2);
            gr_complex** channel_temp = calc_matrix_mult(S_ptr, S_H_ptr, 2, 2, 2);
            channel_temp = calc_matrix_inv(channel_temp, 2);
            channel_temp = calc_matrix_mult(S_H_ptr, channel_temp, 2, 2, 2);
            gr_complex** channel_ptr = calc_matrix_mult(Y_ptr, channel_temp, 2, 2, 2);
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    taps2[i - carr_offset2][j][k] = channel_ptr[j][k];
                }
            }
        }
    }

    //!! As what we've got, this part will use the integer carrier freq corrected pilot
    //!! to do the channel estimation. So we will focus on the changing of this part.
    /*
    for (int i = loop_start; i < loop_end; i++) {
        if ((d_ref_sym[i - carr_offset] != gr_complex(0, 0))) {
            taps[i - carr_offset] = sym[i] / d_ref_sym[i - carr_offset];
        }
    }
    */
    /*
    if (d_interpolate) {
        for (int i = d_first_active_carrier + 1; i < d_last_active_carrier; i += 2) {
            taps[i] = taps[i - 1];
        }
        taps[d_last_active_carrier] = taps[d_last_active_carrier - 1];
    }
    */
    if (d_eq_noise_red_len) {
        // TODO
        // 1) IFFT
        // 2) Set all elements > d_eq_noise_red_len to zero
        // 3) FFT
    }
}


// Operates on a per-frame basis
int ofdm_chanest_vcvc_impl::general_work(int noutput_items,
                                         gr_vector_int& ninput_items,
                                         gr_vector_const_void_star& input_items,
                                         gr_vector_void_star& output_items)
{
    const gr_complex* in1 = (const gr_complex*)input_items[0];
    const gr_complex* in2 = (const gr_complex*)input_items[1];
    gr_complex* out1 = (gr_complex*)output_items[0];
    gr_complex* out2 = (gr_complex*)output_items[1];
    const int framesize = 3 + d_n_data_syms;

    // Channel info estimation
    int carr_offset1 = get_carr_offset(in1, in1 + d_fft_len, d_corr_v1, d_max_pos_carr_offset1, d_max_neg_carr_offset1);
    int carr_offset2 = get_carr_offset(in2, in2 + d_fft_len, d_corr_v2, d_max_pos_carr_offset2, d_max_neg_carr_offset2);
    std::vector<std::vector<std::vector<gr_complex>>> chan_taps1 = std::vector<std::vector<std::vector<gr_complex>>>(d_fft_len, std::vector<std::vector<gr_complex>>(2, std::vector<gr_complex>({gr_complex(0, 0), gr_complex(0, 0)})));
    std::vector<std::vector<std::vector<gr_complex>>> chan_taps2 = std::vector<std::vector<std::vector<gr_complex>>>(d_fft_len, std::vector<std::vector<gr_complex>>(2, std::vector<gr_complex>({gr_complex(0, 0), gr_complex(0, 0)})));
    get_chan_taps(in1 + d_fft_len, in1 + 2 * d_fft_len, in2 + d_fft_len, in2 + 2 * d_fft_len, carr_offset1, carr_offset2, chan_taps1, chan_taps2);
    add_item_tag(0,
                 nitems_written(0),
                 pmt::string_to_symbol("ofdm_sync_carr_offset"),
                 pmt::from_long(carr_offset1));
    add_item_tag(1,
                 nitems_written(1),
                 pmt::string_to_symbol("ofdm_sync_carr_offset"),
                 pmt::from_long(carr_offset2));
    /*
    add_item_tag(0,
                 nitems_written(0),
                 pmt::string_to_symbol("ofdm_sync_chan_taps"),
                 pmt::init_c32vector(d_fft_len, chan_taps));
    */
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            char c[10];
            char s[100] = "ofdm_sync_chan_taps";
            int place = i * 2 + j;
            sprintf(c, "%d", place);
            strcat(s, c);
 
            std::vector<gr_complex> temp1(d_fft_len, gr_complex(0, 0));
            std::vector<gr_complex> temp2(d_fft_len, gr_complex(0, 0));
            for (int k = 0; k < d_fft_len; k++)
            {
                temp1[k] = chan_taps1[k][i][j];
                temp2[k] = chan_taps2[k][i][j];
            }
            add_item_tag(0,
                         nitems_written(0),
                         pmt::string_to_symbol(s),
                         pmt::init_c32vector(d_fft_len, temp1));
            add_item_tag(1,
                         nitems_written(1),
                         pmt::string_to_symbol(s),
                         pmt::init_c32vector(d_fft_len, temp2));
        }
    }

    // Copy data symbols
    /*
    if (output_items.size() == 2) {
        gr_complex* out_chantaps = ((gr_complex*)output_items[1]);
        memcpy((void*)out_chantaps, (void*)&chan_taps[0], sizeof(gr_complex) * d_fft_len);
        produce(1, 1);
    }
    */
    memcpy((void*)out1,
           (void*)&in1[3 * d_fft_len],
           sizeof(gr_complex) * d_fft_len * d_n_data_syms);
    memcpy((void*)out2,
           (void*)&in2[3 * d_fft_len],
           sizeof(gr_complex) * d_fft_len * d_n_data_syms);

    // Propagate tags
    std::vector<gr::tag_t> tags;
    get_tags_in_range(tags, 0, nitems_read(0), nitems_read(0) + framesize);
    for (unsigned t = 0; t < tags.size(); t++) {
        int offset = tags[t].offset - nitems_read(0);
        if (offset < 3) {
            offset = 0;
        } else {
            offset -= 3;
        }
        tags[t].offset = offset + nitems_written(0);
        add_item_tag(0, tags[t]);
    }
    get_tags_in_range(tags, 1, nitems_read(1), nitems_read(1) + framesize);
    for (unsigned t = 0; t < tags.size(); t++) {
        int offset = tags[t].offset - nitems_read(1);
        if (offset < 3) {
            offset = 0;
        } else {
            offset -= 3;
        }
        tags[t].offset = offset + nitems_written(1);
        add_item_tag(1, tags[t]);
    }

    produce(0, d_n_data_syms);
    produce(1, d_n_data_syms);
    consume_each(framesize);
    return WORK_CALLED_PRODUCE;
}

//#########################################################################
//Auxiliary functions
//#########################################################################

//########
//[function function]: This function is used to 
//transfer vector type complex vector vector<gr_complex>
//into gr_complex*.
//[function name]: calc_vector_vector
//[parameter]: $name: vector; %type: vector<complex>
//             $name: len; $type: int
//[output]: $name: output; $type: complex*
gr_complex* calc_vector_vector(std::vector<gr_complex> vector, int len)
{
  gr_complex* output = (gr_complex*)malloc(sizeof(gr_complex) * len);
  for (int i = 0; i < len; i++)
  {
    output[i] = vector[i];
  }
  return output;
}



//########
//[function function]: This function is used to 
//transfer vector type complex matrix vector<vector<gr_complex>>
//into gr_complex**.
//[function name]: calc_matrix_matrix
//[parameter]: $name: matrix; %type: vector<vector<complex>>
//             $name: M; $type: int
//             $name: N; $type: int
//[output]: $name: output; $type: complex**
gr_complex** calc_matrix_matrix(std::vector<std::vector<gr_complex>> matrix, int M, int N)
{
  gr_complex** output = (gr_complex**)malloc(sizeof(gr_complex*) * M);
  for (int i = 0; i < M; i++)
  {
    output[i] = (gr_complex*)malloc(sizeof(gr_complex) * N);
  }
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < N; j++)
    {
      output[i][j] = matrix[i][j];
    }
  }
  return output;
}



//########
//[function function]: This function is used to 
//transfer gr_complex**
//into vector type complex matrix vector<vector<gr_complex>>.
//[function name]: calc_matrix_vector
//[parameter]: $name: matrix; %type: gr_complex**
//             $name: M; $type: int
//             $name: N; $type: int
//[output]: $name: output; $type: vector<vector<gr_complex>>
std::vector<std::vector<gr_complex>> calc_matrix_vector(gr_complex** matrix, int M, int N)
{
  std::vector<std::vector<gr_complex>> output(M, std::vector<gr_complex>(N));
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < N; j++)
    {
      output[i][j] = matrix[i][j];
    }
  }
  return output;
}


//########
//[function function]: This function is used to
//concatenate two vector<complex> into one matrix**
//[function name]: calc_vector_concat
//[parameter]: $name: in1; %type: vector<complex>
//             $name: in2; $type: vector<complex>
//             $name: len; $type: int
//[output]: $name: output; $type: complex**
gr_complex** calc_vector_concat(std::vector<gr_complex> in1, std::vector<gr_complex> in2, int len)
{
  gr_complex** output = (gr_complex**)malloc(sizeof(gr_complex*) * 2);
  for (int i = 0; i < 2; i++)
  {
    output[i] = (gr_complex*)malloc(sizeof(gr_complex) * len);
  }
  for (int i = 0; i < len; i++)
  {
    output[0][i] = in1[i];
    output[1][i] = in2[i];
  }
  return output;
}



gr_complex** calc_vector_concat(gr_complex* in1, gr_complex* in2, int len)
{
  gr_complex** output = (gr_complex**)malloc(sizeof(gr_complex*) * 2);
  for (int i = 0; i < 2; i++)
  {
    output[i] = (gr_complex*)malloc(sizeof(gr_complex) * len);
  }
  for (int i = 0; i < len; i++)
  {
    output[0][i] = in1[i];
    output[1][i] = in2[i];
  }
  return output;
}


//########
//[function function]: This function is used to 
//get the A^H of A, used for Hermite.
//[function name]: calc_matrix_H
//[parameter]: $name: matrix_input; %type: complex**
//             $name: M; $type: int, standing for row
//             $name: N; $type: int, standing for col
//[output]: $name: output; $type: complex**
gr_complex** calc_matrix_H(gr_complex** matrix_input, int M, int N)
{
  gr_complex** output = (gr_complex**)malloc(sizeof(gr_complex*) * N);
  for (int i = 0; i < N; i++)
  {
    output[i] = (gr_complex*)malloc(sizeof(gr_complex) * M);
  }
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < N; j++)
    {
      output[j][i].real(matrix_input[i][j].real());
      output[j][i].imag(-1*matrix_input[i][j].imag());
    }
  }
  return output;
}



//########
//[function function]: transform the vector-type matrix
//into true matrix.
//[function name]: calc_matrix_2d
//[parameter]: $name: A; %type: complex**
//             $name: M; $type: int, standing for row
//             $name: N; $type: int, standing for col
//[output]: $name: output; $type: complex**
gr_complex** calc_matrix_2d(gr_complex* A, int M, int N)
{
  gr_complex** output = (gr_complex**)malloc(sizeof(gr_complex*) * M);
  for (int i = 0; i < M; i++) 
  {
    output[i] = (gr_complex*)malloc(sizeof(gr_complex) * N);
  }
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < N; j++)
    {
      output[i][j] = A[i * N + j];
    }
  }
  //free(A);
  return output;
}



//########
//[function function]: This function is used to 
//calculate A*B
//[function name]: calc_matrix_mult
//[parameter]: $name: A; $type: complex**
//             $name: B; $type: complex**
//             $name: A_row; $type: int, standing for row A
//             $name: B_row; $type: int, standing for row B
//             $name: B_col; $type: int, standing for col B
//[output]: $name: output; $type: complex**
gr_complex** calc_matrix_mult(gr_complex** A, gr_complex** B, int A_row, int B_row, int B_col)
{
  gr_complex** output = (gr_complex**)malloc(sizeof(gr_complex*) * A_row);
  for (int i = 0; i < A_row; i++)
  {
    output[i] = (gr_complex*)malloc(sizeof(gr_complex) * B_col);
  }
  for (int i = 0; i < A_row; i++)
  {
    for (int j = 0; j < B_col; j++)
    {
      for (int k = 0; k < B_row; k++)
      {
        output[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return output;
}



//########
//[function function]: This function is used to 
//calculate A^-1. This is the sub-function of
//calc_matrix_inv().
//[function name]: matrix_inv_LU
//[parameter]: $name: input; $type: float**
//             $name: N; $type: int, size of square matrix A
//[output]: $name: output; $type: float**
float** matrix_inv_LU(float** input, int N)
{
  float** L = (float**)malloc(sizeof(float*) * N);
  float** U = (float**)malloc(sizeof(float*) * N);
  float** L_inv = (float**)malloc(sizeof(float*) * N);
  float** U_inv = (float**)malloc(sizeof(float*) * N);
  float** output = (float**)malloc(sizeof(float*) * N);
  for (int i = 0; i < N; i++)
  {
    L[i] = (float*)malloc(sizeof(float) * N);
    U[i] = (float*)malloc(sizeof(float) * N);
    L_inv[i] = (float*)malloc(sizeof(float) * N);
    U_inv[i] = (float*)malloc(sizeof(float) * N);
    output[i] = (float*)malloc(sizeof(float) * N);
  }
  float s;
  
  for (int i = 0; i < N; i++)
  {
    L[i][i] = 1;
  }

  for (int i = 0; i < N; i++)
  {
    for(int j = i; j < N; j++)
    {
      s = 0; 
      for (int k = 0; k < i; k++)
      {
        s += L[i][k] * U[k][j];
      }
      U[i][j] =  input[i][j] - s;
    }

    for (int j = i + 1; j < N; j++)
    {
      s = 0; 
      for (int k = 0; k < i; k++)
      {
        s += L[j][k] * U[k][i];
      }
      L[j][i] = (input[j][i] - s) / U[i][i];
    }
  }

  for (int i = 0; i < N; i++)
  {
    L_inv[i][i] = 1;
  }
  
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < i; j++)
    {
      s = 0; 
      for (int k = 0; k < i; k++)
      {
        s += L[i][k] * L_inv[k][j];
      }
      L_inv[i][j] = -s;
    }
  }

  for (int i = 0; i < N; i++)
  {
    U_inv[i][i] = 1 / U[i][i];
  }

  for (int i = 0; i < N; i++)
  {
    for (int j = i - 1; j >= 0; j--)
    {
      s = 0;
      for (int k = j + 1; k <= i; k++)
      {
        s += U[j][k] * U_inv[k][i];
      }
      U_inv[j][i] = -s / U[j][j];
    }
  }

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      for (int k = 0; k < N; k++)
      { 
        output[i][j] += U_inv[i][k] * L_inv[k][j];
      }
    }
  }
  /*
  free(L);
  free(U);
  free(L_inv);
  free(U_inv);
  */
  return output;
}



//########
//[function function]: This function is used to 
//calculate A*B (float)
//[function name]: calc_float_matrix_mult
//[parameter]: $name: A; $type: float**
//             $name: B; $type: float**
//             $name: A_row; $type: int, standing for row A
//             $name: B_row; $type: int, standing for row B
//             $name: B_col; $type: int, standing for col B
//[output]: $name: output; $type: float**
float** calc_float_matrix_mult(float** A, float** B, int A_row, int B_row, int B_col)
{
  float** output = (float**)malloc(sizeof(float*) * A_row);
  for (int i = 0; i < A_row; i++)
  {
    output[i] = (float*)malloc(sizeof(float*) * B_col);
  }
  for (int i = 0; i < A_row; i++)
  {
    for (int j = 0; j < B_col; j++)
    {
      for (int k = 0; k < B_row; k++)
      {
        output[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return output;
}



//########
//[function function]: This function is used to 
//calculate A+B (float)
//[function name]: calc_float_matrix_plus
//[parameter]: $name: A; $type: float**
//             $name: B; $type: float**
//             $name: M; $type: int
//             $name: N; $type: int
//[output]: $name: output; $type: float**
float** calc_float_matrix_plus(float** A, float** B, int M, int N)
{
  float** output = (float**)malloc(sizeof(float*) * M);
  for (int i = 0; i < N; i++)
  {
    output[i] = (float*)malloc(sizeof(float) * N);
  }
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < N; j++)
    {
      output[i][j] = A[i][j] + B[i][j];
    }
  }
  return output;
}



//########
//[function function]: This function is used to 
//calculate A-B (float)
//[function name]: calc_float_matrix_minus
//[parameter]: $name: A; $type: float**
//             $name: B; $type: float**
//             $name: M; $type: int
//             $name: N; $type: int
//[output]: $name: output; $type: float**
float** calc_float_matrix_minus(float** A, float** B, int M, int N)
{
  float** output = (float**)malloc(sizeof(float*) * M);
  for (int i = 0; i < N; i++)
  {
    output[i] = (float*)malloc(sizeof(float) * N);
  }
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < N; j++)
    {
      output[i][j] = A[i][j] - B[i][j];
    }
  }
  return output;
}



//########
//[function function]: This function is used to 
//calculate A^-1. 
//[function name]: calc_matrix_inv
//[parameter]: $name: input; $type: complex**
//             $name: N; $type: int, size of square matrix A
//[output]: $name: output; $type: float**
gr_complex** calc_matrix_inv(gr_complex** A, int N)
{
  gr_complex** output = (gr_complex**)malloc(sizeof(gr_complex*) * N);
  float** a_real = (float**)malloc(sizeof(float*) * N);
  float** a_imag = (float**)malloc(sizeof(float*) * N);
  float** Ainv = NULL;
  float** Binv = NULL;
  float** BAinv = NULL;
  float** BAinvB = NULL;
  float** A_P_BAinvB = NULL;
  float** A_P_BAinvB_inv = NULL;
  float** AinvB = NULL;
  float** AinvB_A_P_BAinvB_inv = NULL;

  for (int i = 0; i < N; i++)
  {
    output[i] = (gr_complex*)malloc(sizeof(gr_complex) * N);
    a_real[i] = (float*)malloc(sizeof(float) * N);
    a_imag[i] = (float*)malloc(sizeof(float) * N);
  }

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      a_real[i][j] = A[i][j].real();
      a_imag[i][j] = A[i][j].imag();
    }
  }

  Ainv = matrix_inv_LU(a_real, N);
  Binv = matrix_inv_LU(a_imag, N);
  BAinv = calc_float_matrix_mult(a_imag, Ainv, N, N, N);
  BAinvB = calc_float_matrix_mult(BAinv, a_imag, N, N, N);
  A_P_BAinvB = calc_float_matrix_plus(a_real, BAinvB, N, N);
  A_P_BAinvB_inv = matrix_inv_LU(A_P_BAinvB, N);
  AinvB = calc_float_matrix_mult(Ainv, a_imag, N, N, N);
  AinvB_A_P_BAinvB_inv = calc_float_matrix_mult(AinvB, A_P_BAinvB_inv, N, N, N);

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      AinvB_A_P_BAinvB_inv[i][j] = -AinvB_A_P_BAinvB_inv[i][j];
    }
  }

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      output[i][j].real(A_P_BAinvB_inv[i][j]);
      output[i][j].imag(AinvB_A_P_BAinvB_inv[i][j]);
    }
  }
  /*
  free(a_real);
  free(a_imag);
  free(Ainv);
  free(Binv);
  free(BAinv);
  free(BAinvB);
  free(A_P_BAinvB);
  free(A_P_BAinvB_inv);
  free(AinvB);
  free(AinvB_A_P_BAinvB_inv);
  */
  return output;
}


} /* namespace digital */
} /* namespace gr */
