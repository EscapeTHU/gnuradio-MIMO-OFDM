/* -*- c++ -*- */
/*
 * Copyright 2013 Free Software Foundation, Inc.
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

#ifndef INCLUDED_DIGITAL_OFDM_CHANEST_VCVC_IMPL_H
#define INCLUDED_DIGITAL_OFDM_CHANEST_VCVC_IMPL_H

#include <gnuradio/digital/ofdm_chanest_vcvc.h>

namespace gr {
namespace digital {

class ofdm_chanest_vcvc_impl : public ofdm_chanest_vcvc
{
private:
    int d_fft_len;     //! FFT length
    int d_n_data_syms; //! Number of data symbols following the sync symbol(s)

    //!! The following explanation is written by EscapeTHU (Tengyu Zhang)
    //!! This d_n_sync_syms will NOT represent "Number of sync symbols (1 or 2)"
    //!! But the value will still be limited as 1 or 2 simply as an indicator
    //!! Indicating the operation that "get_carr_offset" and "get_chan_taps" will
    //!! perform.
    int d_n_sync_syms; //! Number of sync symbols (1 or 2)
    //! 0 if no noise reduction is done for the initial channel state estimation.
    //! Otherwise, the maximum length of the channel delay in samples.
    int d_eq_noise_red_len;
    //! Is sync_symbol1 if d_n_sync_syms == 1, otherwise sync_symbol2. Used as a reference
    //! symbol to estimate the channel.
    //!! The following explanation is written by EscapeTHU (Tengyu Zhang)
    //!! 1_1 will be sync_sym2
    //!! 1_2 will be sync_sym3
    //!! 2_1 will be sync_sym4
    //!! 2_2 will be sync_sym5
    std::vector<gr_complex> d_ref_sym1_1;
    std::vector<gr_complex> d_ref_sym1_2;
    std::vector<gr_complex> d_ref_sym2_1;
    std::vector<gr_complex> d_ref_sym2_2;
    //! If d_n_sync_syms == 2 this is used as a differential correlation vector (called
    //! 'v' in [1]).
    std::vector<gr_complex> d_corr_v1;
    std::vector<gr_complex> d_corr_v2;

    //! If d_n_sync_syms == 1 we use this instead of d_corr_v to estimate the coarse freq.
    //! offset
    std::vector<float> d_known_symbol_diffs;
    //! If d_n_sync_syms == 1 we use this instead of d_corr_v to estimate the coarse freq.
    //! offset (temp. variable)
    std::vector<float> d_new_symbol_diffs;

    //! The index of the first carrier with data (index 0 is not DC here, but the lowest
    //! frequency)
    int d_first_active_carrier;
    //! The index of the last carrier with data
    int d_last_active_carrier;
    //! If true, the channel estimation must be interpolated
    bool d_interpolate;

    //! Maximum carrier offset (negative value!)
    int d_max_neg_carr_offset1;
    int d_max_neg_carr_offset2;
    //! Maximum carrier offset (positive value!)
    int d_max_pos_carr_offset1;
    int d_max_pos_carr_offset2;

    //! Calculate the coarse frequency offset in number of carriers
    //!! The following explanation is written by EscapeTHU (Tengyu Zhang)
    //!! This "get_carr_offset" function will be executed by input1 and input2 each;
    //!! AND INPUT1/2 WILL BOTH USE SYNC1 TO ESTIMATE CARRIER OFFSET!!
    int get_carr_offset(const gr_complex* sync_sym1, const gr_complex* sync_sym2, std::vector<gr_complex> d_corr_v, int d_max_pos_carr_offset, int d_max_neg_carr_offset);
    //! Estimate the channel (phase and amplitude offset per carrier)
    //!! The following explanation is written by EscapeTHU (Tengyu Zhang)
    //!! This "get_chan_taps" function will be executed only ONCE, using matrix;
    void get_chan_taps(const gr_complex* sync_sym1_1,
                       const gr_complex* sync_sym1_2,
                       const gr_complex* sync_sym2_1,
                       const gr_complex* sync_sym2_2,
                       int carr_offset1,
                       int carr_offset2,
                       std::vector<std::vector<std::vector<gr_complex>>>& taps1,
                       std::vector<std::vector<std::vector<gr_complex>>>& taps2);

public:
    ofdm_chanest_vcvc_impl(const std::vector<gr_complex>& sync_symbol1,
                           const std::vector<gr_complex>& sync_symbol2,
                           const std::vector<gr_complex>& sync_symbol3,
                           const std::vector<gr_complex>& sync_symbol4,
                           const std::vector<gr_complex>& sync_symbol5,
                           int n_data_symbols,
                           int eq_noise_red_len,
                           int max_carr_offset,
                           bool force_one_sync_symbol);
    ~ofdm_chanest_vcvc_impl();

    void forecast(int noutput_items, gr_vector_int& ninput_items_required);
    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items);
};



gr_complex* calc_vector_vector(std::vector<gr_complex> vector, int len);
gr_complex** calc_matrix_matrix(std::vector<std::vector<gr_complex>> matrix, int M, int N);
std::vector<std::vector<gr_complex>> calc_matrix_vector(gr_complex** matrix, int M, int N);
gr_complex** calc_vector_concat(std::vector<gr_complex> in1, std::vector<gr_complex> in2, int len);
gr_complex** calc_vector_concat(gr_complex* in1, gr_complex* in2, int len);
gr_complex** calc_matrix_H(gr_complex** matrix_input, int M, int N);
gr_complex** calc_matrix_2d(gr_complex* A, int M, int N);
gr_complex** calc_matrix_mult(gr_complex** A, gr_complex** B, int A_row, int B_row, int B_col);
float** matrix_inv_LU(float** input, int N);
float** calc_float_matrix_mult(float** A, float** B, int A_row, int B_row, int B_col);
float** calc_float_matrix_plus(float** A, float** B, int M, int N);
float** calc_float_matrix_minus(float** A, float** B, int M, int N);
gr_complex** calc_matrix_inv(gr_complex** A, int N);


} // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_OFDM_CHANEST_VCVC_IMPL_H */
