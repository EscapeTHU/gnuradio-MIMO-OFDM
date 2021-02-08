/* -*- c++ -*- */
/* Copyright 2012,2018 Free Software Foundation, Inc.
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

#include "ofdm_frame_equalizer_vcvc_impl.h"
#include <gnuradio/expj.h>
#include <gnuradio/io_signature.h>
#include <gnuradio/math.h>

static const pmt::pmt_t CARR_OFFSET_KEY = pmt::mp("ofdm_sync_carr_offset");
static const pmt::pmt_t CHAN_TAPS_KEY1 = pmt::mp("ofdm_sync_chan_taps1");
static const pmt::pmt_t CHAN_TAPS_KEY2 = pmt::mp("ofdm_sync_chan_taps2");
static const pmt::pmt_t CHAN_TAPS_KEY3 = pmt::mp("ofdm_sync_chan_taps3");
static const pmt::pmt_t CHAN_TAPS_KEY4 = pmt::mp("ofdm_sync_chan_taps4");

namespace gr {
namespace digital {

ofdm_frame_equalizer_vcvc::sptr
ofdm_frame_equalizer_vcvc::make(ofdm_equalizer_base::sptr equalizer,
                                int cp_len,
                                const std::string& tsb_key,
                                bool propagate_channel_state,
                                int fixed_frame_len)
{
    return gnuradio::get_initial_sptr(new ofdm_frame_equalizer_vcvc_impl(
        equalizer, cp_len, tsb_key, propagate_channel_state, fixed_frame_len));
}

ofdm_frame_equalizer_vcvc_impl::ofdm_frame_equalizer_vcvc_impl(
    ofdm_equalizer_base::sptr equalizer,
    int cp_len,
    const std::string& tsb_key,
    bool propagate_channel_state,
    int fixed_frame_len)
    : tagged_stream_block(
          "ofdm_frame_equalizer_vcvc",
          io_signature::make2(2, 2, sizeof(gr_complex) * equalizer->fft_len(), sizeof(gr_complex) * equalizer->fft_len()),
          io_signature::make2(2, 2, sizeof(gr_complex) * equalizer->fft_len(), sizeof(gr_complex) * equalizer->fft_len()),
          tsb_key),
      d_fft_len(equalizer->fft_len()),
      d_cp_len(cp_len),
      d_eq(equalizer),
      d_propagate_channel_state(propagate_channel_state),
      d_fixed_frame_len(fixed_frame_len),
      d_channel_state1(std::vector<std::vector<std::vector<gr_complex>>>(equalizer->fft_len(), std::vector<std::vector<gr_complex>>({{gr_complex(0, 0), gr_complex(0, 0)}, {gr_complex(0, 0), gr_complex(0, 0)}}))),
      d_channel_state2(std::vector<std::vector<std::vector<gr_complex>>>(equalizer->fft_len(), std::vector<std::vector<gr_complex>>({{gr_complex(0, 0), gr_complex(0, 0)}, {gr_complex(0, 0), gr_complex(0, 0)}})))
{
    if (tsb_key.empty() && fixed_frame_len == 0) {
        throw std::invalid_argument("Either specify a TSB tag or a fixed frame length!");
    }
    if (d_fixed_frame_len < 0) {
        throw std::invalid_argument("Invalid frame length!");
    }
    if (d_fixed_frame_len) {
        set_output_multiple(d_fixed_frame_len);
    }
    set_relative_rate(1, 1);
    // Really, we have TPP_ONE_TO_ONE, but the channel state is not propagated
    set_tag_propagation_policy(TPP_DONT);
}

ofdm_frame_equalizer_vcvc_impl::~ofdm_frame_equalizer_vcvc_impl() {}

void ofdm_frame_equalizer_vcvc_impl::parse_length_tags(
    const std::vector<std::vector<tag_t>>& tags, gr_vector_int& n_input_items_reqd)
{
    if (d_fixed_frame_len) {
        n_input_items_reqd[0] = d_fixed_frame_len;
        n_input_items_reqd[1] = d_fixed_frame_len;
    } else {
        for (unsigned k = 0; k < tags[0].size(); k++) {
            if (tags[0][k].key == pmt::string_to_symbol(d_length_tag_key_str)) {
                n_input_items_reqd[0] = pmt::to_long(tags[0][k].value);
                n_input_items_reqd[1] = pmt::to_long(tags[0][k].value);
            }
        }
    }
}


int ofdm_frame_equalizer_vcvc_impl::work(int noutput_items,
                                         gr_vector_int& ninput_items,
                                         gr_vector_const_void_star& input_items,
                                         gr_vector_void_star& output_items)
{
    const gr_complex* in1 = (const gr_complex*)input_items[0];
    const gr_complex* in2 = (const gr_complex*)input_items[1];
    gr_complex* out1 = (gr_complex*)output_items[0];
    gr_complex* out2 = (gr_complex*)output_items[1];
    int carrier_offset1 = 0;
    int carrier_offset2 = 0;
    int frame_len = 0;
    if (d_fixed_frame_len) {
        frame_len = d_fixed_frame_len;
    } else {
        frame_len = ninput_items[0];
    }

    std::vector<tag_t> tags1;
    get_tags_in_window(tags1, 0, 0, frame_len);
    for (unsigned i = 0; i < tags1.size(); i++) {
        if (pmt::symbol_to_string(tags1[i].key) == "ofdm_sync_chan_taps0") {
            std::vector<gr_complex> temp = pmt::c32vector_elements(tags1[i].value);
            for (int i = 0; i < d_fft_len; i++)
            {
                d_channel_state1[i][0][0] = temp[i];
            }
        }
        if (pmt::symbol_to_string(tags1[i].key) == "ofdm_sync_chan_taps1") {
            std::vector<gr_complex> temp = pmt::c32vector_elements(tags1[i].value);
            for (int i = 0; i < d_fft_len; i++)
            {
                d_channel_state1[i][0][1] = temp[i];
            }
        }
        if (pmt::symbol_to_string(tags1[i].key) == "ofdm_sync_chan_taps2") {
            std::vector<gr_complex> temp = pmt::c32vector_elements(tags1[i].value);
            for (int i = 0; i < d_fft_len; i++)
            {
                d_channel_state1[i][1][0] = temp[i];
            }
        }
        if (pmt::symbol_to_string(tags1[i].key) == "ofdm_sync_chan_taps3") {
            std::vector<gr_complex> temp = pmt::c32vector_elements(tags1[i].value);
            for (int i = 0; i < d_fft_len; i++)
            {
                d_channel_state1[i][1][1] = temp[i];
            }
        }
        if (pmt::symbol_to_string(tags1[i].key) == "ofdm_sync_carr_offset") {
            carrier_offset1 = pmt::to_long(tags1[i].value);
        }
    }
    std::vector<tag_t> tags2;
    get_tags_in_window(tags2, 1, 0, frame_len);
    for (unsigned i = 0; i < tags2.size(); i++) {
        if (pmt::symbol_to_string(tags2[i].key) == "ofdm_sync_chan_taps0") {
            std::vector<gr_complex> temp = pmt::c32vector_elements(tags2[i].value);
            for (int i = 0; i < d_fft_len; i++)
            {
                d_channel_state2[i][0][0] = temp[i];
            }
        }
        if (pmt::symbol_to_string(tags2[i].key) == "ofdm_sync_chan_taps1") {
            std::vector<gr_complex> temp = pmt::c32vector_elements(tags2[i].value);
            for (int i = 0; i < d_fft_len; i++)
            {
                d_channel_state2[i][0][1] = temp[i];
            }
        }
        if (pmt::symbol_to_string(tags2[i].key) == "ofdm_sync_chan_taps2") {
            std::vector<gr_complex> temp = pmt::c32vector_elements(tags2[i].value);
            for (int i = 0; i < d_fft_len; i++)
            {
                d_channel_state2[i][1][0] = temp[i];
            }
        }
        if (pmt::symbol_to_string(tags2[i].key) == "ofdm_sync_chan_taps3") {
            std::vector<gr_complex> temp = pmt::c32vector_elements(tags2[i].value);
            for (int i = 0; i < d_fft_len; i++)
            {
                d_channel_state2[i][1][1] = temp[i];
            }
        }
        if (pmt::symbol_to_string(tags2[i].key) == "ofdm_sync_carr_offset") {
            carrier_offset2 = pmt::to_long(tags2[i].value);
        }
    }

    // Copy the frame and the channel state vector such that the symbols are shifted to
    // the correct position
    if (carrier_offset1 < 0) {
        memset((void*)out1, 0x00, sizeof(gr_complex) * (-carrier_offset1));
        memcpy((void*)&out1[-carrier_offset1],
               (void*)in1,
               sizeof(gr_complex) * (d_fft_len * frame_len + carrier_offset1));
    } else {
        memset((void*)(out1 + d_fft_len * frame_len - carrier_offset1),
               0x00,
               sizeof(gr_complex) * carrier_offset1);
        memcpy((void*)out1,
               (void*)(in1 + carrier_offset1),
               sizeof(gr_complex) * (d_fft_len * frame_len - carrier_offset1));
    }
    if (carrier_offset2 < 0) {
        memset((void*)out2, 0x00, sizeof(gr_complex) * (-carrier_offset2));
        memcpy((void*)&out2[-carrier_offset2],
               (void*)in2,
               sizeof(gr_complex) * (d_fft_len * frame_len + carrier_offset2));
    } else {
        memset((void*)(out2 + d_fft_len * frame_len - carrier_offset2),
               0x00,
               sizeof(gr_complex) * carrier_offset2);
        memcpy((void*)out2,
               (void*)(in2 + carrier_offset2),
               sizeof(gr_complex) * (d_fft_len * frame_len - carrier_offset2));
    }

    // Correct the frequency shift on the symbols
    gr_complex phase_correction1;
    for (int i = 0; i < frame_len; i++) {
        phase_correction1 =
            gr_expj(-(2.0 * GR_M_PI) * carrier_offset1 * d_cp_len / d_fft_len * (i + 1));
        for (int k = 0; k < d_fft_len; k++) {
            out1[i * d_fft_len + k] *= phase_correction1;
        }
    }
    gr_complex phase_correction2;
    for (int i = 0; i < frame_len; i++) {
        phase_correction2 =
            gr_expj(-(2.0 * GR_M_PI) * carrier_offset2 * d_cp_len / d_fft_len * (i + 1));
        for (int k = 0; k < d_fft_len; k++) {
            out2[i * d_fft_len + k] *= phase_correction2;
        }
    }

    // Do the equalizing
    d_eq->reset();
    d_eq->mimo_equalize(out1, out2, frame_len, d_channel_state1, d_channel_state2);

    // Update the channel state regarding the frequency offset
    phase_correction1 =
        gr_expj((2.0 * GR_M_PI) * carrier_offset1 * d_cp_len / d_fft_len * frame_len);
    for (int k = 0; k < d_fft_len; k++) {
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                d_channel_state1[k][i][j] *= phase_correction1;
            }
        }
    }
    phase_correction2 =
        gr_expj((2.0 * GR_M_PI) * carrier_offset2 * d_cp_len / d_fft_len * frame_len);
    for (int k = 0; k < d_fft_len; k++) {
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                d_channel_state2[k][i][j] *= phase_correction2;
            }
        }
    }

    // Propagate tags (except for the channel state and the TSB tag)
    std::vector<tag_t> tags;
    get_tags_in_window(tags, 0, 0, frame_len);
    for (size_t i = 0; i < tags.size(); i++) {
        if (tags[i].key != CHAN_TAPS_KEY1 &&
            tags[i].key != pmt::mp(d_length_tag_key_str) && tags[i].key != CHAN_TAPS_KEY2 && tags[i].key != CHAN_TAPS_KEY3 && tags[i].key != CHAN_TAPS_KEY4) 
        {
            add_item_tag(0, tags[i]);
        }
    }
    get_tags_in_window(tags, 1, 0, frame_len);
    for (size_t i = 0; i < tags.size(); i++) {
        if (tags[i].key != CHAN_TAPS_KEY1 &&
            tags[i].key != pmt::mp(d_length_tag_key_str) && tags[i].key != CHAN_TAPS_KEY2 && tags[i].key != CHAN_TAPS_KEY3 && tags[i].key != CHAN_TAPS_KEY4) 
        {
            add_item_tag(1, tags[i]);
        }
    }

    // Housekeeping
    /*
    if (d_propagate_channel_state) {
        add_item_tag(0,
                     nitems_written(0),
                     pmt::string_to_symbol("ofdm_sync_chan_taps"),
                     pmt::init_c32vector(d_fft_len, d_channel_state));
    }
    */
    if (d_propagate_channel_state) {
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
                    temp1[k] = d_channel_state1[k][i][j];
                    temp2[k] = d_channel_state2[k][i][j];
                }
                add_item_tag(0, nitems_written(0), pmt::string_to_symbol(s), pmt::init_c32vector(d_fft_len, temp1));
                add_item_tag(1, nitems_written(1), pmt::string_to_symbol(s), pmt::init_c32vector(d_fft_len, temp2));
            }
        }
    }

    if (d_fixed_frame_len && d_length_tag_key_str.empty()) {
        consume_each(frame_len);
    }

    return frame_len;
}

} /* namespace digital */
} /* namespace gr */
