/* -*- c++ -*- */
/* Copyright 2012 Free Software Foundation, Inc.
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

#include <gnuradio/digital/ofdm_equalizer_simpledfe.h>
#include "ofdm_chanest_vcvc_impl.h"

namespace gr {
namespace digital {

ofdm_equalizer_simpledfe::sptr
ofdm_equalizer_simpledfe::make(int fft_len,
                               const gr::digital::constellation_sptr& constellation,
                               const std::vector<std::vector<int>>& occupied_carriers,
                               const std::vector<std::vector<int>>& pilot_carriers,
                               const std::vector<std::vector<gr_complex>>& pilot_symbols,
                               int symbols_skipped,
                               float alpha,
                               bool input_is_shifted)
{
    return ofdm_equalizer_simpledfe::sptr(new ofdm_equalizer_simpledfe(fft_len,
                                                                       constellation,
                                                                       occupied_carriers,
                                                                       pilot_carriers,
                                                                       pilot_symbols,
                                                                       symbols_skipped,
                                                                       alpha,
                                                                       input_is_shifted));
}

ofdm_equalizer_simpledfe::ofdm_equalizer_simpledfe(
    int fft_len,
    const gr::digital::constellation_sptr& constellation,
    const std::vector<std::vector<int>>& occupied_carriers,
    const std::vector<std::vector<int>>& pilot_carriers,
    const std::vector<std::vector<gr_complex>>& pilot_symbols,
    int symbols_skipped,
    float alpha,
    bool input_is_shifted)
    : ofdm_equalizer_1d_pilots(fft_len,
                               occupied_carriers,
                               pilot_carriers,
                               pilot_symbols,
                               symbols_skipped,
                               input_is_shifted),
      d_constellation(constellation),
      d_alpha(alpha)
{
}


ofdm_equalizer_simpledfe::~ofdm_equalizer_simpledfe() {}


void ofdm_equalizer_simpledfe::equalize(gr_complex* frame,
                                        int n_sym,
                                        const std::vector<gr_complex>& initial_taps,
                                        const std::vector<tag_t>& tags)
{
    if (!initial_taps.empty()) {
        d_channel_state = initial_taps;
    }
    gr_complex sym_eq, sym_est;

    for (int i = 0; i < n_sym; i++) {
        for (int k = 0; k < d_fft_len; k++) {
            if (!d_occupied_carriers[k]) {
                continue;
            }
            if (!d_pilot_carriers.empty() && d_pilot_carriers[d_pilot_carr_set][k]) {
                d_channel_state[k] = d_alpha * d_channel_state[k] +
                                     (1 - d_alpha) * frame[i * d_fft_len + k] /
                                         d_pilot_symbols[d_pilot_carr_set][k];
                frame[i * d_fft_len + k] = d_pilot_symbols[d_pilot_carr_set][k];
            } else {
                sym_eq = frame[i * d_fft_len + k] / d_channel_state[k];
                // The `map_to_points` function will treat `sym_est` as an array
                // pointer.  This call is "safe" because `map_to_points` is limited
                // by the dimensionality of the constellation. This class calls the
                // `constellation` class default constructor, which initializes the
                // dimensionality value to `1`. Thus, Only the single `gr_complex`
                // value will be dereferenced.
                d_constellation->map_to_points(d_constellation->decision_maker(&sym_eq),
                                               &sym_est);
                d_channel_state[k] = d_alpha * d_channel_state[k] +
                                     (1 - d_alpha) * frame[i * d_fft_len + k] / sym_est;
                frame[i * d_fft_len + k] = sym_est;
            }
        }
        if (!d_pilot_carriers.empty()) {
            d_pilot_carr_set = (d_pilot_carr_set + 1) % d_pilot_carriers.size();
        }
    }
}

void ofdm_equalizer_simpledfe::mimo_equalize(gr_complex* out1,
                                             gr_complex* out2,
                                             int n_sym,
                                             const std::vector<std::vector<std::vector<gr_complex>>>& initial_taps1,
                                             const std::vector<std::vector<std::vector<gr_complex>>>& initial_taps2)
{
    for (int i = 0; i < n_sym; i++)
    {
        for (int j = 0; j < d_fft_len; j++)
        {
            if (!d_occupied_carriers[j])
            {
                continue;
            }
            gr_complex* tempv1 = (gr_complex*)malloc(sizeof(gr_complex) * 2);
            gr_complex* tempv2 = (gr_complex*)malloc(sizeof(gr_complex) * 2);
            tempv1[0] = out1[i * d_fft_len + j];
            tempv1[1] = out2[i * d_fft_len + j];
            tempv2[0] = out1[i * d_fft_len + j];
            tempv2[1] = out2[i * d_fft_len + j];
            gr_complex** tempv1_2d = calc_matrix_2d(tempv1, 2, 1);
            gr_complex** tempv2_2d = calc_matrix_2d(tempv2, 2, 1);
            gr_complex** channel_ptr1 = calc_matrix_matrix(initial_taps1[j], 2, 2);
            gr_complex** channel_ptr2 = calc_matrix_matrix(initial_taps2[j], 2, 2);
            gr_complex** temp_channel1 = calc_matrix_inv(channel_ptr1, 2);
            gr_complex** temp_channel2 = calc_matrix_inv(channel_ptr2, 2);
            gr_complex** temp_out1 = calc_matrix_mult(temp_channel1, tempv1_2d, 2, 2, 1);
            gr_complex** temp_out2 = calc_matrix_mult(temp_channel2, tempv2_2d, 2, 2, 1);
            out1[i * d_fft_len + j] = (temp_out1[0][0] + temp_out2[0][0]) / gr_complex(2, 0);
            out2[i * d_fft_len + j] = (temp_out1[1][0] + temp_out2[1][0]) / gr_complex(2, 0);
            /*
            free(tempv1);
            free(tempv2);
            free(tempv1_2d);
            free(tempv2_2d);
            free(channel_ptr1);
            free(channel_ptr2);
            free(temp_out1);
            free(temp_out2);
            */
        }
    }
}


} /* namespace digital */
} /* namespace gr */
