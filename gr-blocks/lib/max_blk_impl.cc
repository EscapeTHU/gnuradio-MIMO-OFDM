/* -*- c++ -*- */
/*
 * Copyright 2014,2015,2018 Free Software Foundation, Inc.
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

#include "max_blk_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace blocks {

template <class T>
typename max_blk<T>::sptr max_blk<T>::make(size_t vlen, size_t vlen_out)
{
    return gnuradio::get_initial_sptr(new max_blk_impl<T>(vlen, vlen_out));
}

template <class T>
max_blk_impl<T>::max_blk_impl(size_t vlen, size_t vlen_out)
    : sync_block("max_blk",
                 io_signature::make(1, -1, vlen * sizeof(T)),
                 io_signature::make(1, 1, vlen_out * sizeof(T))),
      d_vlen(vlen),
      d_vlen_out(vlen_out)
{
    assert((vlen_out == vlen) || (vlen_out == 1));
}

template <class T>
max_blk_impl<T>::~max_blk_impl()
{
}

template <class T>
int max_blk_impl<T>::work(int noutput_items,
                          gr_vector_const_void_star& input_items,
                          gr_vector_void_star& output_items)
{
    T* optr = (T*)output_items[0];

    int ninputs = input_items.size();

    if (d_vlen_out == 1)
        for (int i = 0; i < noutput_items; i++) {
            T max = ((T*)input_items[0])[i * d_vlen];

            for (int j = 0; j < (int)d_vlen; j++) {
                for (int k = 0; k < ninputs; k++) {
                    if (((T*)input_items[k])[i * d_vlen + j] > max) {
                        max = ((T*)input_items[k])[i * d_vlen + j];
                    }
                }
            }

            *optr++ = (T)max;
        }

    else // vector mode output
        for (size_t i = 0; i < (size_t)noutput_items * d_vlen_out; i++) {
            T max = ((T*)input_items[0])[i];

            for (int k = 1; k < ninputs; k++) {
                if (((T*)input_items[k])[i] > max) {
                    max = ((T*)input_items[k])[i];
                }
            }

            *optr++ = (T)max;
        }

    return noutput_items;
}

template class max_blk<std::int16_t>;
template class max_blk<std::int32_t>;
template class max_blk<float>;
} /* namespace blocks */
} /* namespace gr */
