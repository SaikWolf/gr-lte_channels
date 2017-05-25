/* -*- c++ -*- */
/* 
 * Copyright 2017 Bill Clark.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_LTE_CHANNELS_RAYLEIGH_CHANNEL_IMPL_H
#define INCLUDED_LTE_CHANNELS_RAYLEIGH_CHANNEL_IMPL_H

#include <lte_channels/rayleigh_channel.h>
#include <lte_channels/doppler_base.h>

#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/random.h>

#include <boost/math/special_functions/sinc.hpp>

#include "jakes_doppler.hpp"
#include "cubic_spline_gaussian.hpp"

namespace gr {
  namespace lte_channels {

    typedef std::complex<float> complexf;

    class rayleigh_channel_impl : public rayleigh_channel
    {
     private:
      double d_samp_rate;
      double d_max_dop;
      std::vector<double> d_delays;
      std::vector<double> d_weights;
      std::vector<double> d_taus;
      bool d_norm;
      bool d_lp;
      int    d_seed;
      size_t d_mode;
      size_t d_dop_res;
      size_t d_path_count;
      size_t d_d_tap_count;

      std::vector<Cubic_Spline_Gaussian*>  d_paths;
      gr::filter::kernel::fir_filter_ccc*   d_filter;

      std::vector<float>  d_taps;

      std::vector<complexf> d_analog_channel;

      std::vector< std::vector<complexf> > d_channel;

      void load_channel(int count);
      void setup_paths();
      void gen_sinc_matrix();
      std::vector<double> d_sinc;

      void sinc_multi(std::vector<complexf> &out, std::vector<complexf> &in);

      Doppler_Base* d_doppler;
      

     public:
      rayleigh_channel_impl(double samp_rate, double maximum_doppler, size_t path_count, const std::vector<double> &delays, const std::vector<double> &weights, size_t doppler_spectrum, size_t doppler_taps, bool normalize, bool lead_phase, int seed);
      ~rayleigh_channel_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);

      int get_tap_len(){ return d_d_tap_count; }
    };

  } // namespace lte_channels
} // namespace gr

#endif /* INCLUDED_LTE_CHANNELS_RAYLEIGH_CHANNEL_IMPL_H */

