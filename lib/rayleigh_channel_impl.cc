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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "rayleigh_channel_impl.h"
#include <iostream>

namespace gr {
  namespace lte_channels {

    rayleigh_channel::sptr
    rayleigh_channel::make(double samp_rate, double maximum_doppler, size_t path_count, const std::vector<double> &delays, const std::vector<double> &weights, size_t doppler_spectrum, size_t doppler_taps, bool normalize, bool lead_phase, int seed)
    {
      return gnuradio::get_initial_sptr
        (new rayleigh_channel_impl(samp_rate, maximum_doppler, path_count, delays, weights, doppler_spectrum, doppler_taps, normalize, lead_phase, seed));
    }

    /*
     * The private constructor
     */
    rayleigh_channel_impl::rayleigh_channel_impl(double samp_rate, double maximum_doppler, size_t path_count, const std::vector<double> &delays, const std::vector<double> &weights, size_t doppler_spectrum, size_t doppler_taps, bool normalize, bool lead_phase, int seed)
      : gr::sync_block("rayleigh_channel",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_mode(doppler_spectrum),
        d_dop_res(doppler_taps),
        d_samp_rate(samp_rate),
        d_max_dop(maximum_doppler),
        d_norm(normalize),
        d_lp(lead_phase),
        d_path_count(path_count),
        d_seed(seed)
    {
      d_delays = std::vector<double>( delays.begin(), delays.end() );
      d_weights = std::vector<double>( weights.begin(), weights.end() );

      double ww = 0.;
      for(size_t idx = 0; idx < d_weights.size(); idx++){
        d_weights[idx] = std::pow(10., d_weights[idx]/10.);
        ww += d_weights[idx];
      }
      if(ww) ww = 1./std::sqrt(ww);
      else ww = 1.;
      for(size_t idx = 0; idx < d_weights.size(); idx++){
        d_weights[idx] *= ww;
      }

      //std::cout << "RC: In constructor, weights("<<d_weights.size()<<"),"
      //          << " delays("<<d_delays.size()<<")"<<std::endl;
      setup_paths();
      gen_sinc_matrix();
      //std::cout << "RC: sinc_mat made\n";

      d_analog_channel = std::vector<complexf>(d_path_count);
      //std::cout << "RC: dac made\n";
      std::vector<complexf> digital_channel(d_d_tap_count,complexf(0.,0.));
      d_filter = new gr::filter::kernel::fir_filter_ccc(1, digital_channel);
      //std::cout << "RC: filt made\n";
      d_channel = std::vector< std::vector<complexf> >(1, digital_channel);
      //std::cout << "RC: chan init made\n";
      //std::cout << "RC: constructed\n";

      set_history(d_d_tap_count);
    }

    /*
     * Our virtual destructor.
     */
    rayleigh_channel_impl::~rayleigh_channel_impl()
    {
      delete d_doppler;
      for(size_t idx = 0; idx < d_path_count; idx++){
        delete d_paths[idx];
      }
      delete d_filter;
    }

    int
    rayleigh_channel_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      load_channel(noutput_items);
      size_t oo(0);
      while(oo < noutput_items){
        d_filter->set_taps(d_channel[oo]);
        out[oo] = d_filter->filter( &in[oo] );
        oo++;
      }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

    void
    rayleigh_channel_impl::setup_paths()
    {
      float delta;
      //std::cout << "RC: Setting up paths of type ";
      switch (d_mode){
        case 0:{
          d_doppler = new Jakes_Doppler(d_dop_res, 0.5, 1.0);//use default constructor
          if(d_max_dop > 0){
            float coherence_time = 0.423*d_samp_rate/d_max_dop;
            delta = 1/coherence_time;
          }
          else{
            delta = 0.;
          }
          //std::cout << "Jakes with delta("<<delta<<")\n";
          break;
        }
        default:{
          d_doppler = new Jakes_Doppler();
          if(d_max_dop > 0){
            float coherence_time = 0.423*d_samp_rate/d_max_dop;
            delta = 1/coherence_time;
          }
          else{
            delta = 0.;
          }
          //std::cout << "Def:Jakes with delta("<<delta<<")\n";
          break;
        }
      }
      if(delta > 1.0) delta = 1.;
      d_doppler->get_taps(d_taps);
      //std::cout << "RC: Doppler Taps acquired, taps = ["<<d_taps[0];
      //for(size_t idx = 1; idx < d_taps.size(); idx++){
      //  std::cout << ", " << d_taps[idx];
      //}
      //std::cout << "];\n";
      d_paths = std::vector<Cubic_Spline_Gaussian*>(d_path_count);
      int lseed;
      for(size_t idx = 0; idx < d_path_count; idx++){
        lseed = (d_seed < 0) ? d_seed : d_seed+idx;
        d_paths[idx] = new Cubic_Spline_Gaussian(d_taps, delta, lseed);
      }
      //std::cout<<"RC: Paths created\n";
    }

    void
    rayleigh_channel_impl::gen_sinc_matrix()
    {
      d_taus = std::vector<double>(d_delays.size());
      double mintau,maxtau;
      for(size_t idx = 1; idx < d_taus.size(); idx++){
        d_taus[idx] = d_delays[idx]*d_samp_rate;
        if(idx == 0){
          mintau = d_taus[idx];
          maxtau = d_taus[idx];
        }
        else{
          mintau = (d_taus[idx] < mintau) ? d_taus[idx] : mintau;
          maxtau = (d_taus[idx] > maxtau) ? d_taus[idx] : maxtau;
        }
      }
      //std::cout << "RC: taus = ["<<d_taus[0];
      //for(size_t idx = 0; idx < d_taus.size(); idx++){
      //  std::cout << ", " << d_taus[idx];
      //}
      //std::cout << "];\n";
      int startat = (mintau-2.718281828459045 < 0.) ? int(std::floor(mintau-2.718281828459045)) : 0;
      int endat = int(std::ceil(maxtau+2.718281828459045));
      int tap_count = endat-startat+1;
      tap_count = (tap_count%2) ? tap_count : tap_count + 1;
      std::vector<int> N(tap_count);
      for(size_t idx = 0; idx < N.size(); idx++){
        N[idx] = startat + idx;
      }

      d_sinc = std::vector<double>(d_path_count*N.size());
      double maxsinc = 0.;
      //std::cout << "RC: tau("<<d_taus.size()<<"), N("<<N.size()<<"), sinc("<<d_sinc.size()<<")\n";
      for(size_t nidx = 0; nidx < N.size(); nidx++){
        for(size_t tidx = 0; tidx < d_taus.size(); tidx++){
          d_sinc[nidx*d_taus.size() + tidx] = boost::math::sinc_pi(M_PI*(d_taus[tidx] - N[nidx]));
          maxsinc = (d_sinc[nidx*d_taus.size() + tidx] > maxsinc) ? d_sinc[nidx*d_taus.size() + tidx] : maxsinc;
        }
      }
      bool good_col = false;
      int col_g_s(0),col_g_e(N.size()-1);
      while(!good_col){
        for(size_t tidx = 0; tidx < d_taus.size(); tidx++){
          if(maxsinc*0.05 < d_sinc[col_g_s*d_taus.size()+tidx]){
            good_col = true;
          }
        }
        if(!good_col) col_g_s++;
      }
      good_col = false;
      while(!good_col){
        for(size_t tidx = 0; tidx < d_taus.size(); tidx++){
          if(maxsinc*0.01 < d_sinc[col_g_e*d_taus.size()+tidx]){
            good_col = true;
          }
        }
        if(!good_col) col_g_e--;
      }
      if(col_g_e+1 != N.size()){
        tap_count = col_g_e-col_g_s+1;
        tap_count = (tap_count%2) ? tap_count : tap_count + 1;
        N = std::vector<int>(col_g_e-col_g_s+1);
        for(size_t idx = 0; idx < N.size(); idx++){
          N[idx] = startat+idx;
        }

        d_sinc = std::vector<double>(d_path_count*N.size());
        for(size_t nidx = 0; nidx < N.size(); nidx++){
          for(size_t tidx = 0; tidx < d_taus.size(); tidx++){
            d_sinc[nidx*d_taus.size() + tidx] = boost::math::sinc_pi(M_PI*(d_taus[tidx] - N[nidx]));
          }
        }
      }
      d_d_tap_count = N.size();
      //std::cout << "RC: N("<<d_d_tap_count<<")\n";
      //std::cout << "RC: N("<<d_d_tap_count<<"), sinc_mat = [" << d_sinc[0];
      //for(size_t idx = 1; idx < d_sinc.size(); idx++){
      //  std::cout << ", "<<d_sinc[idx];
      //}
      //std::cout << "];\n";
    }

    void
    rayleigh_channel_impl::load_channel(int count)
    {
      //std::cout << "RC: loading channel" << std::endl;
      if(count > d_channel.size()){
        std::vector<complexf> digital_channel(d_d_tap_count,complexf(0.,0.));
        d_channel = std::vector< std::vector<complexf> >(count, digital_channel);
      }
      for(size_t idx = 0; idx < d_channel.size(); idx++){
        for(size_t aidx = 0; aidx < d_path_count; aidx++){
          d_analog_channel[aidx] = d_paths[aidx]->step();
        }
        if(d_lp){
          float larg = std::arg(d_analog_channel[0]);
          for(size_t aidx = 0; aidx < d_path_count; aidx++){
            d_analog_channel[aidx] *= std::exp(complexf(0.,-larg));
          }
        }
        sinc_multi( d_channel[idx], d_analog_channel );
      }
      if((d_norm)&&(d_d_tap_count==1)){
        double wd = 0.;
        size_t minidx = std::min(d_channel.size(), size_t(count));
        for(size_t idx = 0; idx < minidx; idx++){
          wd += (d_channel[idx][0] * std::conj(d_channel[idx][0])).real();
        }
        wd = 1./sqrt(wd);
        for(size_t idx = 0; idx < minidx; idx++){
          d_channel[idx][0] *= float(wd);
        }
      }
    }

    void
    rayleigh_channel_impl::sinc_multi(std::vector<complexf> &out, std::vector<complexf> &in)
    {
      for(size_t oidx = 0; oidx < out.size(); oidx++){
        out[oidx] = complexf(0.,0.);
        for(size_t iidx = 0; iidx < in.size(); iidx++){
          out[oidx] += in[iidx]*float(d_weights[iidx]*d_sinc[oidx*in.size() + iidx]);
        }
      }
      if((d_norm)&&(d_d_tap_count > 1)){
        double wd = 0.;
        for(size_t idx = 0; idx < out.size(); idx++){
          wd += (out[idx] * std::conj(out[idx])).real();
        }
        wd = 1./sqrt(wd);
        for(size_t idx = 0; idx < out.size(); idx++){
          out[idx] *= float(wd);
        }
      }
    }

  } /* namespace lte_channels */
} /* namespace gr */

