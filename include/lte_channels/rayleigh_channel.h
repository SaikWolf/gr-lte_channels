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


#ifndef INCLUDED_LTE_CHANNELS_RAYLEIGH_CHANNEL_H
#define INCLUDED_LTE_CHANNELS_RAYLEIGH_CHANNEL_H

#include <lte_channels/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace lte_channels {

    /*!
     * \brief Simulation of a Rayleigh Time-varying Channel
     * \ingroup lte_channels
     *
     */
    class LTE_CHANNELS_API rayleigh_channel : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<rayleigh_channel> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of lte_channels::rayleigh_channel.
       *
       * To avoid accidental use of raw pointers, lte_channels::rayleigh_channel's
       * constructor is in a private implementation
       * class. lte_channels::rayleigh_channel::make is the public interface for
       * creating new instances.
       */
      static sptr make(double samp_rate, double maximum_doppler, size_t component_count, const std::vector<double> &delays, const std::vector<double> &weights, size_t doppler_spectrum=0, size_t doppler_taps=64, bool normalize=true, bool lead_phase=true, int seed=-1);

      virtual int get_tap_len() = 0;
    };

  } // namespace lte_channels
} // namespace gr

#endif /* INCLUDED_LTE_CHANNELS_RAYLEIGH_CHANNEL_H */

