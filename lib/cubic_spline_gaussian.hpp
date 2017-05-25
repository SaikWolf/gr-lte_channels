

#ifndef INCLUDED_LTE_CHANNELS_CUBIC_SPLINE_GAUSSIAN_H
#define INCLUDED_LTE_CHANNELS_CUBIC_SPLINE_GAUSSIAN_H

#include <vector>
#include <complex>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include <boost/random/random_device.hpp>
#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/random.h>
#include <volk/volk.h>


namespace gr{
  namespace lte_channels{

    typedef std::complex<float> complexf;

    class Cubic_Spline_Gaussian
    {
     private:
      boost::random_device    d_rd;
      gr::random             *d_rng;

      std::vector<float>      d_taps;
      float                   d_delta;
      int                     d_seed;

      float                   d_last;

      int                     d_align;
      std::vector<complexf>   d_history;
      std::vector<complexf>   d_active;
      complexf               *d_buff;
      gr::filter::kernel::fir_filter_ccf *d_fir;

      void load_fir();
      void gen_history();

      void load_active();
      void load_next();

     public:
      Cubic_Spline_Gaussian(std::vector<float> &taps, float delta, int seed);
      ~Cubic_Spline_Gaussian();

      complexf step();
      void stepN(complexf* out, size_t N);
    };
  }//lte_channels
}//gr


#endif /* INCLUDED_LTE_CHANNELS_CUBIC_SPLINE_GAUSSIAN_HPP */
