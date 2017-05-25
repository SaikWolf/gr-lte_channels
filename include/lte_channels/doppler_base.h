

#ifndef INCLUDED_LTE_CHANNELS_DOPPLER_BASE_H
#define INCLUDED_LTE_CHANNELS_DOPPLER_BASE_H

#include <complex>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gnuradio/fft/fft.h>


namespace gr{
  namespace lte_channels{

    typedef std::complex<float> complexf;

    class Doppler_Base
    {
     protected:
      double d_fm;
      double d_sr;
      size_t d_rs;

      fft::fft_complex* d_fft;

      std::vector<float> d_f;

      void gen_frequency();

     public:
      Doppler_Base(size_t resolution, double fm, double sample_rate);
      virtual ~Doppler_Base() = 0;

      virtual void set_resolution(size_t resolution) = 0;
      virtual void set_max_doppler(double fm) = 0;
      virtual void set_sample_rate(double samp_rate) = 0;

      virtual void get_taps(std::vector<float>& taps) = 0;
    };

    inline Doppler_Base::~Doppler_Base()
    {}

  }//lte_channels
}//gr

#endif /* INCLUDED_LTE_CHANNELS_DOPPLER_BASE_HPP */
