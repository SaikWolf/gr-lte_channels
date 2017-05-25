

#include "jakes_doppler.hpp"


namespace gr{
  namespace lte_channels{

    Jakes_Doppler::Jakes_Doppler()
    : Doppler_Base(64,.5,1.)
    {
      Doppler_Base::gen_frequency();
    }

    Jakes_Doppler::Jakes_Doppler(size_t resolution, double fm, double samp_rate)
    : Doppler_Base(resolution,fm,samp_rate)
    {
      Doppler_Base::gen_frequency();
    }

    Jakes_Doppler::~Jakes_Doppler()
    {
    }

    void
    Jakes_Doppler::set_resolution(size_t resolution)
    {
      d_rs = resolution;
      Doppler_Base::gen_frequency();
    }

    void
    Jakes_Doppler::set_max_doppler(double fm)
    {
      d_fm = fm;
    }

    void
    Jakes_Doppler::set_sample_rate(double samp_rate)
    {
      d_sr = samp_rate;
    }

    void
    Jakes_Doppler::get_taps(std::vector<float>& taps)
    {
      taps = std::vector<float>(d_rs,0.);

      d_fft = new gr::fft::fft_complex(d_rs,false);
      complexf *ib = d_fft->get_inbuf();
      complexf *ob = d_fft->get_outbuf();

      double lfm = d_fm/d_sr;
      double lfm2 = lfm*lfm;

      size_t offset = int(std::floor(float(d_rs)/2.));

      for(size_t idx = 0; idx < d_rs; idx++){
        if(std::abs(d_f[idx]) < lfm){
          ib[(idx+offset)%d_rs] = complexf(1./std::sqrt(M_PI*(lfm2 - d_f[idx]*d_f[idx])),0.);
        }
        else{
          ib[(idx+offset)%d_rs] = complexf(0.,0.);
        }
      }
      /*std::cout<< "JD: inbuf = ["<< ib[0].real();
      for(size_t idx = 1; idx < d_rs; idx++){
        std::cout << ", " << ib[idx].real();
      }
      std::cout << "];\n";*/

      d_fft->execute();

      int starter = int(std::ceil(float(d_rs)/2.));
      for(size_t idx = 0; idx < d_rs; idx++){
        taps[idx] = ob[(idx+starter)%d_rs].real();
      }

      double wt = 0.;
      for(size_t idx = 0; idx < taps.size(); idx++){
        wt += taps[idx]*taps[idx];
      }
      wt = 1./std::sqrt(wt);
      for(size_t idx = 0; idx < taps.size(); idx++){
        taps[idx] *= wt;
      }

      delete d_fft;
    }

  }//lte_channels
}//gr

