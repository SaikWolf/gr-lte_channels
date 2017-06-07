
#include "cubic_spline_gaussian.hpp"
#include <iostream>

namespace gr{
  namespace lte_channels{

    Cubic_Spline_Gaussian::Cubic_Spline_Gaussian(std::vector<float> &taps, float delta, int seed)
    : d_delta(delta),
      d_seed(seed)
    {
      d_rd = new boost::random_device();
      //std::cout << "CSG: constructing a cubic spline\n";
      d_taps = std::vector<float>(taps.begin(), taps.end());
      if(d_seed < 0){
        d_seed = (*d_rd)();
      }
      delete d_rd;
      d_rng = new gr::random(d_seed, 0, 1);
      //std::cout << "CSG: using seed "<<d_seed<<std::endl;

      d_align = volk_get_alignment();

      load_fir();
      gen_history();

      d_last = -delta;

      load_active();

    }

    Cubic_Spline_Gaussian::~Cubic_Spline_Gaussian()
    {
      delete d_rng;
      delete d_fir;
    }

    void
    Cubic_Spline_Gaussian::load_fir()
    {
      d_fir = new gr::filter::kernel::fir_filter_ccf(1,d_taps);
    }

    void
    Cubic_Spline_Gaussian::gen_history()
    {
      //std::cout << "CSG: building a history of length ";
      d_history = std::vector<complexf>(d_taps.size()+3);
      //std::cout << d_history.size() << std::endl;
      for(size_t idx = 0; idx < d_history.size(); idx++){
        d_history[idx] = d_rng->rayleigh_complex();
      }
    }

    complexf
    Cubic_Spline_Gaussian::step()
    {
      complexf out;
      float n = d_last + d_delta;
      while(n >= 1.0){
        load_next();
        n = n-1.;
      }

      float n3 = n*n*n;
      float n2 = n*n;

      out = (2.f*n3 - 3.f*n2 + 1.f)*d_active[1]
          + 0.5f*(n3 - 2.f*n2 + n)*(d_active[2]-d_active[0])
          + (-2.f*n3 + 3.f*n2)*d_active[2]
          + 0.5f*(n3 - n2)*(d_active[3]-d_active[1]);

      d_last = n;
      return out;
    }

    void
    Cubic_Spline_Gaussian::stepN(complexf* out, size_t N)
    {
      for(size_t idx = 0; idx < N; idx++){
        out[idx] = step();
      }
    }

    void
    Cubic_Spline_Gaussian::load_active()
    {
      d_buff = (complexf*) volk_malloc( d_history.size()*sizeof(complexf), d_align );
      memcpy( d_buff, &d_history[0], d_history.size()*sizeof(complexf) );
      d_active = std::vector<complexf>(4);
      for(size_t idx = 0; idx < 4; idx++){
        d_active[idx] = d_fir->filter(&d_buff[idx]);
      }
      volk_free(d_buff);
      //std::cout << "CSG: The active buffer contains: active = [" << d_active[0]
      //          << ", " << d_active[1] << ", " << d_active[2] << ", " << d_active[3]
      //          << "];\n";
    }

    void
    Cubic_Spline_Gaussian::load_next()
    {//clearly this can be simplified but later
      for(size_t idx = 0; idx < 3; idx++){
        d_active[idx] = d_active[idx+1];
      }
      d_buff = (complexf*) volk_malloc( d_history.size()*sizeof(complexf), d_align );
      memcpy( d_buff, &d_history[1], (d_history.size()-1)*sizeof(complexf) );
      d_buff[d_history.size()-1] = d_rng->rayleigh_complex();
      d_active[3] = d_fir->filter(&d_buff[3]);
      d_history = std::vector<complexf>( d_buff, d_buff+d_history.size() );
      volk_free(d_buff);
    }


  }//lte_channels
}//gr
