
#include <lte_channels/doppler_base.h>


gr::lte_channels::Doppler_Base::Doppler_Base(size_t resolution, double fm, double sample_rate)
{
  d_rs = resolution;
  d_fm = fm;
  d_sr = sample_rate;
}

void
gr::lte_channels::Doppler_Base::gen_frequency()
{
  d_f = std::vector<float>(d_rs);
  float o = (1./float(d_rs))*(float(d_rs%2));
  float p = std::pow(2.,float(d_rs%2));
  float a = -.5 + o/p;
  float b = .5 - 1./float(d_rs)/p;
  float c = (b-a)/float(d_rs-1);
  for(size_t idx = 0; idx < d_rs; idx++){
    d_f[idx] = a+idx*c;
  }
}
