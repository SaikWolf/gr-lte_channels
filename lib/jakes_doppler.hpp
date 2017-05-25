

#ifndef INCLUDED_LTE_CHANNELS_JAKES_DOPPLER_H
#define INCLUDED_LTE_CHANNELS_JAKES_DOPPLER_H

#include <lte_channels/doppler_base.h>


namespace gr{
  namespace lte_channels{

    class Jakes_Doppler : public Doppler_Base
    {
     private:

     public:
      Jakes_Doppler();
      Jakes_Doppler(size_t resolution, double fm, double samp_rate);
      ~Jakes_Doppler();

      void set_resolution(size_t resolution);
      void set_max_doppler(double fm);
      void set_sample_rate(double samp_rate);

      void get_taps(std::vector<float>& taps);
    };
  }//lte_channels
}//gr


#endif /* INCLUDED_LTE_CHANNELS_JAKES_DOPPLER_HPP */
