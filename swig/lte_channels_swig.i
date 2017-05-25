/* -*- c++ -*- */

#define LTE_CHANNELS_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "lte_channels_swig_doc.i"

%{
#include "lte_channels/rayleigh_channel.h"
%}


%include "lte_channels/rayleigh_channel.h"
GR_SWIG_BLOCK_MAGIC2(lte_channels, rayleigh_channel);
