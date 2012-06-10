#include <gbpLib.h>
#include <gbpFFT.h>

void set_FFT_padding_state(field_info *FFT,int mode){
  FFT->flag_padded=mode;
}

