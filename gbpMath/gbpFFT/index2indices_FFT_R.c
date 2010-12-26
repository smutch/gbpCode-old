#include <gbpLib.h>
#include <gbpFFT.h>

void index2indices_FFT_R(field_info *FFT,size_t index_out,int *i_R){
  int    i_d;
  size_t remainder;
  for(i_d=FFT->n_d-1,remainder=index_out;i_d>=0;i_d--){
    i_R[i_d]=remainder%FFT->n[i_d];
    remainder-=i_R[i_d];
    remainder/=FFT->n[i_d];
  }
}
