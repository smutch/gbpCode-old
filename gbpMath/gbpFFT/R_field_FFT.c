#include <gbpLib.h>
#include <gbpFFT.h>

void R_field_FFT(field_info *FFT,int *i_R,double *R_field){
  int    i_d;
  for(i_d=0;i_d<FFT->n_d;i_d++)
    R_field[i_d]=FFT->R_field[i_d][FFT->i_R_start_local[i_d]+i_R[i_d]];
}
