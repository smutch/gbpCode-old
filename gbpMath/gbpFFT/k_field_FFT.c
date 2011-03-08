#include <gbpLib.h>
#include <gbpFFT.h>

void k_field_FFT(field_info *FFT,int *i_k,double *k_field){
  int    i_d,j_d;
  for(j_d=0;j_d<FFT->n_d;j_d++){
    i_d=j_d;
#if USE_MPI
    if(j_d==1)      i_d=0;
    else if(j_d==0) i_d=1;
#endif
    k_field[i_d]=FFT->k_field[i_d][i_k[i_d]];
  }
}

