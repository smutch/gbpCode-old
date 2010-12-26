#include <gbpLib.h>
#include <gbpFFT.h>

void index2indices_FFT_k(field_info *FFT,size_t index,int *i_k){
  int    i_d,j_d;
  size_t remainder;
  for(j_d=FFT->n_d-1,remainder=index;j_d>=0;j_d--){
    i_d=j_d;
#ifdef USE_MPI
    if(j_d==0)      i_d=1;
    else if(j_d==1) i_d=0;
#endif
    i_k[i_d]  =remainder%FFT->n_k_local[i_d];
    remainder-=i_k[i_d];
    remainder/=FFT->n_k_local[i_d];
  }
}
