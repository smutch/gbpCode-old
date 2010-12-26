#include <gbpLib.h>
#include <gbpFFT.h>

size_t index_local_FFT_k(field_info *FFT,int *i_k){
  int    i_d,j_d;
  size_t index;
  for(j_d=0,index=0;j_d<FFT->n_d;j_d++){
    i_d=j_d;
#ifdef USE_MPI
    if(j_d==1)      i_d=0;
    else if(j_d==0) i_d=1;
#endif
    switch(i_d==FFT->n_d-1 && FFT->flag_padded){
      case FALSE:
        index*=FFT->n_k_local[i_d];
        break;
      case TRUE:
        index*=(FFT->n_k_local[i_d]+1);
        break;
    }
    index+=(i_k[i_d]-FFT->i_k_start_local[i_d]);
  }
  return(index);
}
