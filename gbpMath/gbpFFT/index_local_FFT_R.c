#include <gbpLib.h>
#include <gbpFFT.h>

size_t index_local_FFT_R(field_info *FFT,int *i_R){
  int    i_d;
  size_t index;
  for(i_d=0,index=0;i_d<FFT->n_d;i_d++){
    switch(i_d==FFT->n_d-1 && FFT->flag_padded){
      case FALSE:
        index*=FFT->n_R_local[i_d];
        break;
      case TRUE:
        index*=2*(FFT->n_R_local[i_d]/2+1);
        break;
    }
    if(i_R[i_d]<0 || i_R[i_d]>FFT->i_R_stop_local[i_d])
      SID_trap_error("Index (%d;i_d=%d) out of local slab's range (%d->%d).",ERROR_LOGIC,
                     i_R[i_d],i_d,FFT->i_R_start_local[i_d],FFT->i_R_stop_local[i_d]);
    index+=(i_R[i_d]-FFT->i_R_start_local[i_d]);
  }
  return(index);
}
