#include <gbpLib.h>
#include <gbpFFT.h>

void remove_buffer_FFT_k(field_info *FFT){
  size_t i_FFT,j_FFT,index;
  if(FFT->flag_padded){
     SID_log("Removing FFTW padding...",SID_LOG_OPEN);
    for(i_FFT=0;i_FFT<FFT->n_field_k_local;i_FFT++){
      index=pad_index_FFT_k(FFT,i_FFT);
      if(i_FFT!=index){
        FFT->cfield_local[i_FFT].re=FFT->cfield_local[index].re;
        FFT->cfield_local[i_FFT].im=FFT->cfield_local[index].im;
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
  FFT->flag_padded=FALSE;
}

