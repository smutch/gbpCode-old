#include <gbpLib.h>
#include <gbpFFT.h>

void add_buffer_FFT_k(field_info *FFT){
  size_t index;
  size_t i_FFT;
  if(!(FFT->flag_padded)){
    SID_log("Adding FFTW padding...",SID_LOG_OPEN);
    for(i_FFT=FFT->n_field_k_local-1;i_FFT>0;i_FFT--){
      index=pad_index_FFT_k(FFT,i_FFT);
      if(index!=i_FFT){
        FFT->cfield_local[index].re=FFT->cfield_local[i_FFT].re;
        FFT->cfield_local[index].im=FFT->cfield_local[i_FFT].im;
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
  FFT->flag_padded=TRUE;
}

