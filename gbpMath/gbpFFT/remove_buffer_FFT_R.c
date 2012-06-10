#include <gbpLib.h>
#include <gbpFFT.h>

void remove_buffer_FFT_R(field_info *FFT){
  size_t i_FFT,index;
  if(!(FFT->flag_padded)){
    SID_log("Removing FFTW padding...",SID_LOG_OPEN);
    for(i_FFT=0;i_FFT<FFT->n_field_R_local;i_FFT++){
      index=pad_index_FFT_R(FFT,i_FFT);
      if(i_FFT!=index)
        FFT->field_local[i_FFT]=FFT->field_local[index];
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
}

