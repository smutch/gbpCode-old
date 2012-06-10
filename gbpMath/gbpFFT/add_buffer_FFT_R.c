#include <gbpLib.h>
#include <gbpFFT.h>

int add_buffer_FFT_R(field_info *FFT){
  int r_val=FALSE;
  if(!(FFT->flag_padded)){
    SID_log("Adding FFTW padding...",SID_LOG_OPEN);
    size_t i_FFT; // n.b.: Don't run this to the zero'th element, 
                  //       i_FFT is unsigned and the 'i_FFT--' causes problems.
                  //       There's no need since the first run of the last (fastest running)
                  //       dimension is obviously not padded.
    for(i_FFT=FFT->n_field_R_local-1;i_FFT>=FFT->n_R_local[FFT->n_d-1];i_FFT--){
      size_t index_pad;
      index_pad=pad_index_FFT_R(FFT,i_FFT);
      if(index_pad!=i_FFT){
        FFT->field_local[index_pad]=FFT->field_local[i_FFT];
      }
    }
    r_val=TRUE;
    SID_log("Done.",SID_LOG_CLOSE);
  }
  return(r_val);
}

