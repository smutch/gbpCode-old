#include <gbpLib.h>
#include <gbpFFT.h>

void compute_FFT(field_info *FFT){
  int flag_padded;
  flag_padded=FFT->flag_padded;
  add_buffer_FFT_R(FFT);
  if(FFT->flag_padded)
    SID_log("Performing FFT...",SID_LOG_OPEN|SID_LOG_TIMER);
  #if USE_MPI
    rfftwnd_mpi(FFT->plan,1,FFT->field_local,NULL,FFTW_TRANSPOSED_ORDER);
  #else
    rfftwnd_one_real_to_complex(FFT->plan,FFT->field_local,NULL);
  #endif
  if(!flag_padded)
    SID_log("Done.",SID_LOG_CLOSE);
  if(!flag_padded)
    remove_buffer_FFT_k(FFT);
}
