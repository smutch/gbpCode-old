#include <gbpLib.h>
#include <gbpFFT.h>

void compute_FFT(field_info *FFT){

  // Add the FFTW padding if we need to.  The log message
  //   is only needed if we do add a buffer.
  int flag_log_message=FALSE;
  if(add_buffer_FFT_R(FFT)){
     SID_log("Performing FFT...",SID_LOG_OPEN|SID_LOG_TIMER);
     flag_log_message=TRUE;
  }

  // Perform the FFT
  #if USE_MPI
    rfftwnd_mpi(FFT->plan,1,FFT->field_local,NULL,FFTW_TRANSPOSED_ORDER);
  #else
    rfftwnd_one_real_to_complex(FFT->plan,FFT->field_local,NULL);
  #endif

  if(flag_log_message)
     SID_log("Done.",SID_LOG_CLOSE);
}

