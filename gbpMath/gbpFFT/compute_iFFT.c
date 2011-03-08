#include <gbpLib.h>
#include <gbpFFT.h>

void compute_iFFT(field_info *FFT){
  int i_k;
  int flag_padded;

  flag_padded=FFT->flag_padded;
  add_buffer_FFT_k(FFT);
  if(!flag_padded)
    SID_log("Performing iFFT...",SID_LOG_OPEN|SID_LOG_TIMER);
  #if USE_MPI
    rfftwnd_mpi(FFT->iplan,1,FFT->field_local,NULL,FFTW_TRANSPOSED_ORDER);
  #else
    rfftwnd_one_complex_to_real(FFT->iplan,FFT->cfield_local,NULL);
  #endif

  // Divide-out the FFTW scaling with N
  for(i_k=0;i_k<FFT->total_local_size;i_k++)
    FFT->field_local[i_k]/=(fftw_real)FFT->n_field;
  if(!flag_padded)
    SID_log("Done.",SID_LOG_CLOSE);  
  if(!flag_padded)
    remove_buffer_FFT_R(FFT);
}
