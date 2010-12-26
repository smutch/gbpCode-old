#include <gbpLib.h>
#include <gbpDomain.h>

void clear_field(field_info *FFT){
  size_t i_fft;
  for(i_fft=0;i_fft<FFT->total_local_size;i_fft++)
    FFT->field_local[i_fft]=0.;
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

