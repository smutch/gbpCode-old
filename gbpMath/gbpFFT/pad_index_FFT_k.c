#include <gbpLib.h>
#include <gbpFFT.h>

size_t pad_index_FFT_k(field_info *FFT,size_t index){
#if USE_MPI
  if(FFT->n_d==2)
    return(index+(size_t)FFT->pad_size_k*(index/(size_t)FFT->n_k_local[0]));
  else
    return(index+(size_t)FFT->pad_size_k*(index/(size_t)FFT->n_k_local[FFT->n_d-1]));
#else
  return(index+(size_t)FFT->pad_size_k*(index/(size_t)FFT->n_k_local[FFT->n_d-1]));
#endif
}
