#include <gbpLib.h>
#include <gbpFFT.h>

size_t pad_index_FFT_R(field_info *FFT,size_t index_in){
  return(index_in+(size_t)FFT->pad_size_R*(index_in/(size_t)FFT->n_R_local[FFT->n_d-1]));
}
