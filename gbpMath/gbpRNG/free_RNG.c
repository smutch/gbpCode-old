#include <gbpCommon.h>
#include <gbpRNG.h>

void free_RNG(RNG_info *RNG){
  #ifdef USE_MPI
    #ifdef USE_SPRNG
      free_sprng(RNG->stream);
    #endif
  #endif
  RNG->initialized=FALSE;
}
