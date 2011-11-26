#include <gbpCommon.h>
#include <gbpRNG.h>

void free_RNG(RNG_info *RNG){
  #if USE_MPI
    #if USE_SPRNG
      free_sprng(RNG->stream);
    #endif
  #endif
  RNG->initialized=FALSE;
}
