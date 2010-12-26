#include <gbpCommon.h>
#include <gbpRNG.h>

void free_RNG(RNG_info *RNG){
#ifdef USE_MPI
  free_sprng(RNG->stream);
#endif
  RNG->initialized=FALSE;
}
