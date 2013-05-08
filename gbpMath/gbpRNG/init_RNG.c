#include <gbpLib.h>
#include <gbpRNG.h>

void init_RNG(int *seed,RNG_info *RNG,int mode){
  if((*seed)<=0)
    init_seed_from_clock(seed);
  RNG->seed=(*seed);
  #if USE_SPRNG
    #if USE_MPI
      if(check_mode_for_flag(mode,RNG_GLOBAL))
        RNG->global=TRUE;
      else
        RNG->global=FALSE;
      RNG->stream=init_sprng(SID.My_rank,SID.n_proc,RNG->seed,SPRNG_DEFAULT);
    #else
      RNG->stream=init_sprng(RNG->seed,SPRNG_DEFAULT);
    #endif
  #else
    #if USE_MPI
      if(check_mode_for_flag(mode,RNG_GLOBAL))
        RNG->global=TRUE;
      else
        RNG->global=FALSE;
    #endif
    if((*seed)<0)
      RNG->stream= (long)(RNG->seed-(long)(173*SID.My_rank));
    else
      RNG->stream=-(long)(RNG->seed-(long)(173*SID.My_rank));
  #endif
  RNG->IGauss     =0;
  RNG->GaussBak   =0.;
  RNG->initialized=TRUE;
}
