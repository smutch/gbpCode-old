#include <gbpLib.h>
#include <gbpRNG.h>

void init_RNG(int *seed,RNG_info *RNG,int mode){
  if((*seed)==0)
    init_seed_from_clock(seed);
  RNG->seed  =(*seed);
  #ifdef USE_MPI
    if(check_mode_for_flag(mode,RNG_GLOBAL)){
      RNG->seed =MIN(RNG->seed,(-1)*RNG->seed);
      RNG->global=TRUE;
    }
    else{
      RNG->stream=init_sprng(SID.My_rank,SID.n_proc,RNG->seed,SPRNG_DEFAULT);
      RNG->global=FALSE;
    }
  #else
    RNG->stream=init_sprng(RNG->seed,SPRNG_DEFAULT);
  #endif
  RNG->IGauss     =0;
  RNG->GaussBak   =0.;
  RNG->initialized=TRUE;
}
