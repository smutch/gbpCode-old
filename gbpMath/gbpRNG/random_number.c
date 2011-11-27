#include <gbpLib.h>
#include <gbpRNG.h>

GBPREAL random_number(RNG_info *RNG){
  GBPREAL random_local;
  if(RNG->initialized){
    #if USE_MPI
      if(RNG->global){
        if(SID.I_am_Master){
          #if USE_SPRNG
            random_local=(GBPREAL)sprng(RNG->stream);
          #else
            random_local=(GBPREAL)ran1(&(RNG->stream));
          #endif
        }
        MPI_Bcast(&random_local,1,MPI_REAL,MASTER_RANK,SID.COMM_WORLD->comm);
        return(random_local);
      }
      else{
        #if USE_SPRNG
          return((GBPREAL)(sprng(RNG->stream)));
        #else
          return((GBPREAL)ran1(&(RNG->stream)));
        #endif
      }
    #else
      #if USE_SPRNG
        return((GBPREAL)(sprng()));
      #else
        return((GBPREAL)ran1(&(RNG->stream)));
      #endif
    #endif
  }
  else
    SID_trap_error("RNG_info not initialized in call to random_number.",ERROR_LOGIC);
}
