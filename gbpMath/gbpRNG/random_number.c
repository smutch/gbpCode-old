#include <gbpLib.h>
#include <gbpRNG.h>

REAL random_number(RNG_info *RNG){
  REAL random_local;
  if(RNG->initialized){
    #if USE_MPI
      if(RNG->global){
        if(SID.I_am_Master){
          #if USE_SPRNG
            random_local=(REAL)sprng(RNG->stream);
          #else
            random_local=(REAL)ran1(&(RNG->stream));
          #endif
        }
        MPI_Bcast(&random_local,1,MPI_REAL,MASTER_RANK,SID.COMM_WORLD->comm);
        return(random_local);
      }
      else{
        #if USE_SPRNG
          return((REAL)(sprng(RNG->stream)));
        #else
          return((REAL)ran1(&(RNG->stream)));
        #endif
      }
    #else
      #if USE_SPRNG
        return((REAL)(sprng()));
      #else
        return((REAL)ran1(&(RNG->stream)));
      #endif
    #endif
  }
  else
    SID_trap_error("RNG_info not initialized in call to random_number.",ERROR_LOGIC);
}
