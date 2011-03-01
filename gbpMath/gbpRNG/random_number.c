#include <gbpLib.h>
#include <gbpRNG.h>

REAL random_number(RNG_info *RNG){
  REAL random_local;
  if(RNG->initialized){
    #ifdef USE_MPI
      if(RNG->global){
        if(SID.I_am_Master){
          #ifdef USE_SPRNG
            random_local=(REAL)sprng(RNG->stream);
          #else
            random_local=(REAL)ran1(&(RNG->stream));
          #endif
        }
        MPI_Bcast(&random_local,1,MPI_REAL,MASTER_RANK,MPI_COMM_WORLD);
        return(random_local);
      }
      else{
        #ifdef USE_SPRNG
          return((REAL)(sprng(RNG->stream)));
        #else
          return((REAL)ran1(&(RNG->stream)));
        #endif
      }
    #else
      #ifdef USE_SPRNG
        return((REAL)(sprng()));
      #else
        return((REAL)ran1(&(RNG->stream)));
      #endif
    #endif
  }
  else
    SID_trap_error("RNG_info not initialized in call to random_number.",ERROR_LOGIC);
}
