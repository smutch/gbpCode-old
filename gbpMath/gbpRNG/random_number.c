#include <gbpLib.h>
#include <gbpRNG.h>

REAL random_number(RNG_info *RNG){
  REAL random_local;
  if(RNG->initialized){
    #ifdef USE_MPI
      if(RNG->global){
        if(SID.I_am_Master)
          random_local=(REAL)sprng(RNG->stream);
        MPI_Bcast(&random_local,1,MPI_REAL,MASTER_RANK,MPI_COMM_WORLD);
        return(random_local);
      }
      else{
        return((REAL)(sprng(RNG->stream)));
      }
    #else
      return((REAL)(sprng()));
    #endif
  }
  else
    SID_trap_error("RNG_info not initialized in call to random_number.",ERROR_LOGIC);
}
