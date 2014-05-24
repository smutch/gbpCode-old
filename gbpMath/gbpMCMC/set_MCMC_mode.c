#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void set_MCMC_mode(MCMC_info *MCMC,int mode){
  MCMC->mode=mode;

  // Set the local chain number
  if(check_mode_for_flag(MCMC->mode,MCMC_MODE_PARALLEL)){
    MCMC->my_chain=SID.My_rank;
    MCMC->n_chains=SID.n_proc;
  }
  else{
    MCMC->my_chain=MASTER_RANK;
    MCMC->n_chains=1;
  }
  MCMC->flag_no_map_write=check_mode_for_flag(mode,MCMC_MODE_NO_MAP_WRITE);

}

