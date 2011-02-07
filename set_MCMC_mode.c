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
  if(check_mode_for_flag(MCMC->mode,MCMC_MODE_PARALLEL))
    MCMC->my_chain=SID.My_rank;
  else
    MCMC->my_chain=MASTER_RANK;

}

