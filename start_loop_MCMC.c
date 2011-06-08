#include <stdio.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void start_loop_MCMC(MCMC_info *MCMC){
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
  restart_MCMC(MCMC);
  set_MCMC_mode(MCMC,MCMC->mode|MCMC_MODE_NO_MAP_WRITE|MCMC_MODE_MINIMIZE_IO);
  rm_MCMC_directory(MCMC);
}

