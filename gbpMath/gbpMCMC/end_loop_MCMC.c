#include <stdio.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void end_loop_MCMC(MCMC_info *MCMC){
  rm_MCMC_directory(MCMC);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
}

