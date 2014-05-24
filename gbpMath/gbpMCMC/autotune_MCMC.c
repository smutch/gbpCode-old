#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMCMC.h>

void autotune_MCMC(MCMC_info *MCMC){
  int    i_tune;
  double temp_init;
  temp_init=MCMC->temperature;
  SID_log("Perform autotuning...",SID_LOG_OPEN|SID_LOG_TIMER);
  if(MCMC->n_autotune_randomize>0){
    SID_log("Performing initial state randomization...",SID_LOG_OPEN|SID_LOG_TIMER);
    if(check_mode_for_flag(MCMC->mode,MCMC_MODE_AUTOTUNE_REINIT_TEMP))
       set_MCMC_temperature(MCMC,temp_init);
    autotune_MCMC_temperature(MCMC);
    autotune_MCMC_randomize(MCMC);
    SID_log("Done.",SID_LOG_CLOSE);
  }
  for(i_tune=0;i_tune<MCMC->n_autotune;i_tune++){
    SID_log("Autotune iteration %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_tune+1,MCMC->n_autotune);
    if(check_mode_for_flag(MCMC->mode,MCMC_MODE_AUTOTUNE_REINIT_TEMP))
       set_MCMC_temperature(MCMC,temp_init);
    autotune_MCMC_temperature(MCMC);
    autotune_MCMC_covariance(MCMC);
    SID_log("Done.",SID_LOG_CLOSE);
  }
  if(check_mode_for_flag(MCMC->mode,MCMC_MODE_AUTOTUNE_REINIT_TEMP))
     set_MCMC_temperature(MCMC,temp_init);
  autotune_MCMC_temperature(MCMC);
  SID_log("Done.",SID_LOG_CLOSE);
}
