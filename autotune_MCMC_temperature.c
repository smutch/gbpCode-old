#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void autotune_MCMC_temperature(MCMC_info *MCMC){
  double success_target;
  double success_threshold;
  int    n_autotune;
  double success;
  double temperature_hi;
  double temperature_lo;
  double temperature;
  int    flag_raise_range;
  int    i_autotune;
  double ln_likelihood_last;
  int    i_iteration;
  int    n_success;

  SID_log("Autotuning the temperature (target=%6.2lf%%, threshold=%6.2lf%%)...",SID_LOG_OPEN|SID_LOG_TIMER,MCMC->success_target,MCMC->success_threshold);

  // Initialize
  success_target   =MCMC->success_target;
  success_threshold=MCMC->success_threshold;
  n_autotune       =MCMC->n_autotune_temperature;
  success          =0.;
  temperature_hi   =MCMC->temperature;
  temperature_lo   =0.;
  temperature      =temperature_hi;
  flag_raise_range =TRUE;
  i_autotune       =0;

  // Iterate with a bisection algorythm until convergence
  while(fabs(success-success_target)>success_threshold){
    // Perform bisection
    if(i_autotune>0){
      // If the success is too high, the temperature has to go up
      if(success>success_target){
        temperature_lo=temperature;
        if(flag_raise_range)
          temperature_hi*=4.;
      }
      // If the success is too low, the temperature has to go down
      else{
        temperature_hi  =temperature;
        flag_raise_range=FALSE;
      }
      temperature=0.5*(temperature_hi+temperature_lo);
    }
    set_MCMC_temperature(MCMC,temperature);

    // Start from the initial conditions for each iteration
    MCMC->flag_init_chain=TRUE;

    // Generate a success rate for the current trial temperature
    for(i_iteration=0,n_success=0;i_iteration<n_autotune;i_iteration++){
      if(generate_MCMC_proposition(MCMC)) 
        n_success++;
    }
    success=100.*(double)(n_success)/(double)n_autotune;
    i_autotune++;
    SID_log("Iteration #%04d: Temperature=%le Success=%5.1lf%%",SID_LOG_COMMENT,i_autotune,temperature,success);
  }

  SID_log("Done.",SID_LOG_CLOSE);
}

