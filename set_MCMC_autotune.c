#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void set_MCMC_autotune(MCMC_info *MCMC,
                       double     success_target,
                       double     success_threshold,
                       double     covariance_threshold,
                       int        n_autotune,
                       int        n_autotune_temperature,
                       int        n_autotune_covariance){

  SID_log("Setting autotune parameters...",SID_LOG_OPEN);
  // Set the desired success rate
  if(success_target>0.)
    MCMC->success_target=success_target;
  else
    MCMC->success_target=MCMC_DEFAULT_SUCCESS_TARGET;

  // Set the success rate convergence tollerance
  if(success_threshold>0.)
    MCMC->success_threshold=success_threshold;
  else
    MCMC->success_threshold=MCMC_DEFAULT_SUCCESS_THRESH;

  // Set the covariance matrix convergence threshold
  if(covariance_threshold>0.)
    MCMC->covariance_threshold=covariance_threshold;
  else
    MCMC->covariance_threshold=MCMC_DEFAULT_COVARIANCE_THRESH;

  // Set the number of autotune iterations
  if(n_autotune>=0)
    MCMC->n_autotune=n_autotune;
  else
    MCMC->n_autotune=MCMC_DEFAULT_N_AUTOTUNE;

  // Set the number of iterations used to compute success rates
  //  for temperature autotuning
  if(n_autotune_temperature>0)
    MCMC->n_autotune_temperature=n_autotune_temperature;
  else
    MCMC->n_autotune_temperature=MCMC_DEFAULT_N_AUTOTUNE_TEMPERATURE;

  // Set the number of iterations used before checking/re-checking
  //   for the covariance matrix's convergence
  if(n_autotune_covariance>0)
    MCMC->n_autotune_covariance=n_autotune_covariance;
  else
    MCMC->n_autotune_covariance=MCMC_DEFAULT_N_AUTOTUNE_COVMTX;

  // If we're setting the autotune parameters, we must want to use it
  MCMC->mode=MCMC->mode|MCMC_MODE_AUTOTUNE;

  SID_log("success target        =%le",SID_LOG_COMMENT,MCMC->success_target);
  SID_log("success threshold     =%le",SID_LOG_COMMENT,MCMC->success_threshold);
  SID_log("covariance threshold  =%le",SID_LOG_COMMENT,MCMC->covariance_threshold);
  SID_log("n_autotune            =%d", SID_LOG_COMMENT,MCMC->n_autotune);
  SID_log("n_autotune_temperature=%d", SID_LOG_COMMENT,MCMC->n_autotune_temperature);
  SID_log("n_autotune_covariance =%d", SID_LOG_COMMENT,MCMC->n_autotune_covariance);

  SID_log("Done.",SID_LOG_CLOSE);
}

