#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void set_MCMC_likelihood_function(MCMC_info *MCMC,void (*likelihood_function)(MCMC_info *,double **,double *,double *)){
  MCMC->compute_MCMC_ln_likelihood=likelihood_function;
}

