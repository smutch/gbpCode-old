#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void set_MCMC_iterations(MCMC_info *MCMC,int n_avg,int n_avg_covariance,int n_thin,int n_burn,int n_integrate){
  MCMC->n_avg            =n_avg;
  MCMC->n_avg_covariance =n_avg_covariance;
  MCMC->n_thin           =n_thin;
  MCMC->n_iterations_burn=n_burn;
  MCMC->n_iterations     =n_burn+n_integrate;
}

