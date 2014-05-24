#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void set_MCMC_coverage_size(MCMC_info *MCMC,int coverage_size){
  MCMC->coverage_size=coverage_size;
}

