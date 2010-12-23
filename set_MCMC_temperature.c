#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void set_MCMC_temperature(MCMC_info *MCMC,double temperature){
  MCMC->temperature=temperature;
}

