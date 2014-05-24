#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void set_MCMC_directory(MCMC_info *MCMC,const char *directory){
  strcpy(MCMC->filename_output_dir,directory);
}

