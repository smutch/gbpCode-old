#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void rm_MCMC_directory(MCMC_info *MCMC){
  char command_text[2*MAX_FILENAME_LENGTH];
  sprintf(command_text,"rm -rf %s",MCMC->filename_output_dir);
  system(command_text);
}

