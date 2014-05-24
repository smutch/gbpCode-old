#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void rm_MCMC_directory(MCMC_info *MCMC){
  char command_text[2*MAX_FILENAME_LENGTH];
  SID_log("Removing MCMC directory {%s}...",SID_LOG_OPEN,MCMC->filename_output_dir);
  sprintf(command_text,"rm -rf %s > /dev/null",MCMC->filename_output_dir);
  system(command_text);
  SID_log("Done.",SID_LOG_CLOSE);
}

