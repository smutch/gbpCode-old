#define  _MAIN
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpMCMC.h>

int main(int argc, char *argv[]){
  char      filename_root[MAX_FILENAME_LENGTH];
  int       coverage_size=100;
  MCMC_info MCMC;
  
  SID_init(&argc,&argv,NULL);

  strcpy(filename_root,argv[1]);
  if(argc>2)
    coverage_size=atoi(argv[2]);
 
  SID_log("Performing analysis of MCMC dataset {%s}...",SID_LOG_OPEN,filename_root);

  read_MCMC_state(&MCMC,filename_root);
  set_MCMC_coverage_size(&MCMC,coverage_size);
  analyze_MCMC(&MCMC);

  // Clean-up
  SID_log("Cleaning-up...",SID_LOG_OPEN);
  free_MCMC(&MCMC);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
