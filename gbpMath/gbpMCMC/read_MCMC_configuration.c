#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>

void read_MCMC_configuration(MCMC_info *MCMC,char *filename_output_dir,int chain){
  SID_log("Reading the configuration of chain #%d from {%s}...",SID_LOG_OPEN,chain);
  
  // Set the names of the files we have to read
  char filename_chain_dir[MAX_FILENAME_LENGTH];
  char filename_run[MAX_FILENAME_LENGTH];
  char filename_chain_config[MAX_FILENAME_LENGTH];
  sprintf(filename_chain_dir,   "%s/chains/",              filename_output_dir);
  sprintf(filename_run,         "%s/run.dat",              filename_output_dir);
  sprintf(filename_chain_config,"%s/chain_config_%06d.dat",filename_chain_dir,chain);

  // Read the run file to do some sanity checks
  int   n_chains;
  int   n_avg;
  int   flag_autocor_on;
  int   flag_no_map_write;
  int   n_P;
  FILE *fp;
  fp=fopen(filename_run,"r");
  fread(&n_chains,         sizeof(int),1,fp);
  fread(&n_avg,            sizeof(int),1,fp);
  fread(&flag_autocor_on,  sizeof(int),1,fp);
  fread(&flag_no_map_write,sizeof(int),1,fp);
  fread(&n_P,              sizeof(int),1,fp);
  fclose(fp);
  if(n_P!=MCMC->n_P)
    SID_trap_error("The number of parameters for the two runs does not match (ie. %d!=%d)",ERROR_LOGIC,n_P,MCMC->n_P);
  if(chain>=n_chains)
    SID_trap_error("Invalid chain number (ie. %d>%d)",ERROR_LOGIC,chain,n_chains);

  // Read the configuration file
  int     n_iterations_file_total;
  int     n_iterations_file_burn;
  double *V;
  double  temperature;
  V =(double *)SID_malloc(sizeof(double)*n_P*n_P);
  fp=fopen(filename_chain_config,"r");
  fread(&n_iterations_file_total,sizeof(int),   1,      fp);
  fread(&n_iterations_file_burn, sizeof(int),   1,      fp);
  fread(V,                       sizeof(double),n_P*n_P,fp);
  fread(&temperature,            sizeof(double),1,      fp);
  fclose(fp);

  // Set the new state
  set_MCMC_covariance(MCMC, V);
  set_MCMC_temperature(MCMC,temperature);

  // Clean-up
  SID_free(SID_FARG V);
  SID_log("Done.",SID_LOG_CLOSE);
}

