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

void read_MCMC_covariance(MCMC_info *MCMC,char *filename){
  FILE   *fp;
  double *V;
  int     n_P,i_P,j_P;
  if(MCMC->V==NULL)
    SID_log("Initializing the covariance matrix from file {%s}...",SID_LOG_OPEN,filename);
  else
    SID_log("Updating the covariance matrix from file {%s}...",SID_LOG_OPEN,filename);
  fp=fopen(filename,"r");
  fread(&n_P,sizeof(int),1,fp);
  V=(double *)SID_malloc(sizeof(double)*n_P*n_P);
  fread(V,sizeof(double),n_P*n_P,fp);
  fclose(fp);
  set_MCMC_covariance(MCMC,V);
  SID_free(SID_FARG V);
  SID_log("Done.",SID_LOG_CLOSE);
}

