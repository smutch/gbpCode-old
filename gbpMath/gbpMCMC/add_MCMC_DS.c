#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMCMC.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>

void add_MCMC_DS(MCMC_info *MCMC,const char *name,int n_D,int *D,double *DS,double *dDS,void *params,int n_arrays,...){
  int           i_array;
  int           i_D;
  int           n_M;
  MCMC_DS_info *new_DS;
  va_list       vargs;
  va_start(vargs,n_arrays);

  SID_log("Adding data set {%s} to MCMC structure...",SID_LOG_OPEN,name);

  #if USE_CFITSIO==0
  if(n_D>1)
     SID_trap_error("You can only use 2D datasets if you compile with USE_CFITSIO=1 at the moment.",ERROR_LOGIC);
  #endif

  // Initialize new dataset
  init_MCMC_DS(&new_DS,name,n_D,D,DS,dDS,params,n_arrays,vargs);

  // Add new dataset to linked list
  if(MCMC->DS==NULL)
    MCMC->DS=new_DS;
  else
    MCMC->last->next=new_DS;
  MCMC->last=new_DS;

  MCMC->n_DS++;

  SID_log("Done.",SID_LOG_CLOSE);

  va_end(vargs);
}

