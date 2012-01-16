#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void autotune_MCMC_randomize(MCMC_info *MCMC){
  double ln_likelihood_last;
  int    i_iteration;
  time_t time_start, time_stop;
  double time_diff;

  SID_log("Creating faux chain (%d iterations)...",SID_LOG_OPEN|SID_LOG_TIMER,MCMC->n_autotune_randomize);

  // Start from the initial conditions
  MCMC->flag_init_chain=TRUE;

  // Start timer
  time(&time_start);

  // Create a faux chain from the given starting point to create a new starting point
  for(i_iteration=0;i_iteration<MCMC->n_autotune_randomize;i_iteration++)
     generate_MCMC_chain(MCMC);

  // Create arrays for the new state, the old state and the change
  int     i_P;
  double *P_new;
  double *P_old;
  double *P_diff;
  P_new =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
  P_old =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
  P_diff=(double *)SID_malloc(sizeof(double)*MCMC->n_P);
  memcpy(P_new,MCMC->P_last,(size_t)MCMC->n_P*sizeof(double));
  memcpy(P_old,MCMC->P_init,(size_t)MCMC->n_P*sizeof(double));
  for(i_P=0;i_P<MCMC->n_P;i_P++)
     P_diff[i_P]=100.*(P_new[i_P]-P_old[i_P])/P_old[i_P];

  // Set the new starting state
  memcpy(MCMC->P_init,P_new,(size_t)MCMC->n_P*sizeof(double));

  // Report the results
  char   format_string[128];
  if(MCMC->my_chain==SID.My_rank){
     SID_log("New initial parameters:",SID_LOG_ALLRANKS|SID_LOG_OPEN);
     sprintf(format_string,"%s = %%13.6le (was %%13.6le; change=%%7.2lf%%%)",MCMC->P_name_format);
     for(i_P=0;i_P<MCMC->n_P;i_P++)
        SID_log(format_string,SID_LOG_ALLRANKS|SID_LOG_COMMENT,MCMC->P_names[i_P],P_new[i_P],P_old[i_P],P_diff[i_P]);
     SID_log("",SID_LOG_ALLRANKS|SID_LOG_NOPRINT|SID_LOG_CLOSE);
  }

  // Stop timer and report model call speed
  time(&time_stop);
  time_diff  = difftime(time_stop,time_start);
  time_diff /= (double)MCMC->n_autotune_randomize;
  if(time_diff>0.1)
     SID_log("Mean time for single model call: %.1f", SID_LOG_COMMENT, time_diff);

  // Clean-up
  SID_free(SID_FARG P_new);
  SID_free(SID_FARG P_old);
  SID_free(SID_FARG P_diff);

  SID_log("Done.",SID_LOG_CLOSE);
}

