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

void init_MCMC_arrays(MCMC_info *MCMC){
  int           i_DS;
  MCMC_DS_info *current_DS;
  MCMC_DS_info *next_DS;

  if(MCMC->n_M==NULL){
    SID_log("Initializing MCMC target arrays for %d dataset(s)...",SID_LOG_OPEN,MCMC->n_DS);

    // Initialize chain arrays
    MCMC->n_M      =(int      *)SID_malloc(sizeof(int      )*MCMC->n_DS);
    MCMC->M_new    =(double  **)SID_malloc(sizeof(double  *)*MCMC->n_DS);
    MCMC->M_last   =(double  **)SID_malloc(sizeof(double  *)*MCMC->n_DS);
    MCMC->n_M_total=0;
    i_DS           =0;
    current_DS     =MCMC->DS;
    while(current_DS!=NULL){
      next_DS           =current_DS->next;
      MCMC->n_M[i_DS]   =current_DS->n_M;
      MCMC->M_new[i_DS] =(double *)SID_malloc(sizeof(double)*MCMC->n_M[i_DS]);
      MCMC->M_last[i_DS]=(double *)SID_malloc(sizeof(double)*MCMC->n_M[i_DS]);
      MCMC->n_M_total  +=MCMC->n_M[i_DS];
      current_DS        =next_DS;
      i_DS++;
    }

    // Initialize results arrays
    MCMC->P_min  =(double  *)SID_malloc(sizeof(double)*MCMC->n_P);
    MCMC->P_max  =(double  *)SID_malloc(sizeof(double)*MCMC->n_P);
    MCMC->P_avg  =(double  *)SID_malloc(sizeof(double)*MCMC->n_P);
    MCMC->dP_avg =(double  *)SID_malloc(sizeof(double)*MCMC->n_P);
    MCMC->P_best =(double  *)SID_malloc(sizeof(double)*MCMC->n_P);
    MCMC->P_peak =(double  *)SID_malloc(sizeof(double)*MCMC->n_P);
    MCMC->P_lo_68=(double  *)SID_malloc(sizeof(double)*MCMC->n_P);
    MCMC->P_hi_68=(double  *)SID_malloc(sizeof(double)*MCMC->n_P);
    MCMC->P_lo_95=(double  *)SID_malloc(sizeof(double)*MCMC->n_P);
    MCMC->P_hi_95=(double  *)SID_malloc(sizeof(double)*MCMC->n_P);

    MCMC->ln_likelihood_DS     =(double *)SID_malloc(sizeof(double)*MCMC->n_DS);
    MCMC->ln_likelihood_DS_best=(double *)SID_malloc(sizeof(double)*MCMC->n_DS);
    MCMC->ln_likelihood_DS_peak=(double *)SID_malloc(sizeof(double)*MCMC->n_DS);
    MCMC->n_DoF_DS             =(int    *)SID_malloc(sizeof(int)*MCMC->n_DS);
    MCMC->n_DoF_DS_best        =(int    *)SID_malloc(sizeof(int)*MCMC->n_DS);
    MCMC->n_DoF_DS_peak        =(int    *)SID_malloc(sizeof(int)*MCMC->n_DS);

    // If MCMC_MODE_MINIMIZE_IO is on, we need to allocate buffers to store the chain
    if(check_mode_for_flag(MCMC->mode,MCMC_MODE_MINIMIZE_IO)){
      MCMC->flag_success_buffer     =(char   *)SID_malloc(sizeof(char)  *MCMC->n_iterations*MCMC->n_avg);
      MCMC->ln_likelihood_new_buffer=(double *)SID_malloc(sizeof(double)*MCMC->n_iterations*MCMC->n_avg);
      MCMC->P_new_buffer            =(double *)SID_malloc(sizeof(double)*MCMC->n_iterations*MCMC->n_avg*MCMC->n_P);
      MCMC->M_new_buffer            =(double *)SID_malloc(sizeof(double)*MCMC->n_iterations*MCMC->n_avg*MCMC->n_M_total);
    }

    SID_log("Done.",SID_LOG_CLOSE);
  }
}


