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

void free_MCMC_arrays(MCMC_info *MCMC){
  int           i_DS;
  MCMC_DS_info *current_DS;
  MCMC_DS_info *next_DS;

  if(MCMC->n_M!=NULL){
    SID_log("Freeing MCMC target arrays for %d dataset(s)...",SID_LOG_OPEN,MCMC->n_DS);

    // Free chain arrays
    i_DS           =0;
    current_DS     =MCMC->DS;
    while(current_DS!=NULL){
      next_DS           =current_DS->next;
      SID_free(SID_FARG MCMC->M_new[i_DS]);
      SID_free(SID_FARG MCMC->M_last[i_DS]);
      current_DS        =next_DS;
      i_DS++;
    }
    SID_free(SID_FARG MCMC->n_M);
    SID_free(SID_FARG MCMC->M_new);
    SID_free(SID_FARG MCMC->M_last);
    MCMC->n_M_total=0;

    // Free results arrays
    SID_free(SID_FARG MCMC->P_min);
    SID_free(SID_FARG MCMC->P_max);
    SID_free(SID_FARG MCMC->P_avg);
    SID_free(SID_FARG MCMC->dP_avg);
    SID_free(SID_FARG MCMC->P_best);
    SID_free(SID_FARG MCMC->P_peak);
    SID_free(SID_FARG MCMC->P_lo_68);
    SID_free(SID_FARG MCMC->P_hi_68);
    SID_free(SID_FARG MCMC->P_lo_95);
    SID_free(SID_FARG MCMC->P_hi_95);

    SID_free(SID_FARG MCMC->ln_likelihood_DS);
    SID_free(SID_FARG MCMC->ln_likelihood_DS_best);
    SID_free(SID_FARG MCMC->ln_likelihood_DS_peak);
    SID_free(SID_FARG MCMC->n_DoF_DS);
    SID_free(SID_FARG MCMC->n_DoF_DS_best);
    SID_free(SID_FARG MCMC->n_DoF_DS_peak);

    // These only need to be deallocated if we are using MCMC_MODE_MINIMIZE_IO
    if(check_mode_for_flag(MCMC->mode,MCMC_MODE_MINIMIZE_IO)){
      SID_free(SID_FARG MCMC->flag_success_buffer);
      SID_free(SID_FARG MCMC->ln_likelihood_new_buffer);
      SID_free(SID_FARG MCMC->P_new_buffer);
      SID_free(SID_FARG MCMC->M_new_buffer);
    }

    SID_log("Done.",SID_LOG_CLOSE);
  }
}


