#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMCMC.h>

void free_MCMC_DS(MCMC_info *MCMC){
  int           i_array;
  int           i_D;
  int           i_DS;
  int           n_M;
  MCMC_DS_info *current_DS;
  MCMC_DS_info *next_DS;

  SID_log("Freeing data sets...",SID_LOG_OPEN);

  current_DS=MCMC->DS;
  i_DS=0;
  while(current_DS!=NULL){
    next_DS=current_DS->next;
    SID_free(SID_FARG current_DS->M_target);
    SID_free(SID_FARG current_DS->dM_target);
    SID_free(SID_FARG current_DS->M_best);
    SID_free(SID_FARG current_DS->M_best_parameters);
    SID_free(SID_FARG current_DS->M_peak_parameters);
    SID_free(SID_FARG current_DS->M_min);
    SID_free(SID_FARG current_DS->M_max);
    SID_free(SID_FARG current_DS->M_lo_68);
    SID_free(SID_FARG current_DS->M_hi_68);
    SID_free(SID_FARG current_DS->M_lo_95);
    SID_free(SID_FARG current_DS->M_hi_95);
    SID_free(SID_FARG current_DS->M_avg);
    SID_free(SID_FARG current_DS->dM_avg);
    if(MCMC->n_arrays>0){
      for(i_array=0;i_array<current_DS->n_arrays;i_array++){
        SID_free(SID_FARG current_DS->array[i_array]);
        SID_free(SID_FARG current_DS->array_name[i_array]);
      }
      SID_free(SID_FARG current_DS->array);
      SID_free(SID_FARG current_DS->array_name);
      SID_free(SID_FARG current_DS->D);
    }
    SID_free(SID_FARG current_DS);
    current_DS=next_DS;
    i_DS++;
  }
  MCMC->n_DS=0;
  MCMC->DS  =NULL;
  MCMC->last=NULL;

  SID_log("Done.",SID_LOG_CLOSE);
}

