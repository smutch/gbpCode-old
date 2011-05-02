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

int generate_MCMC_chain(MCMC_info *MCMC){
  int              flag_success;
  int              i_P,i_DS;
  static char      format_string[256];
  FILE            *fp_report_props;
  char             filename_report_props[MAX_FILENAME_LENGTH];
  static int       n_P;
  static double   *P_new;
  static double  **M_new;
  static double   *P_last;
  static double   *P_chain;
  static double  **M_last;
  static double   *P_best;
  static double   *P_limit_min;
  static double   *P_limit_max;
  static RNG_info *RNG;
  static int       flag_report_props;
  static int       n_DS;
  static int      *n_M;

  // Initialize a few things on the first call
  switch(MCMC->first_chain_call){
    case TRUE:
      n_P        =MCMC->n_P;
      P_new      =MCMC->P_new;
      M_new      =MCMC->M_new;
      P_last     =MCMC->P_last;
      M_last     =MCMC->M_last;
      P_chain    =MCMC->P_chain;
      P_best     =MCMC->P_best;
      P_limit_min=MCMC->P_limit_min;
      P_limit_max=MCMC->P_limit_max;
      RNG        =MCMC->RNG;
      n_DS       =MCMC->n_DS;
      n_M        =MCMC->n_M;
      flag_report_props     =check_mode_for_flag(MCMC->mode,MCMC_MODE_REPORT_PROPS);
      MCMC->first_chain_call=FALSE;
      sprintf(filename_report_props,"%s/chains/report_props_%06d.dat",MCMC->filename_output_dir,MCMC->my_chain);
      break;
  }

  // We need to compute the likelihood of the initial conditions
  //   if we are just starting a chain
  if(MCMC->flag_init_chain)
    generate_MCMC_proposition(MCMC,TRUE);

  // Set the last new state to the new last state
  MCMC->ln_likelihood_last=MCMC->ln_likelihood_new;
  MCMC->ln_Pr_last        =MCMC->ln_Pr_new;
  memcpy(P_last,P_new,n_P*sizeof(double));
  for(i_DS=0;i_DS<n_DS;i_DS++)
    memcpy(M_last[i_DS],M_new[i_DS],n_M[i_DS]*sizeof(double));

  // Keep generating parameter sets until the mapping function is satisfied
  generate_MCMC_proposition(MCMC,FALSE);

  // Decide if this is a successful proposition or not...
  if(MCMC->ln_likelihood_new>MCMC->ln_likelihood_chain){
    MCMC->ln_Pr_new=0.;
    flag_success=TRUE;
  }
  else if(MCMC->my_chain==SID.My_rank){
    MCMC->ln_Pr_new=MCMC->ln_likelihood_new-MCMC->ln_likelihood_chain;
    if((double)random_number(RNG)<=exp(MCMC->ln_Pr_new))
      flag_success=TRUE;
    else
      flag_success=FALSE;
  }
  if(!check_mode_for_flag(MCMC->mode,MCMC_MODE_PARALLEL))
    SID_Bcast(&flag_success,sizeof(int),MASTER_RANK,SID.COMM_WORLD);

  // ... if it is, then update the chain ...
  if(flag_success){
    memcpy(P_chain,P_new,(size_t)n_P*sizeof(double));
    MCMC->ln_likelihood_chain=MCMC->ln_likelihood_new;
    MCMC->ln_Pr_chain        =MCMC->ln_Pr_new;
    MCMC->n_success++;
  }
  else
    MCMC->n_fail++;
  MCMC->n_propositions++;

  // Report the proposition if asked to
  if(flag_report_props && MCMC->my_chain==SID.My_rank){
     fp_report_props=fopen(filename_report_props,"w");
     if(flag_success)
        fprintf(fp_report_props,"Proposal #%09d: SUCCEEDED\n\n",MCMC->n_propositions);
     else
        fprintf(fp_report_props,"Proposal #%09d: REJECTED\n\n",MCMC->n_propositions);
     sprintf(format_string,"   %s = %%13.6le (was %%13.6le) (best is %%13.6le) (flat prior =%%13.6le -> %%13.6le)\n",MCMC->P_name_format);
     for(i_P=0;i_P<n_P;i_P++)
        fprintf(fp_report_props,format_string,
           MCMC->P_names[i_P],
           P_new[i_P],
           P_last[i_P],
           P_best[i_P],
           P_limit_min[i_P],
           P_limit_max[i_P]);
     fprintf(fp_report_props,"\n");
     fprintf(fp_report_props,"ln(likelihood) =%13.6le+constant (was %13.6le+constant) (best is %13.6le+constant)\n",
        MCMC->ln_likelihood_new,
        MCMC->ln_likelihood_last,
        MCMC->ln_likelihood_best);
     fprintf(fp_report_props,"ln(probability)=%13.6le\n",MCMC->ln_Pr_new);
     fprintf(fp_report_props,"n_success      =%09d\n",   MCMC->n_success);
     fprintf(fp_report_props,"n_fail         =%09d\n",   MCMC->n_fail);
     fprintf(fp_report_props,"success rate   =%.2f %%\n",(double)MCMC->n_success/(double)(MCMC->n_fail+MCMC->n_success)*100.0);
     fclose(fp_report_props);
  }

  return(flag_success);
}

