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

int generate_MCMC_proposition(MCMC_info *MCMC){
  int              flag_success;
  int              i_P,i_DS;
  char             format_string[64];
  static int       n_P;
  static double   *P_new;
  static double  **M_new;
  static double   *P_last;
  static double   *P_chain;
  static double  **M_last;
  static RNG_info *RNG;
  static int       flag_report_props;
  static int       n_DS;
  static int      *n_M;

  // Initialize a few things on the first call
  switch(MCMC->first_gen_proposition_call){
    case TRUE:
      n_P        =MCMC->n_P;
      P_new      =MCMC->P_new;
      M_new      =MCMC->M_new;
      P_last     =MCMC->P_last;
      M_last     =MCMC->M_last;
      P_chain    =MCMC->P_chain;
      RNG        =MCMC->RNG;
      n_DS       =MCMC->n_DS;
      n_M        =MCMC->n_M;
      flag_report_props=check_mode_for_flag(MCMC->mode,MCMC_REPORT_PROPS);
      MCMC->first_gen_proposition_call=FALSE;
      break;
  }

  // We need to compute the likelihood of the initial conditions
  //   if we are just starting a chain
  if(MCMC->flag_init_chain){
    memcpy(P_new,  MCMC->P_init,(size_t)n_P*sizeof(double));
    memcpy(P_last, MCMC->P_init,(size_t)n_P*sizeof(double));
    memcpy(P_chain,MCMC->P_init,(size_t)n_P*sizeof(double));
    while(MCMC->map_P_to_M(P_new,MCMC,M_new)){
      MCMC->first_map_call=FALSE;
      MCMC->n_map_calls++;
      generate_MCMC_parameters(MCMC);
    }
    MCMC->first_map_call=FALSE;
    MCMC->n_map_calls++;
    MCMC->compute_MCMC_ln_likelihood(MCMC,M_new,P_new,&(MCMC->ln_likelihood_new));
    MCMC->ln_likelihood_chain  =MCMC->ln_likelihood_new;
    MCMC->ln_likelihood_last   =MCMC->ln_likelihood_new;
    MCMC->ln_Pr_chain          =0.;
    MCMC->ln_Pr_new            =0.;
    MCMC->first_likelihood_call=FALSE;
    MCMC->flag_init_chain      =FALSE;
  }

  // Set the last new state to the new last state
  MCMC->ln_likelihood_last=MCMC->ln_likelihood_new;
  MCMC->ln_Pr_last        =MCMC->ln_Pr_new;
  memcpy(P_last,P_new,n_P*sizeof(double));
  for(i_DS=0;i_DS<n_DS;i_DS++)
    memcpy(M_last[i_DS],M_new[i_DS],n_M[i_DS]*sizeof(double));

  // Keep generating parameter sets until the mapping function is satisfied
  generate_MCMC_parameters(MCMC);
  while(MCMC->map_P_to_M(P_new,MCMC,M_new)){
    MCMC->n_map_calls++;
    MCMC->first_map_call=FALSE;
    generate_MCMC_parameters(MCMC);
  }
  MCMC->first_map_call=FALSE;
  MCMC->n_map_calls++;
  MCMC->compute_MCMC_ln_likelihood(MCMC,M_new,P_new,&(MCMC->ln_likelihood_new));
  MCMC->first_likelihood_call=FALSE;

  // Decide if this is a successful proposition or not...
  if(MCMC->ln_likelihood_new>MCMC->ln_likelihood_chain){
    MCMC->ln_Pr_new=0.;
    flag_success=TRUE;
  }
  else{
    MCMC->ln_Pr_new=MCMC->ln_likelihood_new-MCMC->ln_likelihood_chain;
    if((double)random_number(RNG)<=exp(MCMC->ln_Pr_new))
      flag_success=TRUE;
    else
      flag_success=FALSE;
  }
fprintf(stderr,"%13.6le %13.6le %13.6le %d\n",MCMC->ln_likelihood_new,MCMC->ln_likelihood_chain,MCMC->ln_Pr_new,flag_success);
if(MCMC->n_map_calls>20) exit(1);
  // ... if it is, then update the chain ...
  if(flag_success){
    memcpy(P_chain,P_new,(size_t)n_P*sizeof(double));
    MCMC->ln_likelihood_chain=MCMC->ln_likelihood_new;
    MCMC->ln_Pr_chain        =MCMC->ln_Pr_new;
  }

  // If we're integrating, then count proposals
  if(MCMC->flag_integrate_on){
    if(flag_success)
      MCMC->n_success++;
    else
      MCMC->n_fail++;
    MCMC->i_proposal++;
  }

  // Report the proposition if asked to
  if(flag_report_props && MCMC->my_chain==SID.My_rank){
    if(MCMC->flag_integrate_on){
      if(flag_success)
        SID_log("Proposal #%09d: SUCCEEDED",SID_LOG_ALLRANKS|SID_LOG_OPEN,MCMC->i_proposal);
      else
        SID_log("Proposal #%09d: REJECTED",SID_LOG_ALLRANKS|SID_LOG_OPEN,MCMC->i_proposal);
    }
    sprintf(format_string,"%s = %%13.6le (was %%13.6le)",MCMC->P_name_format);
    for(i_P=0;i_P<n_P;i_P++)
      SID_log(format_string,SID_LOG_ALLRANKS|SID_LOG_COMMENT,MCMC->P_names[i_P],P_new[i_P],P_last[i_P]);
    SID_log("ln(likelihood) =%le+constant (was %le+constant)",SID_LOG_ALLRANKS|SID_LOG_COMMENT,MCMC->ln_likelihood_new,MCMC->ln_likelihood_last);
    SID_log("ln(probability)=%le", SID_LOG_ALLRANKS|SID_LOG_COMMENT,MCMC->ln_Pr_new);
    if(MCMC->flag_integrate_on){
      SID_log("n_success      =%09d",SID_LOG_ALLRANKS|SID_LOG_COMMENT,MCMC->n_success);
      SID_log("n_fail         =%09d",SID_LOG_ALLRANKS|SID_LOG_COMMENT,MCMC->n_fail);
      SID_log("success rate   =%.2f %%",SID_LOG_ALLRANKS|SID_LOG_COMMENT,(double)MCMC->n_success/(double)(MCMC->n_fail+MCMC->n_success)*100.0);
    }
    SID_log("",SID_LOG_ALLRANKS|SID_LOG_NOPRINT|SID_LOG_CLOSE);
  }

  return(flag_success);
}

