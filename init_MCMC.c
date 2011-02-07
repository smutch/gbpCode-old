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

void init_MCMC(MCMC_info *MCMC,const char *problem_name,void *params,int (*f)(double *,MCMC_info *,double **),int n_P,double *P_init,char **P_names,double *P_limit_min,double *P_limit_max,int n_arrays,...){
  int     i_P;
  int     i_array;
  va_list vargs;
  va_start(vargs,n_arrays);

  SID_log("Initializing MCMC structure...",SID_LOG_OPEN);

  // Set defaults to bare minimums
  sprintf(MCMC->filename_output_dir,"./%s_MCMC/",SID.My_binary);
  MCMC->n_avg                   =100;
  MCMC->n_iterations_burn       =4;
  MCMC->n_iterations            =8;
  MCMC->n_thin                  =1;
  MCMC->coverage_size           =10;
  MCMC->flag_autocor_on         =FALSE;
  MCMC->flag_integrate_on       =FALSE;
  MCMC->flag_analysis_on        =TRUE;
  MCMC->first_map_call          =TRUE;
  MCMC->first_link_call         =TRUE;
  MCMC->flag_init_chain         =TRUE;
  MCMC->first_gen_parameter_call  =TRUE; 
  MCMC->first_gen_proposition_call=TRUE; 
  MCMC->first_likelihood_call     =TRUE; 
  MCMC->ln_likelihood_last        =0.; 
  MCMC->ln_likelihood_new         =0.; 
  MCMC->ln_likelihood_chain       =0.;
  MCMC->P_init                  =NULL;
  MCMC->P_new                   =NULL;
  MCMC->P_last                  =NULL;
  MCMC->P_chain                 =NULL;
  MCMC->P_limit_min             =NULL;
  MCMC->P_limit_max             =NULL;
  MCMC->n_M                     =NULL;
  MCMC->M_new                   =NULL;
  MCMC->M_last                  =NULL;
  MCMC->DS                      =NULL;
  MCMC->last                    =NULL;
  MCMC->V                       =NULL;
  MCMC->m                       =NULL;
  MCMC->b                       =NULL;
  MCMC->RNG                     =NULL;
  MCMC->params                  =NULL;
  MCMC->temperature             =1.0;
  MCMC->n_DS                    =0;
  MCMC->n_M_total               =0;
  MCMC->n_fail                  =0;
  MCMC->n_success               =0;
  MCMC->i_proposal              =0;
  MCMC->n_map_calls             =0;

  // Process the passed arguments
  MCMC->map_P_to_M                =f;
  MCMC->compute_MCMC_ln_likelihood=compute_MCMC_ln_likelihood_default;
  MCMC->params                    =params;
  MCMC->n_P                       =n_P;
  MCMC->problem_name=(char *)SID_malloc(sizeof(char)*MCMC_NAME_SIZE);
  sprintf(MCMC->problem_name,"%s\0",problem_name);
  MCMC->P_names      =(char **)SID_malloc(sizeof(char *)*MCMC->n_P);
  MCMC->P_name_length=0;
  for(i_P=0;i_P<n_P;i_P++){
    MCMC->P_names[i_P]=(char *)SID_malloc(sizeof(char)*MCMC_NAME_SIZE);
    sprintf(MCMC->P_names[i_P],"%s\0",P_names[i_P]);
    MCMC->P_name_length=MAX(MCMC->P_name_length,strlen(MCMC->P_names[i_P]));
  }
  sprintf(MCMC->P_name_format,"%%-%ds",MCMC->P_name_length);

  // Initialize the covariance matrix
  set_MCMC_covariance(MCMC,NULL);

  // Initialize the MCMC mode and set things associated with it
  set_MCMC_mode(MCMC,MCMC_MODE_DEFAULT); // MCMC->my_chain is set here

  // Set parameter arrays and limits
  MCMC->P_init       =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
  MCMC->P_new        =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
  MCMC->P_last       =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
  MCMC->P_chain      =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
  MCMC->P_limit_min  =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
  MCMC->P_limit_max  =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
  if(P_limit_min==NULL){
    for(i_P=0;i_P<n_P;i_P++)
      MCMC->P_limit_min[i_P]=DBL_MIN*1e3;
  }
  else{
    for(i_P=0;i_P<n_P;i_P++)
      MCMC->P_limit_min[i_P]=P_limit_min[i_P];
  }
  if(P_limit_max==NULL){
    for(i_P=0;i_P<n_P;i_P++)
      MCMC->P_limit_max[i_P]=DBL_MAX*1e-3;
  }
  else{
    for(i_P=0;i_P<n_P;i_P++)
      MCMC->P_limit_max[i_P]=P_limit_max[i_P];
  }

  // Set parameter initial values
  memcpy(MCMC->P_init, P_init,(size_t)MCMC->n_P*sizeof(double));
  memcpy(MCMC->P_new,  P_init,(size_t)MCMC->n_P*sizeof(double));
  memcpy(MCMC->P_last, P_init,(size_t)MCMC->n_P*sizeof(double));
  memcpy(MCMC->P_chain,P_init,(size_t)MCMC->n_P*sizeof(double));

  // Set arrays
  MCMC->n_arrays=n_arrays;
  if(n_arrays>0){
    MCMC->array     =(double **)SID_malloc(sizeof(double *)*MCMC->n_arrays);
    MCMC->array_name=(char   **)SID_malloc(sizeof(char   *)*MCMC->n_arrays);
    for(i_array=0;i_array<n_arrays;i_array++){
      MCMC->array[i_array]     =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
      MCMC->array_name[i_array]=(char *)SID_malloc(sizeof(char)*MCMC_NAME_SIZE);
      memcpy(MCMC->array[i_array],(double *)va_arg(vargs,double *),(size_t)(MCMC->n_P)*sizeof(double));
      sprintf(MCMC->array_name[i_array],"%s\0",(char *)va_arg(vargs,char *));
    }
  }
  else
    MCMC->array=NULL;

  // Set autotune defaults (if needed)
  if(check_mode_for_flag(MCMC->mode,MCMC_MODE_AUTOTUNE))
    set_MCMC_autotune(MCMC,-1.,-1.,-1.,-1,-1,-1); // Negatives mean use defaults (set in gbpMCMC.h)

  // Initilize the random number generator
  MCMC->RNG=(RNG_info *)SID_malloc(sizeof(RNG_info));
  init_seed_from_clock(&(MCMC->seed));
  init_RNG(&(MCMC->seed),MCMC->RNG,RNG_DEFAULT);

  SID_log("Done.",SID_LOG_CLOSE);

  va_end(vargs);
}

