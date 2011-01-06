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

void init_MCMC(MCMC_info *MCMC,const char *problem_name,void *params,int (*f)(double *,MCMC_info *,double **),int n_P,double *P_init,char **P_names,int *flag_limit,double *P_limit_min,double *P_limit_max,int n_arrays,...){
  int     i_P;
  int     i_array;
  va_list vargs;
  va_start(vargs,n_arrays);

  SID_log("Initializing MCMC structure...",SID_LOG_OPEN);

  MCMC->map_P_to_M                =f;
  MCMC->compute_MCMC_ln_likelihood=compute_MCMC_ln_likelihood_default;
  MCMC->params                    =params;
  MCMC->temperature               =0.75;
  MCMC->n_P                       =n_P;
  MCMC->n_DS                      =0;
  MCMC->n_M_total                 =0;
  MCMC->problem_name=(char *)SID_malloc(sizeof(char)*MCMC_NAME_SIZE);
  sprintf(MCMC->problem_name,"%s\0",problem_name);
  MCMC->P_names     =(char **)SID_malloc(sizeof(char *)*MCMC->n_P);
  for(i_P=0;i_P<n_P;i_P++){
    MCMC->P_names[i_P]=(char *)SID_malloc(sizeof(char)*MCMC_NAME_SIZE);
    sprintf(MCMC->P_names[i_P],"%s\0",P_names[i_P]);
  }
  MCMC->P_init       =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
  MCMC->flag_limit   =(int    *)SID_malloc(sizeof(int)*MCMC->n_P);
  MCMC->P_limit_min  =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
  MCMC->P_limit_max  =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
  if(flag_limit!=NULL){
    for(i_P=0;i_P<n_P;i_P++){
      MCMC->flag_limit[i_P]=flag_limit[i_P];
      if(flag_limit[i_P]){
        MCMC->P_limit_min[i_P]=P_limit_min[i_P];
        MCMC->P_limit_max[i_P]=P_limit_max[i_P];
      }  
    }
  }
  memcpy(MCMC->P_init,P_init,(size_t)MCMC->n_P*sizeof(double));
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

  // Set defaults to bare minimums
  MCMC->n_avg                 =MAX(10,n_P*n_P);
  MCMC->n_avg_covariance      =MAX(10,n_P*n_P);
  MCMC->n_iterations_burn     =4;
  MCMC->n_iterations          =8;
  MCMC->n_thin                =1;
  MCMC->coverage_size         =10;
  MCMC->flag_autocor_on       =FALSE;
  MCMC->flag_integrate_on     =TRUE;
  MCMC->flag_analysis_on      =TRUE;
  MCMC->first_map_call        =TRUE;
  MCMC->mode                  =MCMC_MODE_SERIAL;
  sprintf(MCMC->filename_output_dir,"./%s_MCMC/",SID.My_binary);
  MCMC->V                     =NULL;
  MCMC->m                     =NULL;
  
  MCMC->DS           =NULL;
  MCMC->last         =NULL;

  add_covariance_to_MCMC(MCMC,NULL);

  SID_log("Done.",SID_LOG_CLOSE);

  va_end(vargs);
}

