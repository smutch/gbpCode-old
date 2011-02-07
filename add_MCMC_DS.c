#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMCMC.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>

void add_MCMC_DS(MCMC_info *MCMC,const char *name,int n_M,double *DS,double *dDS,void *params,int n_arrays,...){
  int           i_array;
  MCMC_DS_info *new;
  va_list       vargs;
  va_start(vargs,n_arrays);

  SID_log("Adding data set {%s} to MCMC structure...",SID_LOG_OPEN,name);

  // Create new item
  new=(MCMC_DS_info *)SID_malloc(sizeof(MCMC_DS_info));
  strcpy(new->name,name);
  new->n_M          =n_M;
  new->M_target     =(double *)SID_malloc(sizeof(double)*new->n_M);
  new->dM_target    =(double *)SID_malloc(sizeof(double)*new->n_M);
  memcpy(new->M_target,  DS,(size_t)new->n_M*sizeof(double));
  memcpy(new->dM_target,dDS,(size_t)new->n_M*sizeof(double));
  new->params   =params;
  new->n_arrays =n_arrays;
  if(n_arrays>0){
    new->array     =(double **)SID_malloc(sizeof(double *)*new->n_arrays);
    new->array_name=(char   **)SID_malloc(sizeof(char   *)*new->n_arrays);
    for(i_array=0;i_array<n_arrays;i_array++){
      new->array[i_array]     =(double *)SID_malloc(sizeof(double)*new->n_M);
      new->array_name[i_array]=(char *)SID_malloc(sizeof(char)*MCMC_NAME_SIZE);
      memcpy(new->array[i_array],(double *)va_arg(vargs,double *),(size_t)(new->n_M)*sizeof(double));
      strcpy(new->array_name[i_array],(char *)va_arg(vargs,char *));
    }
  }
  else
    MCMC->array=NULL;
  new->next=NULL;
  
  // Add new item to linked list
  if(MCMC->DS==NULL)
    MCMC->DS=new;
  else
    MCMC->last->next=new;
  MCMC->last=new;

  MCMC->n_DS++;

  SID_log("Done.",SID_LOG_CLOSE);

  va_end(vargs);
}

