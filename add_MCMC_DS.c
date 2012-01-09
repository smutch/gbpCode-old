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

  // Create new item
  new_DS=(MCMC_DS_info *)SID_malloc(sizeof(MCMC_DS_info));
  strcpy(new_DS->name,name);
  if(n_D>0){
    for(i_D=0,n_M=1;i_D<n_D;i_D++)
      n_M*=D[i_D];
  }
  else
    n_M=0;
  new_DS->n_D          =n_D;
  new_DS->D            =(int *)SID_malloc(sizeof(int)*n_D);
  new_DS->n_M          =n_M;
  new_DS->M_target     =(double *)SID_malloc(sizeof(double)*new_DS->n_M);
  new_DS->dM_target    =(double *)SID_malloc(sizeof(double)*new_DS->n_M);
  memcpy(new_DS->D,          D,(size_t)new_DS->n_D*sizeof(int));
  memcpy(new_DS->M_target,  DS,(size_t)new_DS->n_M*sizeof(double));
  memcpy(new_DS->dM_target,dDS,(size_t)new_DS->n_M*sizeof(double));
  new_DS->params   =params;
  new_DS->n_arrays =n_arrays;
  if(n_arrays>0){
    new_DS->array     =(double **)SID_malloc(sizeof(double *)*new_DS->n_arrays);
    new_DS->array_name=(char   **)SID_malloc(sizeof(char   *)*new_DS->n_arrays);
    for(i_array=0;i_array<n_arrays;i_array++){
      new_DS->array[i_array]     =(double *)SID_malloc(sizeof(double)*new_DS->n_M);
      new_DS->array_name[i_array]=(char *)SID_malloc(sizeof(char)*MCMC_NAME_SIZE);
      memcpy(new_DS->array[i_array],(double *)va_arg(vargs,double *),(size_t)(new_DS->n_M)*sizeof(double));
      strcpy(new_DS->array_name[i_array],(char *)va_arg(vargs,char *));
    }
  }
  else
    MCMC->array=NULL;
  new_DS->next=NULL;
  
  // Add new item to linked list
  if(MCMC->DS==NULL)
    MCMC->DS=new_DS;
  else
    MCMC->last->next=new_DS;
  MCMC->last=new_DS;

  MCMC->n_DS++;

  SID_log("Done.",SID_LOG_CLOSE);

  va_end(vargs);
}

