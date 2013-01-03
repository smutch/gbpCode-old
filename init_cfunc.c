#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>

void init_cfunc(cfunc_info *cfunc,int    n_data,  int    F_random,int    n_bits_PHK,
                double redshift,  double box_size,int    n_jack,
                double r_min_l1D, double r_max_1D,double dr_1D,
                double r_min_2D,  double r_max_2D,double dr_2D){
  SID_log("Initializing correlation function...",SID_LOG_OPEN);

  // Initialize flags
  cfunc->initialized    =TRUE;
  cfunc->flag_compute_RR=TRUE;

  // Initialize constants
  cfunc->n_data    =n_data;
  cfunc->n_random  =n_data*F_random;
  cfunc->F_random  =F_random;
  cfunc->redshift  =redshift;
  cfunc->box_size  =box_size;
  cfunc->n_jack    =n_jack;
  cfunc->r_min_l1D =r_min_l1D;
  cfunc->lr_min_l1D=take_log10(r_min_l1D);
  cfunc->r_max_1D  =r_max_1D;
  cfunc->r_min_2D  =r_min_2D;
  cfunc->r_max_2D  =r_max_2D;
  cfunc->r_max     =MAX(cfunc->r_max_1D,cfunc->r_max_2D);
  cfunc->dr_1D     =dr_1D;
  cfunc->dr_2D     =dr_2D;

  // Decide on PHK boundary widths
  cfunc->n_bits_PHK=n_bits_PHK;
  for(cfunc->PHK_width=1;cfunc->PHK_width<20 && (double)cfunc->PHK_width*(cfunc->box_size/pow(2.,(double)(cfunc->n_bits_PHK)))<cfunc->r_max;) cfunc->PHK_width++;
  SID_log("using %d-bit keys and %d-key boundries...",SID_LOG_CONTINUE,cfunc->n_bits_PHK,cfunc->PHK_width);

/*
  // Initialize the bit size of the PHKs
  int i_shift    =PHK_width;
  int n_bits_PHK;
  cfunc->PHK_width=PHK_width;
  for(cfunc->n_bits_PHK=1;(cfunc->box_size/pow(2.,(double)(cfunc->n_bits_PHK+1)))>cfunc->r_max && cfunc->n_bits_PHK<=20;) cfunc->n_bits_PHK++;

  // Add bits to allow for the boundary width
  while(i_shift>1){
     i_shift/=2;
     cfunc->n_bits_PHK++;
  }
*/

  // Initialize the number of bins
  cfunc->n_1D        =(int)(0.5+(cfunc->r_max_1D)/cfunc->dr_1D);
  cfunc->n_2D        =(int)(0.5+(cfunc->r_max_2D-cfunc->r_min_2D)/cfunc->dr_2D);
  cfunc->n_2D_total  =cfunc->n_2D*cfunc->n_2D;
  cfunc->n_jack_total=cfunc->n_jack*cfunc->n_jack*cfunc->n_jack;

  // Initialize logarythmic bin sizes
  cfunc->dr_l1D=(take_log10(cfunc->r_max_1D)-cfunc->lr_min_l1D)/(double)cfunc->n_1D;

  // Initialize arrays
  cfunc->CFUNC_l1D =(double **)SID_malloc(sizeof(double *)*4);
  cfunc->dCFUNC_l1D=(double **)SID_malloc(sizeof(double *)*4);
  cfunc->COVMTX_l1D=(double **)SID_malloc(sizeof(double *)*4);
  cfunc->CFUNC_1D  =(double **)SID_malloc(sizeof(double *)*4);
  cfunc->dCFUNC_1D =(double **)SID_malloc(sizeof(double *)*4);
  cfunc->COVMTX_1D =(double **)SID_malloc(sizeof(double *)*4);
  cfunc->CFUNC_2D  =(double **)SID_malloc(sizeof(double *)*4);
  cfunc->dCFUNC_2D =(double **)SID_malloc(sizeof(double *)*4);
  cfunc->COVMTX_2D =(double **)SID_malloc(sizeof(double *)*4);
  int i_run;
  for(i_run=0;i_run<4;i_run++){
     cfunc->CFUNC_l1D[i_run] =(double *)SID_malloc(sizeof(double)*(cfunc->n_1D));
     cfunc->dCFUNC_l1D[i_run]=(double *)SID_malloc(sizeof(double)*(cfunc->n_1D));
     cfunc->COVMTX_l1D[i_run]=(double *)SID_malloc(sizeof(double)*(cfunc->n_1D*cfunc->n_1D));
     cfunc->CFUNC_1D[i_run]  =(double *)SID_malloc(sizeof(double)*(cfunc->n_1D));
     cfunc->dCFUNC_1D[i_run] =(double *)SID_malloc(sizeof(double)*(cfunc->n_1D));
     cfunc->COVMTX_1D[i_run] =(double *)SID_malloc(sizeof(double)*(cfunc->n_1D*cfunc->n_1D));
     cfunc->CFUNC_2D[i_run]  =(double *)SID_malloc(sizeof(double)*(cfunc->n_2D)*(cfunc->n_2D));
     cfunc->dCFUNC_2D[i_run] =(double *)SID_malloc(sizeof(double)*(cfunc->n_2D)*(cfunc->n_2D));
     cfunc->COVMTX_2D[i_run] =(double *)SID_malloc(sizeof(double)*(cfunc->n_2D_total*cfunc->n_2D_total));
  }

  cfunc->DD_l1D    =(long long **)SID_malloc(sizeof(long long *)*(cfunc->n_jack_total+1));
  cfunc->DR_l1D    =(long long **)SID_malloc(sizeof(long long *)*(cfunc->n_jack_total+1));
  cfunc->RR_l1D    =(long long **)SID_malloc(sizeof(long long *)*(cfunc->n_jack_total+1));
  cfunc->DD_1D     =(long long **)SID_malloc(sizeof(long long *)*(cfunc->n_jack_total+1));
  cfunc->DR_1D     =(long long **)SID_malloc(sizeof(long long *)*(cfunc->n_jack_total+1));
  cfunc->RR_1D     =(long long **)SID_malloc(sizeof(long long *)*(cfunc->n_jack_total+1));
  cfunc->DD_2D     =(long long **)SID_malloc(sizeof(long long *)*(cfunc->n_jack_total+1));
  cfunc->DR_2D     =(long long **)SID_malloc(sizeof(long long *)*(cfunc->n_jack_total+1));
  cfunc->RR_2D     =(long long **)SID_malloc(sizeof(long long *)*(cfunc->n_jack_total+1));
  int i_jack;
  for(i_jack=0;i_jack<=cfunc->n_jack_total;i_jack++){
    cfunc->DD_l1D[i_jack]=(long long *)SID_calloc(sizeof(long long)*cfunc->n_1D);
    cfunc->DR_l1D[i_jack]=(long long *)SID_calloc(sizeof(long long)*cfunc->n_1D);
    cfunc->RR_l1D[i_jack]=(long long *)SID_calloc(sizeof(long long)*cfunc->n_1D);
    cfunc->DD_1D[i_jack] =(long long *)SID_calloc(sizeof(long long)*cfunc->n_1D);
    cfunc->DR_1D[i_jack] =(long long *)SID_calloc(sizeof(long long)*cfunc->n_1D);
    cfunc->RR_1D[i_jack] =(long long *)SID_calloc(sizeof(long long)*cfunc->n_1D);
    cfunc->DD_2D[i_jack] =(long long *)SID_calloc(sizeof(long long)*cfunc->n_2D*cfunc->n_2D);
    cfunc->DR_2D[i_jack] =(long long *)SID_calloc(sizeof(long long)*cfunc->n_2D*cfunc->n_2D);
    cfunc->RR_2D[i_jack] =(long long *)SID_calloc(sizeof(long long)*cfunc->n_2D*cfunc->n_2D);
  }

  // Initialize the cosmology
  init_cosmo_std(&(cfunc->cosmo));

  SID_log("Done.",SID_LOG_CLOSE);
}

