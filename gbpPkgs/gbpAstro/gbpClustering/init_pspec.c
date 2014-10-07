#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>

void init_pspec(pspec_info *pspec,
	        int mass_assignment_scheme,
                double redshift,double box_size,int grid_size,
                double k_min_1D,double k_max_1D,double dk_1D,
                double k_min_2D,double k_max_2D,double dk_2D){

  SID_log("Initializing power spectrum...",SID_LOG_OPEN);

  // Initialize flags
  pspec->mass_assignment_scheme=mass_assignment_scheme;
  pspec->initialized           =TRUE;

  // Initialize constants
  pspec->redshift =redshift;
  pspec->box_size =box_size;
  pspec->grid_size=grid_size;
  pspec->k_min_1D =k_min_1D;
  pspec->k_max_1D =k_max_1D;
  pspec->dk_1D    =dk_1D;
  pspec->k_min_2D =k_min_2D;
  pspec->k_max_2D =k_max_2D;
  pspec->dk_2D    =dk_2D;
  pspec->n_k_1D   =(int)(0.5+(pspec->k_max_1D-pspec->k_min_1D)/pspec->dk_1D);
  pspec->n_k_2D   =(int)(0.5+(pspec->k_max_2D-pspec->k_min_2D)/pspec->dk_2D);

  // Initialize arrays
  pspec->k_1D      =(double  *)SID_malloc(sizeof(double)*(pspec->n_k_1D)); 
  pspec->n_modes_1D=(int     *)SID_malloc(sizeof(int   )*(pspec->n_k_1D)); 
  pspec->n_modes_2D=(int     *)SID_malloc(sizeof(int   )*(pspec->n_k_2D)*(pspec->n_k_2D));
  pspec->P_k_1D    =(double **)SID_calloc(4*sizeof(double *)); 
  pspec->dP_k_1D   =(double **)SID_calloc(4*sizeof(double *)); 
  pspec->P_k_2D    =(double **)SID_calloc(4*sizeof(double *));
  pspec->dP_k_2D   =(double **)SID_calloc(4*sizeof(double *));
  int i;
  for(i=0;i<4;i++){
      pspec->flag_processed[i]=FALSE;
      pspec->P_k_1D[i] =(double *)SID_calloc(sizeof(double)*(pspec->n_k_1D)); 
      pspec->dP_k_1D[i]=(double *)SID_calloc(sizeof(double)*(pspec->n_k_1D)); 
      pspec->P_k_2D[i] =(double *)SID_calloc(sizeof(double)*(pspec->n_k_2D)*(pspec->n_k_2D));
      pspec->dP_k_2D[i]=(double *)SID_calloc(sizeof(double)*(pspec->n_k_2D)*(pspec->n_k_2D));
  }

  // Initialize the grid
  double L[3];
  int    n[3];
  n[0]=grid_size;
  n[1]=n[0];
  n[2]=n[0];
  L[0]=box_size;
  L[1]=L[0];
  L[2]=L[0];
  init_field(3,n,L,&(pspec->FFT));

  // Initialize the cosmology
  init_cosmo_std(&(pspec->cosmo));

  SID_log("Done.",SID_LOG_CLOSE);
}

