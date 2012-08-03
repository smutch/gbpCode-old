#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>

void free_pspec(pspec_info *pspec){
  SID_log("Freeing power spectrum...",SID_LOG_OPEN);
  free_cosmo(&(pspec->cosmo));
  free_field(&(pspec->FFT));
  SID_free(SID_FARG (pspec->k_1D));
  SID_free(SID_FARG (pspec->n_modes_1D));
  SID_free(SID_FARG (pspec->n_modes_2D));
  int i_run;
  for(i_run=0;i_run<4;i_run++){
      SID_free(SID_FARG (pspec->P_k_1D[i_run]));
      SID_free(SID_FARG (pspec->dP_k_1D[i_run]));
      SID_free(SID_FARG (pspec->P_k_2D[i_run]));
      SID_free(SID_FARG (pspec->dP_k_2D[i_run]));
  }
  SID_free(SID_FARG pspec->P_k_1D);
  SID_free(SID_FARG pspec->dP_k_1D);
  SID_free(SID_FARG pspec->P_k_2D);
  SID_free(SID_FARG pspec->dP_k_2D);
  SID_log("Done.",SID_LOG_CLOSE);
}

