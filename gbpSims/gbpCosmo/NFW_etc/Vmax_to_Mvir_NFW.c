#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>
#include <gbpCosmo_NFW_etc.h>
#include <gsl/gsl_sf_expint.h>

void init_Vmax_to_Mvir_NFW(cosmo_info **cosmo,
                           int          mode,
                           double       z){
  interp_info *interp;
  if(!ADaPS_exist((*cosmo),"lVmax_to_lMvir_%.5f_interp",z)){
    SID_log("Initializing Vmax->M_vir interpolation...",SID_LOG_OPEN);
    int     n_k;
    double *lk_P;
    double *lM;
    double *lVmax;
    int     n_M  =201;
    double  lM_lo= 0.;
    double  lM_hi=20.;
    double  dlM  =(lM_hi-lM_lo)/(double)(n_M-1);
    lM     =(double *)SID_malloc(sizeof(double)*n_M);
    lVmax  =(double *)SID_malloc(sizeof(double)*n_M);
    int i_M;
    for(i_M=0;i_M<n_M;i_M++){
      if(i_M==0)            lM[i_M]=lM_lo;
      else if(i_M==(n_M-1)) lM[i_M]=lM_hi;
      else                  lM[i_M]=lM[i_M-1]+dlM;
    }
    for(i_M=0;i_M<n_M;i_M++){
      lM[i_M]+=take_log10(M_SOL);
      lVmax[i_M]=take_log10(V_max_NFW(take_alog10(lM[i_M]),z,mode,cosmo));
    }
    init_interpolate(lVmax,
                     lM,
                     (size_t)n_M,
                     gsl_interp_cspline,
                     &interp);
    ADaPS_store_interp(cosmo,
                       (void *)(interp),
                       "lVmax_to_lMvir_%.5f_interp",z);
    SID_free(SID_FARG lM);
    SID_free(SID_FARG lVmax);
    SID_log("Done.",SID_LOG_CLOSE);
  }
}
double Vmax_to_Mvir_NFW(double       V_max,
                        double       z,
                        int          mode,
                        cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double V2_vir;
  double g_c;
  double r_val=0.;
  if(V_max>0.){
    interp_info *interp;
    init_Vmax_to_Mvir_NFW(cosmo,mode,z);
    interp=(interp_info *)ADaPS_fetch(*cosmo,"lVmax_to_lMvir_%.5f_interp",z);
    r_val =take_alog10(interpolate(interp,take_log10(V_max)));
  }
  return(r_val);
}

