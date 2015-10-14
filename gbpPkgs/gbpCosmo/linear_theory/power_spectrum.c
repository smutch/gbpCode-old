#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

// Initialize the power spectrum at z=0 
void init_power_spectrum_TF(cosmo_info **cosmo){
  char   *line=NULL;
  size_t  line_length=0;
  int     n_k_tmp;
  double  k_P;
  double *lk_P;
  double *lk_P_tmp;
  double *lP_k;
  double *lP_k_dark;
  double *lP_k_gas;
  double *lP_k_tmp;
  double *lP_k_gas_tmp;
  double *lP_k_dark_tmp;
  interp_info *interp;
  interp_info *interp_dark;
  interp_info *interp_gas;
  double  M_8;
  double  sigma2_unnorm;
  double  Omega_b;
  double  Omega_M;
  double  h_Hubble;
  double  norm;
  double  n_spectral;
  double  M_WDM,R_WDM;

  SID_log("Initializing P(k)...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Fetch the transfer function filename (must be set before power_spectrum() is called).
  char filename_TF[MAX_FILENAME_LENGTH];
  if(!ADaPS_exist(*cosmo,"filename_transfer_function"))
     SID_trap_error("Transfer function filename has not been specified prior to calling init_power_spectrum_TF().",ERROR_LOGIC);
  memcpy(filename_TF,ADaPS_fetch((*cosmo),"filename_transfer_function"),MAX_FILENAME_LENGTH*sizeof(char));
  
  n_spectral =((double *)ADaPS_fetch((*cosmo),"n_spectral"))[0];
  h_Hubble   =((double *)ADaPS_fetch((*cosmo),"h_Hubble"))[0];
  Omega_M    =((double *)ADaPS_fetch((*cosmo),"Omega_M"))[0];
  Omega_b    =((double *)ADaPS_fetch((*cosmo),"Omega_b"))[0];

  // Read the transfer function
  int   n_skip=5;
  FILE *fp    =fopen(filename_TF,"r");
  int   n_k_in=count_lines_data(fp);
  int   n_k   =n_k_in/n_skip+1;
  lk_P      =(double *)SID_malloc(sizeof(double)*n_k);
  lP_k      =(double *)SID_malloc(sizeof(double)*n_k);
  lP_k_gas  =(double *)SID_malloc(sizeof(double)*n_k);
  lP_k_dark =(double *)SID_malloc(sizeof(double)*n_k);
  for(int i=0,j=0;j<n_k_in;j++){
    grab_next_line_data(fp,&line,&line_length);
    if(!(j%n_skip)){
       grab_double(line,1,&(lk_P[i]));
       grab_double(line,2,&(lP_k_dark[i]));
       grab_double(line,3,&(lP_k_gas[i]));
       lP_k[i]=((Omega_M-Omega_b)*lP_k_dark[i]+Omega_b*lP_k_gas[i])/Omega_M;
       i++;
    }
    n_k=i;
  }
  fclose(fp);

  // Take the log of lk_P
  for(int i=0;i<n_k;i++)
    lk_P[i]=take_log10(lk_P[i]/(M_PER_MPC/h_Hubble));

  // Supress small-scale power following the
  //    ENS recipe if M_ENS or M_WDM is set > 0
  /*
  M_WDM=-1.;
  if(ADaPS_exist(*cosmo,"M_ENS"))
    M_WDM=((double *)ADaPS_fetch(*cosmo,"M_ENS"))[0];
  else if(ADaPS_exist(*cosmo,"M_WDM"))
    M_WDM=((double *)ADaPS_fetch(*cosmo,"M_WDM"))[0];
  if(M_WDM>0.){
    R_WDM=0.065*pow((M_WDM/(1e11*M_SOL))/Omega_M,ONE_THIRD)*M_PER_MPC;
    for(i=0;i<n_k;i++){
      k_P=take_alog10(lk_P[i]);
      if(R_WDM/R_of_k(k_P)>8.){
         lP_k[i]     =0.;
         lP_k_gas[i] =0.;
         lP_k_dark[i]=0.;
      }
      else{
         lP_k[i]     *=exp(-0.5*k_P*R_WDM-0.5*pow(k_P*R_WDM,2.));
         lP_k_gas[i] *=exp(-0.5*k_P*R_WDM-0.5*pow(k_P*R_WDM,2.));
         lP_k_dark[i]*=exp(-0.5*k_P*R_WDM-0.5*pow(k_P*R_WDM,2.));
      }
    }
  }
  */

  // Compute power spectrum
  for(int i=0;i<n_k;i++){
    lP_k[i]     =2.*take_log10(lP_k[i])     +(n_spectral)*lk_P[i];
    lP_k_gas[i] =2.*take_log10(lP_k_gas[i]) +(n_spectral)*lk_P[i];
    lP_k_dark[i]=2.*take_log10(lP_k_dark[i])+(n_spectral)*lk_P[i];
  }

  // Normalize at large scales
  norm=lP_k[0];
  for(int i=0;i<n_k;i++){
    lP_k[i]     -=norm;
    lP_k_gas[i] -=norm;
    lP_k_dark[i]-=norm;
  }

  // Store P(k) arrays now so that the normalization routine can use them
  ADaPS_store(cosmo,
              (void *)(&n_k),
              "n_k",
              ADaPS_SCALAR_INT);
  ADaPS_store(cosmo,
              (void *)(lk_P),
              "lk_P",
              ADaPS_DEFAULT);
  ADaPS_store(cosmo,
              (void *)(lP_k),
              "lP_k_TF_all",
              ADaPS_DEFAULT);
  ADaPS_store(cosmo,
              (void *)(lP_k_gas),
              "lP_k_TF_gas",
              ADaPS_DEFAULT);
  ADaPS_store(cosmo,
              (void *)(lP_k_dark),
              "lP_k_TF_dark",
              ADaPS_DEFAULT);

  // Compute and store interpolation information for the unnormalized P(k)
  //   (needed by the normalization routine)
  init_interpolate(lk_P,lP_k,(size_t)n_k,gsl_interp_cspline,&interp);
  ADaPS_store_interp(cosmo,(void *)(interp),"lP_k_TF_all_interp");

  // Normalize P(k) to sigma_8
  norm=power_spectrum_normalization(*cosmo,PSPEC_LINEAR_TF,PSPEC_ALL_MATTER);
  for(int i=0;i<n_k;i++){
    lP_k[i]     +=norm;
    lP_k_gas[i] +=norm;
    lP_k_dark[i]+=norm;
  }

  // Create interpolation information for P(k) arrays
  init_interpolate(lk_P,lP_k,     (size_t)n_k,gsl_interp_cspline,&interp);
  init_interpolate(lk_P,lP_k_gas, (size_t)n_k,gsl_interp_cspline,&interp_gas);
  init_interpolate(lk_P,lP_k_dark,(size_t)n_k,gsl_interp_cspline,&interp_dark);

  // Store interpolation information for P(k) arrays
  ADaPS_store_interp(cosmo,
                     (void *)(interp),
                     "lP_k_TF_all_interp");
  ADaPS_store_interp(cosmo,
                    (void *)(interp_gas),
                    "lP_k_TF_gas_interp");
  ADaPS_store_interp(cosmo,
                     (void *)(interp_dark),
                     "lP_k_TF_dark_interp");

  // Compute end-point slopes for extrapolation
  double slope_lo_all =interpolate_derivative(interp,     lk_P[1]);
  double slope_hi_all =interpolate_derivative(interp,     lk_P[n_k-2]);
  double slope_lo_gas =interpolate_derivative(interp_gas, lk_P[1]);
  double slope_hi_gas =interpolate_derivative(interp_gas, lk_P[n_k-2]);
  double slope_lo_dark=interpolate_derivative(interp_dark,lk_P[1]);
  double slope_hi_dark=interpolate_derivative(interp_dark,lk_P[n_k-2]);
  ADaPS_store(cosmo,&slope_lo_all, "lP_k_TF_all_slope_lo", ADaPS_SCALAR_DOUBLE);
  ADaPS_store(cosmo,&slope_hi_all, "lP_k_TF_all_slope_hi", ADaPS_SCALAR_DOUBLE);
  ADaPS_store(cosmo,&slope_lo_gas, "lP_k_TF_gas_slope_lo", ADaPS_SCALAR_DOUBLE);
  ADaPS_store(cosmo,&slope_hi_gas, "lP_k_TF_gas_slope_hi", ADaPS_SCALAR_DOUBLE);
  ADaPS_store(cosmo,&slope_lo_dark,"lP_k_TF_dark_slope_lo",ADaPS_SCALAR_DOUBLE);
  ADaPS_store(cosmo,&slope_hi_dark,"lP_k_TF_dark_slope_hi",ADaPS_SCALAR_DOUBLE);

  SID_free(SID_FARG line);
  SID_log("Done.",SID_LOG_CLOSE);
}

double power_spectrum_normalization(cosmo_info *cosmo,
                                    int         mode,
                                    int         component){
  double     sigma2_unnorm;
  double     sigma_8;
  double     h_Hubble;
  double     norm;
  int        n_k;
  double    *lk_P;
  double     b_z;
  double     coefficient;
  int        n_int=10000;
  double     rel_accuracy=1e-5;
  double     limit_lo;
  double     limit_hi;
  double     r_val;
  double     abs_error;
  sigma2_integrand_params    params;
  gsl_integration_workspace *wspace;
  gsl_function               integrand;

  sigma_8 =((double *)ADaPS_fetch((ADaPS *)(cosmo),"sigma_8"))[0];
  h_Hubble=((double *)ADaPS_fetch((ADaPS *)(cosmo),"h_Hubble"))[0];
  n_k     =((int    *)ADaPS_fetch((cosmo),"n_k"))[0];
  lk_P    =(double  *)ADaPS_fetch((cosmo),"lk_P");

  /***************************************/
  /* Initialize data needed by integrand */
  /***************************************/
  b_z        =1.;
  coefficient=b_z*b_z/(TWO_PI*PI);

  // Initialize data needed by integrand
  limit_lo          =take_alog10(lk_P[0]);
  limit_hi          =take_alog10(lk_P[n_k-1]);
  params.component  =component;
  params.mode       =mode;
  params.z          =0.;
  params.cosmo      =cosmo;
  params.R          =8.*M_PER_MPC/h_Hubble;
  integrand.function=sigma2_integrand;
  integrand.params  =(void *)(&params);
  wspace            =gsl_integration_workspace_alloc(n_int);

  // Integrate spherical top hat
  gsl_integration_qag(&integrand,
            limit_lo,limit_hi,
            0,rel_accuracy,
            n_int,
            GSL_INTEG_GAUSS61,
            wspace,
            &sigma2_unnorm,&abs_error);
  sigma2_unnorm*=coefficient;

  // Normalize power spectrum to sigma_8
  norm=take_log10(sigma_8*sigma_8/sigma2_unnorm);

  // Clean-up
  gsl_integration_workspace_free(wspace);

  return(norm);
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double power_spectrum(double       k_interp,
                      double       redshift,
                      cosmo_info **cosmo,
                      int          mode,
                      int          component){

  double  rval;
  switch(mode){
    case PSPEC_LINEAR_TF:{       // Linear theory from transfer function
      char mode_name[ADaPS_NAME_LENGTH];
      char component_name[ADaPS_NAME_LENGTH];
      pspec_names(mode,component,mode_name,component_name);

      if(!ADaPS_exist(*cosmo,"lP_k_%s_%s",mode_name,component_name))
        init_power_spectrum_TF(cosmo);
      int          n_k   =((int        *)ADaPS_fetch(*cosmo,"n_k"))[0];
      double      *lk_P  =(double      *)ADaPS_fetch(*cosmo,"lk_P");
      double      *lP_k  =(double      *)ADaPS_fetch(*cosmo,"lP_k_%s_%s",       mode_name,component_name);
      interp_info *interp=(interp_info *)ADaPS_fetch(*cosmo,"lP_k_%s_%s_interp",mode_name,component_name);
      // Compute the needed normalization
      double  norm;
      if(redshift!=0.)
        norm=pow(linear_growth_factor(redshift,*cosmo),2.);
      else
        norm=1.;
      // Compute the power spectrum.  Exstrapolate from the
      //    given transfer function if needed.
      double lP_k_0;
      double lk_interp=take_log10(k_interp);
      if(k_interp<=0)
         return(0.);
      else if(lk_interp<lk_P[0]){
         double slope_lo=((double *)ADaPS_fetch(*cosmo,"lP_k_TF_%s_slope_lo",component_name))[0];
         double dl_k    =lk_interp-lk_P[0];
         lP_k_0         =lP_k[0]+dl_k*slope_lo;
      }
      else if(lk_interp>lk_P[n_k-1]){
         double slope_hi=((double *)ADaPS_fetch(*cosmo,"lP_k_TF_%s_slope_hi",component_name))[0];
         double dl_k    =lk_interp-lk_P[n_k-1];
         lP_k_0         =lP_k[n_k-1]+dl_k*slope_hi;
      }
      else
         lP_k_0=interpolate(interp,lk_interp);
      rval=norm*take_alog10(lP_k_0);
      break;
      }
    default:
      SID_trap_error("Unsupported mode (%d) passed to power_spectrum().",ERROR_LOGIC,mode);
      break;
  }
  return(rval);
}

