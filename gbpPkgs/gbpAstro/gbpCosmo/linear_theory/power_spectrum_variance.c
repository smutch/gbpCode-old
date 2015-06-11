#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

void init_power_spectrum_variance(cosmo_info **cosmo,int mode,int component){

  SID_log("Initializing P(k) variance...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Make sure the power spectrum is initialized
  if(!ADaPS_exist((*cosmo),"n_k")){
    if(mode==PSPEC_LINEAR_TF)
      init_power_spectrum_TF(cosmo);
    else
      SID_trap_error("Given mode (%d) not supported in init_power_spectrum_variance().",ERROR_LOGIC,mode);
  }

  // Fetch the k array and its size
  int     n_k =((int    *)ADaPS_fetch((*cosmo),"n_k"))[0];
  double *lk_P= (double *)ADaPS_fetch((*cosmo),"lk_P");

  // Initialize lower integration limit
  double limit_lo_init=0.;
  if(ADaPS_exist((*cosmo),"box_size")){
     double h_Hubble=((double *)ADaPS_fetch((*cosmo),"h_Hubble"))[0];
     double box_size=((double *)ADaPS_fetch((*cosmo),"box_size"))[0];
     SID_log("Using integration limits for a box_size of %.2lf h^-1 [Mpc].",SID_LOG_COMMENT,box_size*h_Hubble/M_PER_MPC);
     limit_lo_init=TWO_PI/box_size;
  }

  // Compute sigma^2 coefficient
  double b_z        =1.;
  double coefficient=b_z*b_z/(TWO_PI*PI);

  // Initialize integral
  int    n_int       =1024;
  double abs_accuracy=0.;
  double rel_accuracy=1e-4;
  sigma2_integrand_params    params;
  gsl_integration_workspace *wspace;
  gsl_function               integrand;
  params.component  =component;
  params.mode       =mode;
  params.z          =0.;
  params.cosmo      =(*cosmo);
  integrand.function=sigma2_integrand;
  integrand.params  =(void *)(&params);
  wspace            =gsl_integration_workspace_alloc(128*n_int);

  // Loop over the whole range of scales
  double *sigma2    =(double *)SID_malloc(sizeof(double)*n_k);
  double  sigma2_min=1e10;
  for(int i=0;i<n_k;i++){
    // Set top-hat size 
    params.R=R_of_k(take_alog10(lk_P[i]));

    // This is an oscillating function.  Integrate
    //    in chunks integral periods long until 
    //    convergence is found.
    double rel_diff   =10.*rel_accuracy;
    double integral   =0.;
    int    i_order    =2;
    int    n_order_max=10000;
    double k_period   =TWO_PI/params.R;
    double limit_lo   =limit_lo_init;
    double limit_hi   =(double)i_order*k_period;
    if(limit_lo<limit_hi){
       while(rel_diff>=0.25*rel_accuracy){
   //fprintf(stderr,"A:%le %le\n",limit_lo,limit_hi);
   //for(double k_i=limit_lo;k_i<limit_hi;k_i+=(limit_hi-limit_lo)*0.01) printf("%le %le\n",k_i,sigma2_integrand(k_i,(void *)(&params)));
   //fprintf(stderr,"B:%le %le\n",limit_lo,limit_hi);
   //SID_exit(0);
   
          // Perform spherical top-hat integral
          //    for this iteration
          double dintegral;
          double abs_error;
          gsl_integration_qag(&integrand,
               limit_lo,limit_hi,
               abs_accuracy,rel_accuracy,
               n_int,
               GSL_INTEG_GAUSS15,
               wspace,
               &dintegral,&abs_error);
          //gsl_integration_qags(&integrand,
          //                     limit_lo,limit_hi,
          //                     abs_accuracy,rel_accuracy,
          //                     n_int,
          //                     wspace,
          //                     &dintegral,&abs_error); 
   
          // Update integral and compute relative change
          //    from this iteration
          if(integral>0.)
             rel_diff=dintegral/integral;
          integral+=dintegral;
          // Add an order to the integral range
          //    and adjust limits accordingly
          i_order++;
          if(i_order>=n_order_max)
             SID_trap_error("Variance integral is not converging for scale R=%le [Mpc]",ERROR_LOGIC,params.R/M_PER_MPC);
          limit_lo=limit_hi;
          limit_hi=(double)i_order*k_period;
       }
    }

    // Apply coefficient
    sigma2[i]=coefficient*integral;

    // Find the smallest non-zero value
    if(sigma2[i]>0.) sigma2_min=MIN(sigma2_min,sigma2[i]);
  }

  // Because we may need 1/sigma in various places, we
  //    can not have sigma=0 anywhere ... so, set all zeros
  //    to the the minimum value.  Zeros can come
  //    about principly due to finite box size effects.
  for(int i=0;i<n_k;i++) sigma2[i]=MAX(sigma2[i],sigma2_min);

  // Initialize interpolation
  interp_info *interp;
  init_interpolate(lk_P,sigma2,(size_t)n_k,gsl_interp_cspline,&interp);

  // Store sigma^2 array and its interpolation
  char mode_name[ADaPS_NAME_LENGTH];
  char component_name[ADaPS_NAME_LENGTH];
  pspec_names(mode,component,mode_name,component_name);
  ADaPS_store(cosmo,
              (void *)(sigma2),
              "sigma2_k_%s_%s",
              ADaPS_DEFAULT,
              mode_name,component_name);
  ADaPS_store_interp(cosmo,
                     (void *)(interp),
                     "sigma2_k_%s_%s_interp",
                     mode_name,component_name);

//double *y_temp=(double *)SID_malloc(sizeof(double)*n_k);
//for(int i_k=0;i_k<n_k;i_k++) y_temp[i_k]=1./sqrt(sigma2[i_k]);
//interp_info *interp2;
//init_interpolate(lk_P,y_temp,(size_t)n_k,gsl_interp_cspline,&interp2);
//
//for(double lk=lk_P[0];lk<lk_P[n_k-1];lk+=0.05) printf("%le %le %le %le %le\n",lk,interpolate(interp,lk),1./sqrt(interpolate(interp,lk)),interpolate_derivative(interp2,lk),interpolate_derivative(interp,lk));
//SID_exit(0);

  // Clean-up
  gsl_integration_workspace_free(wspace);

  SID_log("Done.",SID_LOG_CLOSE);
}

double power_spectrum_variance(double       k_interp,
                               double       redshift,
                               cosmo_info **cosmo,
                               int          mode,
                               int          component){
  // Set array names
  char mode_name[ADaPS_NAME_LENGTH];
  char component_name[ADaPS_NAME_LENGTH];
  char sigma2_name[ADaPS_NAME_LENGTH];
  char d2sigma2_name[ADaPS_NAME_LENGTH];
  pspec_names(mode,component,mode_name,component_name);
  sprintf(sigma2_name,  "sigma2_k_%s_%s",       mode_name,component_name);
  sprintf(d2sigma2_name,"sigma2_k_%s_%s_interp",mode_name,component_name);

  // Initialize arrays if needed
  if(!ADaPS_exist(*cosmo,sigma2_name))
    init_power_spectrum_variance(cosmo,mode,component);

  // Fetch arrays
  int          n_k   =((int   *)ADaPS_fetch(*cosmo,"n_k"))[0];
  double      *lk_P  =(double *)ADaPS_fetch(*cosmo,"lk_P");
  double      *sigma2=(double *)ADaPS_fetch(*cosmo,sigma2_name);
  interp_info *interp=(interp_info *)ADaPS_fetch(*cosmo,d2sigma2_name);

  // Set result
  double norm=pow(linear_growth_factor(redshift,*cosmo),2.);
  double rval=norm*interpolate(interp,take_log10(k_interp));

  return(rval);
}

