#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>

/*****************************************/
/* R, k and M as functions of each other */
/*****************************************/
double R_of_k(double k){
  return(TWO_PI/k);
}
double k_of_R(double R){
  return(TWO_PI/R);
}
double M_of_R(double R,double z,cosmo_info *cosmo){
  double Omega_M;
  double M;
  Omega_M=((double *)ADaPS_fetch((ADaPS *)(cosmo),"Omega_M"))[0];
  M      =Omega_M*FOUR_THIRDS_PI*pow(R,3.)*rho_crit_z(z,cosmo);
  return(M);
}
double M_of_k(double k,double z,cosmo_info *cosmo){
  double R;
  double Omega_M;
  double M;
  R      =R_of_k(k);
  Omega_M=((double *)ADaPS_fetch((ADaPS *)(cosmo),"Omega_M"))[0];
  M      =Omega_M*FOUR_THIRDS_PI*pow(R,3.)*rho_crit_z(z,cosmo);
  return(M);
}
double R_of_M(double M,double z,cosmo_info *cosmo){
  double Omega_M;
  double R3,R;
  Omega_M=((double *)ADaPS_fetch((ADaPS *)(cosmo),"Omega_M"))[0];
  R3     =M/(Omega_M*FOUR_THIRDS_PI*rho_crit_z(z,cosmo));
  R      =pow(R3,ONE_THIRD);
  return(R);
}
double k_of_M(double M,double z,cosmo_info *cosmo){
  double R;
  double k;
  R=R_of_M(M,z,cosmo);
  k=k_of_R(R);
  return(k);
}

void pspec_names(int   mode,
		 int   component,
		 char *mode_name,
		 char *component_name){
  switch(mode){
  case PSPEC_LINEAR_TF:
    sprintf(mode_name,"TF");
    break;
  case PSPEC_LINEAR_BBKS:
    sprintf(mode_name,"BBKS");
    break;
  case PSPEC_NONLINEAR_JAIN:
    sprintf(mode_name,"JAIN");
    break;
  case PSPEC_NONLINEAR_PD:
    sprintf(mode_name,"NL_PD");
    break;
  case PSPEC_NONLINEAR_SMITH:
    sprintf(mode_name,"NL_Smith");
    break;
  default:
    fprintf(stderr,"Invalid PSPEC mode in pspec_names!\n");
    break;
  }

  switch(component){
  case PSPEC_ALL_MATTER:
    sprintf(component_name,"all");
    break;
  case PSPEC_DARK_MATTER:
    sprintf(component_name,"dark");
    break;
  case PSPEC_BARYON:
    sprintf(component_name,"gas");
    break;
  default:
    fprintf(stderr,"Invalid PSPEC component in pspec_names!\n");
    break;
  }
}

/****************************************/
/* Initialize the power spectrum at z=0 */ 
/****************************************/
void init_power_spectrum_TF(cosmo_info **cosmo){
  FILE   *fp;
  char   *line=NULL;
  char    filename_TF[MAX_FILENAME_LENGTH];
  char    filename_nl[MAX_FILENAME_LENGTH];
  size_t  line_length=0;
  int     i;
  int     n_k;
  size_t  n_k_dim;
  int     n_k_tmp;
  double  k_P;
  double *lk_P;
  double *lk_P_tmp;
  double *lP_k;
  double *lP_k_dark;
  double *lP_k_gas;
  double *lP_k_nl;
  double *lP_k_tmp;
  double *lP_k_gas_tmp;
  double *lP_k_dark_tmp;
  double *lP_k_nl_tmp;
  interp_info *interp;
  interp_info *interp_dark;
  interp_info *interp_gas;
  interp_info *interp_nl;
  double  M_8;
  double  sigma2_unnorm;
  double  Omega_b;
  double  Omega_M;
  double  h_Hubble;
  double  norm;
  double  norm_nl;
  double  n_spectral;
  double  M_WDM,R_WDM;

  // Names of the files where the transfer function and non-linear power spectrum are stored
  sprintf(filename_TF,"%s/transfer_function.dat",       GBP_DATA_DIR);
  sprintf(filename_nl,"%s/nonlinear_power_spectrum.dat",GBP_DATA_DIR);
  
  n_spectral =((double *)ADaPS_fetch((*cosmo),"n_spectral"))[0];
  h_Hubble   =((double *)ADaPS_fetch((*cosmo),"h_Hubble"))[0];
  Omega_M    =((double *)ADaPS_fetch((*cosmo),"Omega_M"))[0];
  Omega_b    =((double *)ADaPS_fetch((*cosmo),"Omega_b"))[0];

  // Read the transfer function
  fp       =fopen(filename_TF,"r");
  n_k      =count_lines_data(fp);
  n_k_dim  =(size_t)n_k;
  lk_P     =(double *)SID_malloc(sizeof(double)*n_k);
  lP_k     =(double *)SID_malloc(sizeof(double)*n_k);
  lP_k_gas =(double *)SID_malloc(sizeof(double)*n_k);
  lP_k_dark=(double *)SID_malloc(sizeof(double)*n_k);
  lP_k_nl  =(double *)SID_malloc(sizeof(double)*n_k);
  for(i=0;i<n_k;i++){
    grab_next_line_data(fp,&line,&line_length);
    grab_double(line,1,&(lk_P[i]));
    grab_double(line,2,&(lP_k_dark[i]));
    grab_double(line,3,&(lP_k_gas[i]));
    lP_k[i]=((Omega_M-Omega_b)*lP_k_dark[i]+Omega_b*lP_k_gas[i])/Omega_M;
  }
  fclose(fp);

  // Read tabulated non-linear power spectrum
  fp      =fopen(filename_nl,"r");
  n_k_tmp =count_lines_data(fp);
  lk_P_tmp=(double *)SID_malloc(sizeof(double)*n_k_tmp);
  lP_k_tmp=(double *)SID_malloc(sizeof(double)*n_k_tmp);
  for(i=0;i<n_k_tmp;i++){
    grab_next_line_data(fp,&line,&line_length);
    grab_double(line,1,&(lk_P_tmp[i]));
    grab_double(line,2,&(lP_k_tmp[i]));
  }
  fclose(fp);
  init_interpolate(lk_P_tmp,lP_k_tmp,(size_t)n_k_tmp,gsl_interp_cspline,&interp);
  for(i=0;i<n_k;i++)
    lP_k_nl[i]=interpolate(interp,lk_P[i]);
  SID_free(SID_FARG lk_P_tmp);
  SID_free(SID_FARG lP_k_tmp);
  free_interpolate(&interp);

  // Take the log of lk_P
  for(i=0;i<n_k;i++)
    lk_P[i]=take_log10(lk_P[i]/(M_PER_MPC/h_Hubble));
  
  // Supress small-scale power following the
  //    ENS recipe if M_ENS or M_WDM is set > 0
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
         lP_k_nl[i]  =0.;
      }
      else{
         lP_k[i]     *=exp(-0.5*k_P*R_WDM-0.5*pow(k_P*R_WDM,2.));
         lP_k_gas[i] *=exp(-0.5*k_P*R_WDM-0.5*pow(k_P*R_WDM,2.));
         lP_k_dark[i]*=exp(-0.5*k_P*R_WDM-0.5*pow(k_P*R_WDM,2.));
         lP_k_nl[i]  *=exp( 2.0*(-0.5*k_P*R_WDM-0.5*pow(k_P*R_WDM,2.)));
      }
    }
  }
    
  // Compute power spectrum
  for(i=0;i<n_k;i++){
    lP_k[i]     =2.*take_log10(lP_k[i])     +(n_spectral)*lk_P[i];
    lP_k_gas[i] =2.*take_log10(lP_k_gas[i]) +(n_spectral)*lk_P[i];
    lP_k_dark[i]=2.*take_log10(lP_k_dark[i])+(n_spectral)*lk_P[i];
    lP_k_nl[i]  =   take_log10(lP_k_nl[i]);
  }

  // Normalize at large scales
  norm   =lP_k[0];
  norm_nl=lP_k_nl[0];
  for(i=0;i<n_k;i++){
    lP_k[i]     -=norm;
    lP_k_gas[i] -=norm;
    lP_k_dark[i]-=norm;
    lP_k_nl[i]  -=norm_nl;
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
  ADaPS_store(cosmo,
              (void *)(lP_k_nl),
              "lP_k_NL_Smith_all",
              ADaPS_DEFAULT);

  // Compute and store interpolation information for the unnormalized P(k)
  //   (needed by the normalization routine)
  init_interpolate(lk_P,lP_k,   (size_t)n_k,gsl_interp_cspline,&interp);
  init_interpolate(lk_P,lP_k_nl,(size_t)n_k,gsl_interp_cspline,&interp_nl);
  ADaPS_store_interp(cosmo,
                     (void *)(interp),
                     "lP_k_TF_all_interp");
  ADaPS_store_interp(cosmo,
                     (void *)(interp_nl),
                     "lP_k_NL_Smith_all_interp");

  // Normalize P(k) to sigma_8
  norm=power_spectrum_normalization(*cosmo,PSPEC_LINEAR_TF,PSPEC_ALL_MATTER);
  norm_nl=lP_k[0]-lP_k_nl[0]+norm;
  for(i=0;i<n_k;i++){
    lP_k[i]     +=norm;
    lP_k_gas[i] +=norm;
    lP_k_dark[i]+=norm;
    lP_k_nl[i]  +=norm_nl;
  }

  // Create interpolation information for P(k) arrays
  ADaPS_remove(cosmo,"lP_k_TF_all_interp");
  ADaPS_remove(cosmo,"lP_k_NL_Smith_all_interp");
  init_interpolate(lk_P,lP_k,     (size_t)n_k,gsl_interp_cspline,&interp);
  init_interpolate(lk_P,lP_k_nl,  (size_t)n_k,gsl_interp_cspline,&interp_nl);
  init_interpolate(lk_P,lP_k_gas, (size_t)n_k,gsl_interp_cspline,&interp_gas);
  init_interpolate(lk_P,lP_k_dark,(size_t)n_k,gsl_interp_cspline,&interp_dark);

  // Store interpolation information for P(k) arrays
  ADaPS_store_interp(cosmo,
                     (void *)(interp),
                     "lP_k_TF_all_interp");
  ADaPS_store_interp(cosmo,
                     (void *)(interp_nl),
                     "lP_k_NL_Smith_all_interp");
  ADaPS_store_interp(cosmo,
                    (void *)(interp_gas),
                    "lP_k_TF_gas_interp");
  ADaPS_store_interp(cosmo,
                     (void *)(interp_dark),
                     "lP_k_TF_dark_interp");
  ADaPS_store_interp(cosmo,
                    (void *)(interp_nl),
                    "lP_k_NL_Smith_all_interp");

  SID_free(SID_FARG line);
}

double power_spectrum(double      k_interp,
		      double      redshift,
		      cosmo_info *cosmo,
		      int         mode,
		      int         component){
  int     n_k;
  double *lk_P;
  double *lP_k;
  interp_info *interp;
  double  norm;
  double  rval;
  char    mode_name[ADaPS_NAME_LENGTH];
  char    component_name[ADaPS_NAME_LENGTH];
  char    lP_name[ADaPS_NAME_LENGTH];
  char    d2lP_name[ADaPS_NAME_LENGTH];

  switch(mode){
  case PSPEC_LINEAR_TF:       // Linear theory from transfer function
  case PSPEC_NONLINEAR_JAIN:  // Jain et al correlation function needs linear power spectrum
  case PSPEC_NONLINEAR_SMITH: // Smith et al '02; generated by CAMB
    pspec_names(mode,
		component,
		mode_name,
		component_name);
    sprintf(lP_name,  "lP_k_%s_%s",       mode_name,component_name);
    sprintf(d2lP_name,"lP_k_%s_%s_interp",mode_name,component_name);
    if(!ADaPS_exist(cosmo,lP_name))
      init_power_spectrum_TF(&cosmo);
    n_k   =((int   *)ADaPS_fetch(cosmo,"n_k"))[0];
    lk_P  =(double *)ADaPS_fetch(cosmo,"lk_P");
    lP_k  =(double *)ADaPS_fetch(cosmo,lP_name);
    interp=(interp_info *)ADaPS_fetch(cosmo,d2lP_name);
    if(redshift!=0.)
      norm=pow(linear_growth_factor(redshift,cosmo),2.);
    else
      norm=1.;
    rval=norm*take_alog10(interpolate(interp,take_log10(k_interp)));
    break;
/*
  case PSPEC_LINEAR_BBKS:     // Compute BBKS with Smith et al code
    rval=P_NL(k_interp,redshift,cosmo,0);
    break;  
  case PSPEC_NONLINEAR_PD:    // Peacock and Dodds '96
    rval=P_NL(k_interp,redshift,cosmo,1);
    break;
    rval=P_NL(k_interp,redshift,cosmo,2);
    break;
*/
  }
  return(rval);
}

/******************************/
/* Linear growth factor stuff */
/******************************/
double dlogDplus_dloga(double      a,
		       cosmo_info *cosmo){
  double Ez;
  double h_Hubble;
  double Omega_M;
  double Omega_k;
  h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
  Omega_M =((double *)ADaPS_fetch(cosmo,"Omega_M"))[0];
  Omega_k =((double *)ADaPS_fetch(cosmo,"Omega_k"))[0];
  Ez      =H_z(z_of_a(a),cosmo)/(100.*h_Hubble);
  return((2.5/Dplus(a,cosmo)-1.5*Omega_M/a-Omega_k)/pow(a*Ez,2.));
}
double dDplus_da(double      a,
		 cosmo_info *cosmo){
  double      h_Hubble;
  double      Ez;
  h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
  Ez      =H_z(z_of_a(a),cosmo)/(1e2*h_Hubble);
  return(1.0/pow(a*Ez,3.0));
}
double Dplus_integrand(double  a,
		       void   *cosmo){
  return(dDplus_da(a,(cosmo_info *)cosmo));
}
double Dplus(double a,cosmo_info *cosmo){
  int    n_int=1000;
  double h_Hubble;
  double Omega_M;
  double Ez;
  double rel_accuracy=1e-5;
  double abs_error;
  double limit_lo;
  double limit_hi;
  double integral;
  double r_val;
  gsl_integration_workspace *wspace;
  gsl_function               integrand;

  // Initialize integral
  limit_lo          =0.;
  limit_hi          =a;
  integrand.function=Dplus_integrand;
  integrand.params  =(void *)cosmo;
  wspace            =gsl_integration_workspace_alloc(n_int);

  gsl_integration_qags(&integrand,
		       limit_lo,limit_hi,
		       0,rel_accuracy,
		       n_int,
		       wspace,
		       &r_val,&abs_error); // use qags for singularity at a=0

  h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
  Omega_M =((double *)ADaPS_fetch(cosmo,"Omega_M"))[0];
  Ez      =H_z(z_of_a(a),cosmo)/(100.*h_Hubble);

  r_val  *=2.5*Omega_M*Ez;

  // Clean-up
  gsl_integration_workspace_free(wspace);

  return(r_val);
}
double linear_growth_factor(double       redshift,
			    cosmo_info  *cosmo){
  double        Dplus_a;
  double        Dplus_1;
  static double b_z;
  static double z_last=-42.;
  if(redshift!=z_last){
    Dplus_a=Dplus(a_of_z(redshift),cosmo);
    Dplus_1=Dplus(1.,              cosmo);
    b_z    =Dplus_a/Dplus_1;
    z_last =redshift;
  }
  return(b_z);
}

/*****************/
/* sigma^2 stuff */
/*****************/
double W_k_tophat(double kR){
  return(3.*(sin(kR)-kR*cos(kR))/pow(kR,3.));
}
typedef struct sigma2_integrand_params_struct sigma2_integrand_params;
struct sigma2_integrand_params_struct {
  double      R; 
  cosmo_info *cosmo; 
  double      z; 
  int         mode; 
  int         component;
};
double sigma2_integrand(double  k,
			void   *params_in){
  double P_k_tmp;
  double W_k_tmp;
  double r_val;
  sigma2_integrand_params *params;
  params =(sigma2_integrand_params *)params_in;
  P_k_tmp=power_spectrum(k,params->z,params->cosmo,params->mode,params->component);
  W_k_tmp=W_k_tophat(k*params->R);
  r_val  =k*k*P_k_tmp*W_k_tmp*W_k_tmp;
  return(r_val);
}
void init_power_spectrum_variance(cosmo_info **cosmo,double z,int mode,int component){
  int        i;
  int        n_k;
  size_t     n_k_dim;
  double    *lk_P;
  double    *sigma2;
  interp_info *interp;
  double     b_z;
  double     coefficient;
  char       mode_name[ADaPS_NAME_LENGTH];
  char       component_name[ADaPS_NAME_LENGTH];
  char       sigma2_name[ADaPS_NAME_LENGTH];
  char       d2sigma2_name[ADaPS_NAME_LENGTH];
  double     lk_step;
  double     integral;
  int        n_int=5000;
  double     rel_accuracy=5e-5;
  double     limit_lo;
  double     limit_hi;
  double     r_val;
  double     abs_error;
  sigma2_integrand_params    params;
  gsl_integration_workspace *wspace;
  gsl_function               integrand;

  z=0.;

  // Compute sigma^2 coefficient
  b_z        =1.;
  coefficient=b_z*b_z/(TWO_PI*PI);

  // Make sure k's are initialized
  if(!ADaPS_exist((*cosmo),"n_k")){
    if(mode==PSPEC_LINEAR_TF)
      init_power_spectrum_TF(cosmo);
    else{
      n_k        =200;
      n_k_dim    =(size_t)n_k;
      lk_P       =(double *)SID_malloc(sizeof(double)*n_k_dim);
      lk_P[0]    =take_log10(k_of_R(1e0*M_PER_KPC));
      lk_P[n_k-1]=take_log10(k_of_R(2e4*M_PER_MPC));
      lk_step    =(lk_P[n_k-1]-lk_P[0])/(double)(n_k-1);
      for(i=1;i<n_k-1;i++)
         lk_P[i]=lk_P[i-1]+lk_step;
      ADaPS_store(cosmo,
                  (void *)(&n_k),
                  "n_k",
                  ADaPS_SCALAR_INT);
      ADaPS_store(cosmo,
                  (void *)(lk_P),
                  "lk_P",
                  ADaPS_DEFAULT);
    }
  }
  n_k    =((int   *)ADaPS_fetch((*cosmo),"n_k"))[0];
  lk_P   =(double *)ADaPS_fetch((*cosmo),"lk_P");
  n_k_dim=(size_t)n_k;

  // Initialize integral
  limit_lo          =take_alog10(lk_P[0]);
  limit_hi          =take_alog10(lk_P[n_k-1]);
  integrand.function=sigma2_integrand;
  params.component  =component;
  params.mode       =mode;
  if(mode==PSPEC_LINEAR_TF || mode==PSPEC_NONLINEAR_SMITH)
    params.z        =0.;
  else
    params.z        =z;
  params.cosmo      =(*cosmo);
  integrand.params  =(void *)(&params);
  wspace            =gsl_integration_workspace_alloc(n_int);

  // Loop over the whole range of scales
  sigma2=(double *)SID_malloc(sizeof(double)*n_k_dim);
  for(i=0;i<n_k;i++){

    // Perform spherical top-hat integral
    params.R=R_of_k(take_alog10(lk_P[i]));
    gsl_integration_qag(&integrand,
			limit_lo,limit_hi,
			0,rel_accuracy,
			n_int,
			GSL_INTEG_GAUSS61,
			wspace,
			&integral,&abs_error);
    sigma2[i]=coefficient*integral;
  }
  init_interpolate(lk_P,sigma2,(size_t)n_k,gsl_interp_cspline,&interp);

  pspec_names(mode,component,mode_name,component_name);
  sprintf(sigma2_name,  "sigma2_k_%s_%s",       mode_name,component_name);
  sprintf(d2sigma2_name,"sigma2_k_%s_%s_interp",mode_name,component_name);
  ADaPS_store(cosmo,
              (void *)(sigma2),
              "sigma2_k_%s_%s",
              ADaPS_DEFAULT,
              mode_name,component_name);
  ADaPS_store_interp(cosmo,
                     (void *)(interp),
                     "sigma2_k_%s_%s_interp",
                     mode_name,component_name);

  // Clean-up
  gsl_integration_workspace_free(wspace);
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
  integrand.function=sigma2_integrand;
  params.component  =component;
  params.mode       =mode;
  params.z          =0.;
  params.cosmo      =cosmo;
  params.R          =8.*M_PER_MPC/h_Hubble;
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

  return(norm);
}

double power_spectrum_variance(double      k_interp,
			       double      redshift,
			       cosmo_info *cosmo,
			       int         mode,
			       int         component){
  int     n_k;
  double *lk_P;
  double *sigma2;
  interp_info *interp;
  double  rval=0.;
  char    mode_name[ADaPS_NAME_LENGTH];
  char    component_name[ADaPS_NAME_LENGTH];
  char    sigma2_name[ADaPS_NAME_LENGTH];
  char    d2sigma2_name[ADaPS_NAME_LENGTH];
  static double z_last=-42.;
  static double norm;
  
  pspec_names(mode,component,mode_name,component_name);
  sprintf(sigma2_name,  "sigma2_k_%s_%s",       mode_name,component_name);
  sprintf(d2sigma2_name,"sigma2_k_%s_%s_interp",mode_name,component_name);

  if(!ADaPS_exist(cosmo,sigma2_name))
    init_power_spectrum_variance(&cosmo,redshift,mode,component);

  n_k     =((int   *)ADaPS_fetch(cosmo,"n_k"))[0];
  lk_P    =(double *)ADaPS_fetch(cosmo,"lk_P");
  sigma2  =(double *)ADaPS_fetch(cosmo,sigma2_name);
  interp  =(interp_info *)ADaPS_fetch(cosmo,d2sigma2_name);

  norm=pow(linear_growth_factor(redshift,cosmo),2.);

  rval=norm*interpolate(interp,take_log10(k_interp));

  return(rval);
}

double dln_Inv_sigma_dlogM(cosmo_info *cosmo,
			   double      M_interp,
			   double      z,
			   int         mode,
			   int         component){
  double       r_val;
  double      *ln_Inv_sigma;
  double       derr;
  int          i;
  int          n_k;
  size_t  n_k_dim;
  double      *lk_P;
  double      *lM_k;
  interp_info *interp_ln_Inv_sigma;
  char         mode_name[ADaPS_NAME_LENGTH];
  char         component_name[ADaPS_NAME_LENGTH];
  char         d2ln_Inv_sigma_name[ADaPS_NAME_LENGTH];

  // Fetch/initialize ln(sigma^-1) and its derivative wrt log M
  pspec_names(mode,component,mode_name,component_name);
  sprintf(d2ln_Inv_sigma_name,"ln_Inv_sigma_M_%s_%s_interp",mode_name,component_name);
  if(!ADaPS_exist(cosmo,d2ln_Inv_sigma_name)){
    n_k         =((int   *)ADaPS_fetch(cosmo,"n_k"))[0];
    lk_P        =(double *)ADaPS_fetch(cosmo,"lk_P");
    // Fetch/initialize M(k)
    if(!ADaPS_exist(cosmo,"lM_k")){
      lM_k   =(double *)SID_malloc(sizeof(double)*n_k);
      n_k_dim=(size_t)n_k;
      for(i=0;i<n_k;i++)
         lM_k[i]=take_log10(M_of_k(take_alog10(lk_P[i]),z,cosmo));
      ADaPS_store((&cosmo),
		 "lM_k",
		 (void *)(lM_k),
		 ADaPS_DOUBLE,
		 1,&n_k_dim);
      SID_free(SID_FARG lM_k);
    }
    lM_k        =(double  *)ADaPS_fetch(cosmo,"lM_k");
    ln_Inv_sigma=(double *)SID_malloc(sizeof(double)*n_k);
    n_k_dim     =(size_t)n_k;
    for(i=0;i<n_k;i++)
      ln_Inv_sigma[i]=take_ln(1./sqrt(power_spectrum_variance(take_alog10(lk_P[i]),
                                                              0.,
                                                              cosmo,
                                                              mode,
                                                              component)));
    init_interpolate(lM_k,
                     ln_Inv_sigma,
                     (size_t)n_k,
                     gsl_interp_cspline,
                     &interp_ln_Inv_sigma);
    ADaPS_store_interp((&cosmo),
                       (void *)(interp_ln_Inv_sigma),
                       d2ln_Inv_sigma_name);
    SID_free(SID_FARG ln_Inv_sigma);
  }
  interp_ln_Inv_sigma=(interp_info *)ADaPS_fetch(cosmo,d2ln_Inv_sigma_name);
  r_val=
    interpolate_derivative(interp_ln_Inv_sigma,
		           take_log10(M_interp));
  return(r_val);
}

double dln_sigma_dlnM(cosmo_info *cosmo,
		      double      M_interp,
		      double      z,
		      int         mode,
		      int         component){
  double       r_val;
  double      *ln_sigma;
  double       derr;
  int          i;
  int          n_k;
  size_t       n_k_dim;
  double      *lk_P;
  double      *lM_k;
  interp_info *interp_ln_sigma;
  char         mode_name[ADaPS_NAME_LENGTH];
  char         component_name[ADaPS_NAME_LENGTH];
  char         d2ln_sigma_name[ADaPS_NAME_LENGTH];

  // Fetch/initialize ln(sigma^-1) and its derivative wrt log M
  pspec_names(mode,component,mode_name,component_name);
  sprintf(d2ln_sigma_name,"ln_sigma_lnM_%s_%s_interp",mode_name,component_name);
  if(!ADaPS_exist(cosmo,d2ln_sigma_name)){
    n_k    =((int   *)ADaPS_fetch(cosmo,"n_k"))[0];
    lk_P   =(double *)ADaPS_fetch(cosmo,"lk_P");
    lM_k   =(double *)SID_malloc(sizeof(double)*n_k);
    n_k_dim=(size_t)n_k;
    for(i=0;i<n_k;i++)
      lM_k[i]=take_ln(M_of_k(take_alog10(lk_P[i]),z,cosmo));

    ln_sigma=(double *)SID_malloc(sizeof(double)*n_k);
    n_k_dim =(size_t)n_k;
    for(i=0;i<n_k;i++)
      ln_sigma[i]=take_ln(sqrt(power_spectrum_variance(take_alog10(lk_P[i]),
						       0.,
						       cosmo,
						       mode,
						       component)));
    init_interpolate(lM_k,
                     ln_sigma,
                     n_k_dim,
                     gsl_interp_cspline,
                     &interp_ln_sigma);
    ADaPS_store_interp((&cosmo),
                       (void *)(interp_ln_sigma),
                       d2ln_sigma_name);
    free_interpolate(&interp_ln_sigma);
    SID_free(SID_FARG lM_k);
    SID_free(SID_FARG ln_sigma);
  }
  interp_ln_sigma=(interp_info *)ADaPS_fetch(cosmo,d2ln_sigma_name);
  r_val=
    interpolate_derivative(interp_ln_sigma,
		           take_ln(M_interp));
  return(r_val);
}

/***************************************/
/* Stuff related to spherical collapse */
/***************************************/
double M_sc(double      z,
	    cosmo_info *cosmo,
	    int         mode,
	    int         component){
  int     i;
  int     n_k;
  size_t  n_k_dim;
  double *lM_k;
  double *lk_P;
  double *sigma2;
  interp_info *interp;
  double  delta_sc=1.686;
  double  b_z;
  double  r_val;
  char    mode_name[ADaPS_NAME_LENGTH];
  char    component_name[ADaPS_NAME_LENGTH];
  char    sigma2_name[ADaPS_NAME_LENGTH];
  char    d2sigma2_name[ADaPS_NAME_LENGTH];
  static double M_sc_last;
  static double z_last=-42.;

  if(z!=z_last){
    
    // Set/initialize variance
    pspec_names(mode,component,mode_name,component_name);
    sprintf(sigma2_name,  "sigma2_k_%s_%s",       mode_name,component_name);
    sprintf(d2sigma2_name,"sigma2_k_%s_%s_interp",mode_name,component_name);
    if(!ADaPS_exist(cosmo,sigma2_name))
      init_power_spectrum_variance(&cosmo,z,mode,component);
    sigma2=(double      *)ADaPS_fetch(cosmo,sigma2_name);
    interp=(interp_info *)ADaPS_fetch(cosmo,d2sigma2_name);
    
    b_z=linear_growth_factor(z,cosmo);

    /*
    r_val=bisect_array(interp,
		       delta_sc*delta_sc/(b_z*b_z),
		       1e-4);
    */
    SID_trap_error("Need to implement bisect_array here.",ERROR_LOGIC);
    M_sc_last=M_of_k(take_alog10(r_val),z,cosmo);
    z_last   =z;
  }

  return(M_sc_last);
}
double lk_sc(double      z,
	     cosmo_info *cosmo,
	     int         mode,
	     int         component){
  return(take_log10(k_of_M(M_sc(z,cosmo,mode,component),z,cosmo)));
}
