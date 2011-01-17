#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>

double a_of_z(double z){
  return(1./(1.+z));
}

double z_of_a(double a){
  return((1./a)-1.);
}

double R_Delta_z(double      M_Delta,
		 double      Delta,
		 double      redshift,
		 cosmo_info *cosmo){
  double rho_crit;
  rho_crit=rho_crit_z(redshift,cosmo);
  return(pow(M_Delta/(FOUR_THIRDS_PI*Delta*rho_crit),ONE_THIRD));
}

double Omega_z(double      redshift,
               cosmo_info *cosmo){
  double Omega_M,Omega_k,Omega_Lambda;
  double Ez,one_plus_z_cube;
  Omega_M        =((double *)ADaPS_fetch(cosmo,"Omega_M"))[0];
  Omega_k        =((double *)ADaPS_fetch(cosmo,"Omega_k"))[0];
  Omega_Lambda   =((double *)ADaPS_fetch(cosmo,"Omega_Lambda"))[0];
  Ez             =E_z(Omega_M,Omega_k,Omega_Lambda,redshift);
  one_plus_z_cube=(1.+redshift)*(1.+redshift)*(1.+redshift);
  return(Omega_M*one_plus_z_cube/(Ez*Ez));
}

// Compute Delta_vir(z) following the formula in 
//   Bryan & Norman (ApJ 495, 80, 1998)
double Delta_vir(double      redshift,
                 cosmo_info *cosmo){
  double x;
  double Omega;
  Omega=Omega_z(redshift,cosmo);
  x    =Omega-1.;
  return((18.*PI*PI+82*x-39*x*x)/Omega);
}

double Ha_Ho(double a, cosmo_info *cosmo){
  double Omega_M;
  double Omega_k;
  double Omega_Lambda;
  Omega_M     =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_M"))[0];
  Omega_k     =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_k"))[0];
  Omega_Lambda=((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_Lambda"))[0];
  return(E_z(Omega_M,Omega_k,Omega_Lambda,z_of_a(a))/
	 E_z(Omega_M,Omega_k,Omega_Lambda,0.));
}

double H_z(double redshift,cosmo_info *cosmo){
  double Omega_M;
  double Omega_k;
  double Omega_Lambda;
  double h_Hubble;
  Omega_M     =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_M"))[0];
  Omega_k     =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_k"))[0];
  Omega_Lambda=((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_Lambda"))[0];
  h_Hubble    =((double *)ADaPS_fetch((ADaPS *)cosmo,"h_Hubble"))[0];
//fprintf(stderr,"E(%10.3le)=%10.3le (%10.3le %10.3le %10.3le %10.3le)\n",redshift,E_z(Omega_M,Omega_k,Omega_Lambda,redshift),Omega_M,Omega_k,Omega_Lambda,h_Hubble);
  return(h_Hubble*1e2*E_z(Omega_M,Omega_k,Omega_Lambda,redshift));
}

double H_convert(double Hz){
  return(Hz*1e3/M_PER_MPC);
}

/* Returns dln(a)_dtau=da/(Ho*dtau) where tau is conformal time */
double dlna_dtau(double      a,
		 cosmo_info *cosmo){
  return(a*a*Ha_Ho(a,cosmo));
}

double rho_crit_z(double redshift, cosmo_info *cosmo){
  double Hz;
  Hz   =H_convert(H_z(redshift,cosmo));
  return(Hz*Hz/(2.*FOUR_THIRDS_PI*G_NEWTON));
}

double rho_crit_z_strip(double redshift,
                        double h_Hubble,
                        double Omega_M,
                        double Omega_Lambda){
  double Hz;
  Hz=H_convert(h_Hubble*1e2*E_z(Omega_M,(1.-Omega_M-Omega_Lambda),Omega_Lambda,redshift));
  return(Hz*Hz/(2.*FOUR_THIRDS_PI*G_NEWTON));
}

double D_Hubble(double h_Hubble){
  return(C_VACUUM*M_PER_MPC*1e-3/(100.0*h_Hubble));
}

double E_z(double Omega_M, double Omega_k, double Omega_Lambda, double z){
  double one_plus_z;
  double one_plus_z_sq;
  double one_plus_z_cu;
  double result;
  one_plus_z   =1.+z;
  one_plus_z_sq=one_plus_z*one_plus_z;
  one_plus_z_cu=one_plus_z_sq*one_plus_z;
  result=sqrt(Omega_M*one_plus_z_cu+Omega_k*one_plus_z_sq+Omega_Lambda);
  return(result);
}

double D_comove_integrand(double z, void *cosmo){
  double Omega_M;
  double Omega_k;
  double Omega_Lambda;
  double Ez;
  double r_val;
  Omega_M     =((double *)ADaPS_fetch((cosmo_info *)cosmo,"Omega_M"))[0];
  Omega_k     =((double *)ADaPS_fetch((cosmo_info *)cosmo,"Omega_k"))[0];
  Omega_Lambda=((double *)ADaPS_fetch((cosmo_info *)cosmo,"Omega_Lambda"))[0];
  if(z<=0.)
    r_val=1.;
  else{
    Ez   =E_z(Omega_M,Omega_k,Omega_Lambda,z);
    r_val=1.0/Ez;
  }
  return(r_val);
}
double D_comove(double z,cosmo_info *cosmo){
  int    n_int=1000;
  double h_Hubble;
  double rel_accuracy=1e-5;
  double abs_error;
  double limit_lo;
  double limit_hi;
  double integral;
  double r_val;
  gsl_integration_workspace *wspace;
  gsl_function               integrand;

  if (z<=0.)
    return(0.0);

  h_Hubble=((double *)ADaPS_fetch((ADaPS *)cosmo,"h_Hubble"))[0];

  // Initialize integral
  limit_lo          =0.;
  limit_hi          =z;
  integrand.function=D_comove_integrand;
  integrand.params  =(void *)cosmo;
  wspace            =gsl_integration_workspace_alloc(n_int);

  gsl_integration_qag(&integrand,
		      limit_lo,limit_hi,
		      0,rel_accuracy,
		      n_int,
		      GSL_INTEG_GAUSS61,
		      wspace,
		      &integral,&abs_error);
  r_val=D_Hubble(h_Hubble)*integral;

  // Clean-up
  gsl_integration_workspace_free(wspace);

  return(D_Hubble(h_Hubble)*integral);
}
double D_comove_transverse(double z,cosmo_info *cosmo){
  double D_c;
  double D_H;
  double D_c_t;
  double Omega_M;
  double Omega_k;
  double Omega_Lambda;
  double h_Hubble;
  Omega_M     =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_M"))[0];
  Omega_k     =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_k"))[0];
  Omega_Lambda=((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_Lambda"))[0];
  h_Hubble    =((double *)ADaPS_fetch((ADaPS *)cosmo,"h_Hubble"))[0];
  D_c=D_comove(z,cosmo);
  D_H=D_Hubble(h_Hubble);
  if(Omega_k>0.)
    D_c_t=D_H*sinh(sqrt(Omega_k)*D_c/D_H)/sqrt(Omega_k);
  else if(Omega_k==0.)
    D_c_t=D_c;
  else
    D_c_t=D_H*sin(sqrt(fabs(Omega_k))*D_c/D_H)/sqrt(fabs(Omega_k));
  return(D_c_t);
}

double D_angular(double z,cosmo_info *cosmo){
  return(D_comove_transverse(z,cosmo)/(1.+z));
}
double D_angular_1to2(double z_1,double z_2,cosmo_info *cosmo){
  double D_c_t_1;
  double D_c_t_2;
  double D_H;
  double h_Hubble;
  double Omega_k;
  h_Hubble=((double *)ADaPS_fetch((ADaPS *)cosmo,"h_Hubble"))[0];
  Omega_k =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_k"))[0];
  D_c_t_1 =D_comove_transverse(z_1,cosmo);
  D_c_t_2 =D_comove_transverse(z_2,cosmo);
  D_H     =D_Hubble(h_Hubble);
  return((D_c_t_2*sqrt(1.+Omega_k*pow(D_c_t_1/D_H,2.0))-
	  D_c_t_1*sqrt(1.+Omega_k*pow(D_c_t_2/D_H,2.0)))/(1.+z_2));
}
double D_luminosity(double z,cosmo_info *cosmo){
  return(D_comove_transverse(z,cosmo)*(1.+z));
}

void init_cosmo_std(cosmo_info **cosmo){
  double Omega_M;
  double Omega_k;
  double Omega_Lambda;
  double Omega_b;
  double f_gas;
  double h_Hubble;
  double sigma_8;
  double n_spectral;
  // WMAP-5 cosmology
  Omega_Lambda=0.727;
  Omega_M     =0.273;
  Omega_k     =0.;
  Omega_b     =0.0456;
  f_gas       =Omega_b/Omega_M;
  h_Hubble    =0.705;
  sigma_8     =0.812;
  n_spectral  =0.960;
  init_cosmo(cosmo,
	     Omega_Lambda,
	     Omega_M,
	     Omega_k,
             Omega_b,
             f_gas,
	     h_Hubble,
	     sigma_8,
	     n_spectral);
}

void init_cosmo(cosmo_info **cosmo,
		double       Omega_Lambda,
		double       Omega_M,
		double       Omega_k,
                double       Omega_b,
                double       f_gas,
		double       h_Hubble,
		double       sigma_8,
		double       n_spectral){
  ADaPS_init((ADaPS **)cosmo);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&Omega_Lambda),
             "Omega_Lambda",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&Omega_M),
             "Omega_M",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&Omega_k),
             "Omega_k",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&Omega_b),
             "Omega_b",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&f_gas),
             "f_gas",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&h_Hubble),
             "h_Hubble",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&sigma_8),
             "sigma_8",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&n_spectral),
             "n_spectral",
             ADaPS_SCALAR_DOUBLE);
}

void free_cosmo(cosmo_info **cosmo){
  ADaPS_free((ADaPS **)cosmo);
}
