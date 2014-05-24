#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gsl/gsl_sf_expint.h>

void set_NFW_params(double       M,
		    double       z,
		    int          mode,
		    cosmo_info **cosmo,
		    double      *c_vir,
		    double      *R_vir){
  double M_o;
  double Delta;
  double h_Hubble;
  double Omega_M;

  if(mode!=NFW_MODE_DEFAULT)
     SID_trap_error("Unknown mode (%d) in set_NFW_params()",ERROR_LOGIC,mode);

  switch(ADaPS_exist(*cosmo,"M_WDM")){
  case FALSE:
    {
    Omega_M =((double *)ADaPS_fetch(*cosmo,"Omega_M"))[0];
    h_Hubble=((double *)ADaPS_fetch(*cosmo,"h_Hubble"))[0];
    M_o     =M_sc(z,cosmo,PSPEC_LINEAR_TF,PSPEC_ALL_MATTER);

    // Mass-concentration from Munoz-Cuartas et al 2010
    //(*c_vir)=(11./(1.+z))*pow(M/M_o,-0.13);     // Bullock et al '01 & Zehavi et al '04
    double w    =   0.029;
    double m    =   0.097;
    double alpha=-110.001;
    double beta =2469.720;
    double gamma=  16.885;
    double a_z  =w*z-m;
    double b_z  =alpha/(z+gamma)+beta/pow(z+gamma,2.);
    (*c_vir)    =take_alog10(a_z*take_log10(M/(M_SOL/h_Hubble))+b_z);

    Delta   =Delta_vir(z,*cosmo);
    Delta=200.;
    (*R_vir)=R_Delta_z(M,Delta,z,*cosmo);       // Bullock et al '01
    }
    break;
  case TRUE:
    SID_trap_error("ENS not working.",ERROR_LOGIC);
    //(*c_vir)=c_ENS(M,z,*cosmo);          // Eke, Navarro and Steinmetz
    //(*R_vir)=R_Delta_z(M,200.,z,*cosmo); // R_200
    break;
  }
}

double R_vir_NFW(double       M_vir,
		 double       z,
		 int          mode,
		 cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);
  return(R_vir);
}

double c_vir_NFW(double       M_vir,
		 double       z,
		 int          mode,
		 cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);
  return(c_vir);
}

double rho_NFW(double       r,
	       double       M_vir,
	       double       z,
	       int          mode,
	       cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double g_c;
  double rho_o;
  double x;
  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);
  g_c  =1./(log(1.+c_vir)-c_vir/(1.+c_vir));
  rho_o=M_vir*g_c/(FOUR_PI*pow(R_vir/c_vir,3.));
  x    =r*c_vir/R_vir;
  return(rho_o/(x*pow(1.+x,2.)));
}

/* FFT of NFW profile from White '01 */
double rho_NFW_fft(double       k,
		   double       M_vir,
		   double       z,
		   int          mode,
		   cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double r_s;
  double z_k;
  double g_c;
  double rho_o;
  double x1;
  double x2;
  double Ci1;
  double Ci2;
  double Si1;
  double Si2;
  double r_val;
  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);
  r_s  =R_vir/c_vir;
  z_k  =k*r_s;
  g_c  =1./(log(1.+c_vir)-c_vir/(1.+c_vir));
  rho_o=M_vir*g_c/(FOUR_PI*r_s*r_s*r_s);
  x1   =1.+c_vir;
  x2   =c_vir*z_k;
  Ci1=gsl_sf_Ci((1.+c_vir)*z_k);
  Si1=gsl_sf_Si((1.+c_vir)*z_k);
  Ci2=gsl_sf_Ci(           z_k);
  Si2=gsl_sf_Si(           z_k);
  r_val=(FOUR_PI*rho_o*r_s*r_s*r_s/M_vir)*(cos(z_k)*(Ci1-Ci2)+
				           sin(z_k)*(Si1-Si2)-
				           sin(c_vir*z_k)/(z_k*(1.+c_vir)));
  return(r_val);
}
 
/**************************************/
/* Cole and Lacey (1996) M(r) profile */
/**************************************/
double M_r_NFW(double       r,
	       double       M_vir,
	       double       z,
	       int          mode,
	       cosmo_info **cosmo){
  double r_o;
  double x;
  double g_c;
  double c_vir;
  double R_vir;
  double r_val=0.;

  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);
  r_o=R_vir/c_vir;
  x=r/r_o;
  if(x>0.){
    g_c  =1./(log(1.+c_vir)-c_vir/(1.+c_vir));
    r_val=M_vir*g_c*(log(1.+x)-x/(1.+x));
  }
  return(r_val);
}

/****************************************/
/* Cole and Lacey (1996) V_c(r) profile */
/****************************************/
double V_circ_NFW(double       r,
		  double       M_vir,
		  double       z,
		  int          mode,
		  cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double V2_c_characteristic;
  double g_c;
  double r_o;
  double x;
  double r_val=0.;

  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);

  r_o=R_vir/c_vir;
  x=r/r_o;
  if(x>0.){
    V2_c_characteristic=G_NEWTON*M_vir/R_vir;
    g_c=1./(log(1.+c_vir)-c_vir/(1.+c_vir));
    r_val=sqrt(V2_c_characteristic*g_c*(log(1.+x)-x/(1.+x))*c_vir/x);
  }
  return(r_val);
}

double V_circ_vir_NFW(double       M_vir,
		      double       z,
		      int          mode,
		      cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double r;
  double V2_c_characteristic;
  double g_c;
  double r_o;
  double x;
  double r_val=0.;

  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);

  r_val=V_circ_NFW(R_vir,M_vir,z,mode,cosmo);

  return(r_val);
}

/* From Alam et al '02 */
double V_max_NFW(double       M_vir,
		 double       z,
		 int          mode,
		 cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double V2_vir;
  double g_c;
  double r_val=0.;

  if(M_vir>0.){
    set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);
    
    V2_vir=V2_circ(M_vir,R_vir);
    g_c   =1./(log(1.+c_vir)-c_vir/(1.+c_vir));
    r_val =sqrt(0.2*V2_vir*c_vir*g_c);
  }

  return(r_val);
}

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

/* From Alam et al '02 */
double R_half_V_max_NFW(double       M_vir,
			double       z,
			int          mode,
			cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double r_o;
  double r_val=0.;

  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);

  r_o  =R_vir/c_vir;  
  r_val=0.13*r_o;

  return(r_val);
}

// Approximate R(V=V_max) from NFW '97
double R_V_max_NFW(double       M_vir,
                   double       z,
                   int          mode,
                   cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double r_o;
  double r_val=0.;

  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);

  r_o  =R_vir/c_vir;
  r_val=2.*r_o;

  return(r_val);
}

/* From Alam et al '02 */
double Delta_half_V_max_NFW(double       M_vir,
			    double       z,
			    int          mode,
			    cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double g_c;
  double r_val=0.;

  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);

  g_c  =1./(log(1.+c_vir)-c_vir/(1.+c_vir));
  r_val=409.*c_vir*c_vir*c_vir*g_c;

  return(r_val);
}

double V2_circ(double M, 
	       double r){
  return(G_NEWTON*M/r);
}

