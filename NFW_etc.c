#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
//#include <gsl/gsl_sf_expint.h>

double Delta_vir_z(double z, cosmo_info *cosmo,int select){
  double Omega_M;
  double x;
  double Delta;
  switch(select){
    case DELTA_BRYAN_NORMAN: // Bryan & Norman '98
      Omega_M =((double *)ADaPS_fetch(cosmo,"Omega_M"))[0];
      x       =Omega_M-1.;
      Delta   =(18.*PI*PI+82.*x-39.*x*x)/Omega_M;
      break;
    case DELTA_LINEAR:
      Delta=1.686;
      break;
  }
  return(Delta);
}

void set_NFW_params(double      M,
		    double      z,
		    int         mode,
		    cosmo_info *cosmo,
		    double     *c_vir,
		    double     *R_vir){
  double M_o;
  double Delta,x;
  double h_Hubble;
  double Omega_M;
  switch(ADaPS_exist(cosmo,"M_WDM")){
  case FALSE:
    Omega_M =((double *)ADaPS_fetch(cosmo,"Omega_M"))[0];
    h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
    M_o     =M_sc(z,cosmo,PSPEC_LINEAR_TF,PSPEC_ALL_MATTER);
M_o=4.364e12*M_SOL/h_Hubble;
    (*c_vir)=(11./(1.+z))*pow(M/M_o,-0.13);     // Bullock et al '01 & Zehavi et al '04
    x       =Omega_M-1.;
    Delta   =Delta_vir_z(z,cosmo,DELTA_BRYAN_NORMAN);
Delta=200.;
    (*R_vir)=R_Delta_z(M,Delta,z,cosmo);        // Bullock et al '01
(*R_vir)=R_Delta_z(M,Delta,0.,cosmo);        // Bullock et al '01

    break;
  case TRUE:
    (*c_vir)=c_ENS(M,z,cosmo);          // Eke, Navarro and Steinmetz
    (*R_vir)=R_Delta_z(M,200.,z,cosmo); // R_200
    break;
  }
}

double R_vir_NFW(double      M_vir,
		 double      z,
		 int         mode,
		 cosmo_info *cosmo){
  double c_vir;
  double R_vir;
  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);
  return(R_vir);
}

double c_vir_NFW(double      M_vir,
		 double      z,
		 int         mode,
		 cosmo_info *cosmo){
  double c_vir;
  double R_vir;
  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);
  return(c_vir);
}

double rho_NFW(double      r,
	       double      M_vir,
	       double      z,
	       int         mode,
	       cosmo_info *cosmo){
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
double rho_NFW_fft(double      k,
		   double      M_vir,
		   double      z,
		   int         mode,
		   cosmo_info *cosmo){
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
double M_r_NFW(double      r,
	       double      M_vir,
	       double      z,
	       int         mode,
	       cosmo_info *cosmo){
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
double V_circ_NFW(double      r,
		  double      M_vir,
		  double      z,
		  int         mode,
		  cosmo_info *cosmo){
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

double V_circ_vir_NFW(double      M_vir,
		      double      z,
		      int         mode,
		      cosmo_info *cosmo){
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
double V_max_NFW(double      M_vir,
		 double      z,
		 int         mode,
		 cosmo_info *cosmo){
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

/* From Alam et al '02 */
double R_half_V_max_NFW(double      M_vir,
			double      z,
			int         mode,
			cosmo_info *cosmo){
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
double R_V_max_NFW(double      M_vir,
                   double      z,
                   int         mode,
                   cosmo_info *cosmo){
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
double Delta_half_V_max_NFW(double      M_vir,
			    double      z,
			    int         mode,
			    cosmo_info *cosmo){
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

