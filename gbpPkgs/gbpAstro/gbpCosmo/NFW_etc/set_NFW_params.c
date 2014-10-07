#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>
#include <gbpCosmo_NFW_etc.h>
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

