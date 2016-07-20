#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_mass_functions.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// Scaled mass functions
double scaled_mass_function(double sigma,int mode,double *P){
  double  rval;
  if(check_mode_for_flag(mode,MF_WATSON)){ // Watson et al. (2013) universal FoF mass function
    double A     = 0.282;
    double alpha = 2.163;
    double beta  = 1.406;
    double gamma = 1.210;
    if(check_mode_for_flag(mode,MF_PASS_PARAMS)){
       A     = P[0];
       alpha = P[1];
       beta  = P[2];
       gamma = P[3];
    }
    rval=A*(pow(beta/sigma,alpha)+1.)*exp(-gamma/(sigma*sigma));
  }
  else if(check_mode_for_flag(mode,MF_TIAMAT)){ // Poole et al. (2013) universal FoF mass function for Tiamat
    double A     =  0.03331;
    double alpha =  1.153;
    double beta  = 12.33;
    double gamma =  1.009;
    if(check_mode_for_flag(mode,MF_PASS_PARAMS)){
       A     = P[0];
       alpha = P[1];
       beta  = P[2];
       gamma = P[3];
    }
    rval=A*(pow(beta/sigma,alpha)+1.)*exp(-gamma/(sigma*sigma));
  }
  else if(check_mode_for_flag(mode,MF_JENKINS)){ // Jenkins et al. (2001)
    double A=0.315;
    double B=0.61;
    double C=3.8;
    if(check_mode_for_flag(mode,MF_PASS_PARAMS)){
       A=P[0];
       B=P[1];
       C=P[2];
    }
    rval =A*exp(-pow((double)fabs((float)(take_ln(1./sigma)+B)),C));
  }
  else if(check_mode_for_flag(mode,MF_ST)){ // Sheth-Torman (1999)
    double delta_k=1.686;
    double A      =0.3222;
    double B      =0.707;
    double C      =0.3;
    if(check_mode_for_flag(mode,MF_PASS_PARAMS)){
       A      =P[0];
       B      =P[1];
       C      =P[2];
    }
    rval =A*sqrt(2.*B/PI)*(delta_k/sigma)*exp(-B*delta_k*delta_k/(2.*sigma*sigma))*(1.+pow(sigma*sigma/(B*delta_k*delta_k),C));
  }
  else if(check_mode_for_flag(mode,MF_PS)){ // Press-Schechter (1974)
    double delta_k=1.686;
    if(check_mode_for_flag(mode,MF_PASS_PARAMS)){
       delta_k=P[0];
    }
    rval =sqrt(2./PI)*(delta_k/sigma)*exp(-delta_k*delta_k/(2.*sigma*sigma));
  }
  else
    SID_trap_error("A valid mass function was not specified with mode (%d) in scaled_mass_function().\n",ERROR_LOGIC,mode);
  return(rval);
}

