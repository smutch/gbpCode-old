#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_mass_functions.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// Scaled mass functions
double scaled_mass_function(double sigma,int select_flag){
  double delta_k;
  double rval;
  switch(select_flag){
  case MF_WATSON:{ // Watson et al. (2013) universal FoF mass function
    double A     = 0.282;
    double alpha = 2.163;
    double beta  = 1.406;
    double gamma = 1.210;
    rval         = A*(pow(beta/sigma,alpha)+1.)*exp(-gamma/(sigma*sigma));
    break;
    }
  case MF_JENKINS: // Jenkins et al. (2001)
    delta_k=1.686;
    rval   =0.315*exp(-pow((double)fabs((float)(take_ln(1./sigma)+0.61)),3.8));
    break;
  case MF_PS: // Press-Schechter (1974)
    delta_k=1.686;
    rval   =sqrt(2./PI)*(delta_k/sigma)*exp(-delta_k*delta_k/(2.*sigma*sigma));
    break;
  case MF_ST: // Sheth-Torman (1999)
    delta_k=1.686;
    rval   =0.3222*sqrt(2.*0.707/PI)*(delta_k/sigma)*exp(-0.707*delta_k*delta_k/(2.*sigma*sigma))*
      (1.+pow(sigma*sigma/(0.707*delta_k*delta_k),0.3));
    break;
  default:
    fprintf(stderr,"ERROR in scaled_mass_function: select_flag=%d unknown.\n",select_flag);
    rval=0.;
    break;
  }
  return rval;
}

