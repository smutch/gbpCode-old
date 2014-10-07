#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double dDplus_da(double      a,
       cosmo_info *cosmo){
  double      h_Hubble;
  double      Ez;
  h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
  Ez      =H_z(z_of_a(a),cosmo)/(1e2*h_Hubble);
  return(1.0/pow(a*Ez,3.0));
}

