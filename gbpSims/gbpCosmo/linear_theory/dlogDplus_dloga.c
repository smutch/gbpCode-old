#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

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

