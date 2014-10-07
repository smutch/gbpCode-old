#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

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

