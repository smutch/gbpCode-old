#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>
#include <gbpCosmo_NFW_etc.h>
#include <gsl/gsl_sf_expint.h>

double R_vir_NFW(double       M_vir,
                 double       z,
                 int          mode,
                 cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);
  return(R_vir);
}

