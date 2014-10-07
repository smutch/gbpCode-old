#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double lk_sc(double       z,
             cosmo_info **cosmo,
             int          mode,
             int          component){
  double M_sc_val=M_sc(z,cosmo,mode,component);
  return(take_log10(k_of_M(M_sc_val,z,*cosmo)));
}

