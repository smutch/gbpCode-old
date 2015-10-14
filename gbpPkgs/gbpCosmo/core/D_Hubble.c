#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double D_Hubble(double h_Hubble){
  return(C_VACUUM*M_PER_MPC*1e-3/(100.0*h_Hubble));
}

