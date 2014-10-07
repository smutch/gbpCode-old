#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>
#include <gbpCosmo_NFW_etc.h>
#include <gsl/gsl_sf_expint.h>

double V2_circ(double M, 
               double r){
  return(G_NEWTON*M/r);
}

