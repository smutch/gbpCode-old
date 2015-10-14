#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double H_convert(double Hz){
  return(Hz*1e3/M_PER_MPC);
}

