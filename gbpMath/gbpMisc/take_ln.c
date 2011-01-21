#include <math.h>
#include <gbpLib.h>

double take_ln(double val){
  double rval;
  if(val>0.0)
    rval=(double)log((double)val);
  else
    rval=(double)LOG_ZERO;
  return((double)rval);
}
