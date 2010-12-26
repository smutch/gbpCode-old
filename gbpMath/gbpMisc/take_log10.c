#include <gbpCommon.h>
#include <math.h>

double take_log10(double val){
  double rval;
  if(val>0.0)
    rval=(double)log10((double)val);
  else
    rval=(double)LOG_ZERO;
  return((double)rval);
}
