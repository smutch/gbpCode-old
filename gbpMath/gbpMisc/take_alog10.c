#include <gbpCommon.h>
#include <math.h>

double take_alog10(double val){
  double rval;
  if(val>2.*LOG_ZERO)
    rval=(double)pow(10.,(double)val);
  else
    rval=0.;
  return((double)rval);
}
