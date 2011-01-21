#include <math.h>
#include <gbpLib.h>

double take_aln(double val){
  double rval;
  if(val>2.*LOG_ZERO)
    rval=(double)exp((double)val);
  else
    rval=0.;
  return((double)rval);
}
