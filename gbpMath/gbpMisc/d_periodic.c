#include <gbpLib.h>
#include <math.h>

double d_periodic(double d,double box_size){
  double d_min;
  d_min=d;
  if(fabs(d-box_size)<fabs(d_min))
    d_min=d-box_size;
  if(fabs(d+box_size)<fabs(d_min))
    d_min=d+box_size;
  return(d_min);
}

