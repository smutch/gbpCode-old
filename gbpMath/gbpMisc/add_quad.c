#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <gbpLib.h>
double add_quad(int n_d, ...){
  int      i;
  double   t_val;
  double   r_val;
  va_list  vargs;
  va_start(vargs,n_d);
  r_val=0.;
  for(i=0;i<n_d;i++){
    t_val =(double)va_arg(vargs,double);
    r_val+=t_val*t_val;
  }
  va_end(vargs);
  return(sqrt(r_val));
}

