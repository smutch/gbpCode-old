#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gbpInterpolate.h>

double interpolate_derivative(interp_info *interp,
	                      double       x){
  double       r_val;
  r_val=gsl_interp_eval_deriv(((interp_info *)interp)->interp,
                              ((interp_info *)interp)->x,
                              ((interp_info *)interp)->y,
                              x,
                              ((interp_info *)interp)->accel);

  return(r_val);
}

