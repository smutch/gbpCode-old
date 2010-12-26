#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h> 
#include <gbpLib.h>
#include <gbpInterpolate.h>

double interpolate_integral(interp_info *interp,
	                    double       x_lo,
                            double       x_hi){
  double x_lo_interp;
  double x_hi_interp;
  double x_lo_tmp;
  double x_hi_tmp;
  double x_min;
  double x_max;
  double y_x_min;
  double y_x_max;
  double r_val;
  double alpha;
  int    flag_error=FALSE;

  x_min   =((interp_info *)interp)->x[0];
  x_max   =((interp_info *)interp)->x[((interp_info *)interp)->n-1];
  y_x_min =((interp_info *)interp)->y[0];
  y_x_max =((interp_info *)interp)->y[((interp_info *)interp)->n-1];

  r_val      =0.;
  x_lo_interp=x_lo;
  x_hi_interp=x_hi;

  if(x_hi>x_lo){
  // Extrapolate integral to low-x assuming dI prop. to x^alpha
  if(x_lo<x_min){
    x_lo_tmp=x_lo;
    x_hi_tmp=MIN(x_hi,x_min);
    alpha=interpolate_derivative(interp,((interp_info *)interp)->x[0]);
    if(alpha<0. && x_lo_tmp*x_hi_tmp==0.){
      fprintf(stderr,"ERROR: integration extrapolation to low-x is divergent! (alpha=%lf)\n",alpha);
      flag_error=TRUE;
    }
    else
      r_val+=(y_x_min/(alpha+1.))*(pow(x_hi_tmp/x_min,alpha+1.)-pow(x_lo_tmp/x_min,alpha+1.));
    x_lo_interp=x_min;
  }

  // Extrapolate integral to high-x assuming dI prop. to x^alpha
  if(x_hi>x_max){
    x_lo_tmp=MAX(x_max,x_lo);
    x_hi_tmp=x_hi;
    alpha=interpolate_derivative(interp,((interp_info *)interp)->x[((interp_info *)interp)->n-1]);
    if(alpha<0. && x_lo_tmp*x_hi_tmp==0.){
      fprintf(stderr,"ERROR: integration extrapolation to high-x is divergent! (alpha=%lf)\n",alpha);
      flag_error=TRUE;
    }
    else
      r_val+=(y_x_min/(alpha+1.))*(pow(x_hi_tmp/x_min,alpha+1.)-pow(x_lo_tmp/x_min,alpha+1.));
    x_hi_interp=x_max;
  }

  r_val+=gsl_interp_eval_integ(((interp_info *)interp)->interp,
                               ((interp_info *)interp)->x,
                               ((interp_info *)interp)->y,
                               x_lo_interp,x_hi_interp,
                               ((interp_info *)interp)->accel);
  }
  return(r_val);
}

