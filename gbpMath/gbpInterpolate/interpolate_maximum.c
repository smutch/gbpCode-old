#include <stdio.h>
#include <math.h>
#include <gbpLib.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h> 
#include <gsl/gsl_min.h> 
#include <gbpInterpolate.h>

double interpolate_maximum_function(double x, void * params);
double interpolate_maximum_function(double x, void * params){ 
  return(-interpolate((interp_info *)params,x)); 
} 

void   interpolate_maximum(interp_info *interp,
			   double       x_lo_in,
			   double       x_guess_in,
			   double       x_hi_in,
			   double       threshold,
                           double      *x_maxima,
                           double      *y_maxima){
  double                   x_lo;
  double                   x_hi;
  double                   x_guess;
  double                   r_val;
  const gsl_min_fminimizer_type *T; 
  gsl_min_fminimizer      *s;
  gsl_function             F; 
  int                      status; 
  int                      iter=0;
  int                      max_iter; 

  x_lo      =x_lo_in;
  x_guess   =x_guess_in;
  x_hi      =x_hi_in;
  max_iter  =MAX(100,interp->n);
  F.function=interpolate_maximum_function;
  F.params  =(void *)interp;
  T         =gsl_min_fminimizer_brent; 
  s         =gsl_min_fminimizer_alloc(T); 
  gsl_min_fminimizer_set(s,&F,x_guess,x_lo,x_hi);

  do { 
    iter++; 
    status =gsl_min_fminimizer_iterate(s); 
    x_lo   =gsl_min_fminimizer_x_lower(s); 
    x_guess=gsl_min_fminimizer_x_minimum(s); 
    x_hi   =gsl_min_fminimizer_x_upper(s); 
    status =gsl_min_test_interval(x_lo,x_hi,0.,threshold); 
  } 
  while(status==GSL_CONTINUE && iter<max_iter);
  gsl_min_fminimizer_free(s);
  (*x_maxima)=x_guess;
  (*y_maxima)=interpolate(interp,x_guess);
}

