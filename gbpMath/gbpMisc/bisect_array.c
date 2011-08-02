#include <stdio.h>
#include <stdlib.h>
#include <gbpLib.h>
#include <gbpMisc.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

double bisect_array_function(double  x_interp,
			     void   *params){
  double interp;
  double r_val;
  interp=
    interpolate(((bisect_af_params *)params)->interp,x_interp);
  r_val=interp-((bisect_af_params *)params)->value;
  return(r_val);
}

double bisect_array(interp_info *interp, 
		    double       value,
		    double       threshold) {
  double *x;
  double *y;
  double *d2y;
  double  x_lo;
  double  x_hi;
  double  r_val;
  int     i;
  int     alloc_flag=FALSE;
  int     iter      =0;
  int     max_iter  =500;
  int     status    =GSL_CONTINUE;
  bisect_af_params             params;
  gsl_function                 F;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver            *s;

  params.interp=interp;
  params.value =value;
  F.function   =bisect_array_function;
  F.params     =(void *)(&params);

  x_lo=interp->x[0];
  x_hi=interp->x[interp->n-1];

  if(fabs(bisect_array_function(x_lo,&params))<threshold)
    return(x_lo);
  if(fabs(bisect_array_function(x_hi,&params))<threshold)
    return(x_hi);

  T=gsl_root_fsolver_brent;
  s=gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s,
		       &F,
		       x_lo,
		       x_hi);

  while(status==GSL_CONTINUE && iter < max_iter){
    status=gsl_root_fsolver_iterate(s);
    r_val =gsl_root_fsolver_root(s);
    x_lo  =gsl_root_fsolver_x_lower(s);
    x_hi  =gsl_root_fsolver_x_upper(s);
    status=gsl_root_test_interval(x_lo, 
				  x_hi,
				  0,
				  threshold);
    iter++;
  }
     
  // Clean-up
  gsl_root_fsolver_free(s);

  return(r_val);
}

