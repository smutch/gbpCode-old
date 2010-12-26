#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gbpInterpolate.h>

void init_interpolate(double                 *x, 
	              double                 *y, 
	              size_t                  n,
                      const gsl_interp_type  *T,
	              interp_info           **interp){
  size_t  i;

  (*interp)   =(interp_info *)malloc(sizeof(interp_info));
  (*interp)->n=n;
  (*interp)->x=(double *)malloc(sizeof(double)*n);
  (*interp)->y=(double *)malloc(sizeof(double)*n);
  (*interp)->T=T;

  // If the x-array is not in ascending order, switch it
  if(x[0]<x[n-1]){
    for(i=0;i<n;i++){
      (*interp)->x[i]=x[i];
      (*interp)->y[i]=y[i];
    }
  }
  else{
    for(i=0;i<n;i++){
      (*interp)->x[i]=x[n-1-i];
      (*interp)->y[i]=y[n-1-i];
    }
  }
  (*interp)->accel =gsl_interp_accel_alloc();
  (*interp)->interp=gsl_interp_alloc(T,n);
  gsl_interp_init((*interp)->interp,
                  (const double *)(*interp)->x,
                  (const double *)(*interp)->y,
                  (*interp)->n);
}

