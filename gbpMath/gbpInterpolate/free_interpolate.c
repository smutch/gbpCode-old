#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gbpInterpolate.h>

void free_interpolate(interp_info **interp){
  if((*interp)!=NULL){
    free((*interp)->x);
    free((*interp)->y);
    gsl_interp_free((*interp)->interp);
    gsl_interp_accel_free((*interp)->accel);
    free(*interp);
  }
}

