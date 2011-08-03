#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gbpInterpolate.h>

void free_interpolate(void **interp){
  if((*interp)!=NULL){
    free(((interp_info *)(*interp))->x);
    free(((interp_info *)(*interp))->y);
    gsl_interp_free(((interp_info *)(*interp))->interp);
    gsl_interp_accel_free(((interp_info *)(*interp))->accel);
    free(*interp);
  }
}

