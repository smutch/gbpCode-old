#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gbpInterpolate.h>

// params is not used at the moment but we must allow for
//   it to meet the required ADaPS function definition
void free_interpolate(void **interp,void *params){
  if((*interp)!=NULL){
    SID_free(SID_FARG ((interp_info *)(*interp))->x);
    SID_free(SID_FARG ((interp_info *)(*interp))->y);
    gsl_interp_free(((interp_info *)(*interp))->interp);
    gsl_interp_accel_free(((interp_info *)(*interp))->accel);
    SID_free(SID_FARG *interp);
  }
}

