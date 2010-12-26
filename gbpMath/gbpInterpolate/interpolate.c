#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gbpInterpolate.h>

double interpolate(interp_info *interp, 
		   double       x) {
  return(gsl_interp_eval(interp->interp,
			 interp->x,
			 interp->y,
			 x,
			 interp->accel));
}

