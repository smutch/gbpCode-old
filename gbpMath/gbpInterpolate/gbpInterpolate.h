#ifndef GBPINTERPOLATE_AWAKE
#define GBPINTERPOLATE_AWAKE
#include <gbpLib.h>
#include <gsl/gsl_interp.h>

typedef struct interp_struct interp_info;
struct interp_struct{
  gsl_interp            *interp;
  gsl_interp_accel      *accel;
  size_t                 n;
  const gsl_interp_type *T;
  double                *x;
  double                *y;
};

void free_interpolate(void **interp);
void init_interpolate(double                 *x, 
	              double                 *y, 
	              size_t                  n,
                      const gsl_interp_type  *T,
	              interp_info           **interp);
void ADaPS_store_interp(ADaPS  **list,
                        void    *data,
                        char    *name,
                        ...);
double interpolate(interp_info *interp, 
		   double       x) ;
double interpolate_derivative(interp_info *interp,
	                      double       x);
double interpolate_integral(interp_info *interp,
	                    double       x_lo,
                            double       x_hi);
double interpolate_maximum_function(double x, void * params); 
void   interpolate_maximum(interp_info *interp,
			   double       x_lo_in,
			   double       x_guess_in,
			   double       x_hi_in,
			   double       threshold,
                           double      *x_maxima,
                           double      *y_maxima);
double interpolate_minimum_function(double x, void * params); 
void   interpolate_minimum(interp_info *interp,
			   double       x_lo_in,
			   double       x_guess_in,
			   double       x_hi_in,
			   double       threshold,
                           double      *x_minima,
                           double      *y_minima);

#endif

