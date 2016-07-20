#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

void init_deltat_a(cosmo_info **cosmo){
  int     n_int;
  int     i_int;
  double *x_int;
  double *y_int;
  double  a_min;
  double  a_max;
  double  da;
  double  a;
  interp_info *interp;
  n_int=250;
  x_int=(double *)SID_malloc(sizeof(double)*n_int);
  y_int=(double *)SID_malloc(sizeof(double)*n_int);
  a_min=DELTAT_A_MIN_A;
  a_max=1.;
  da   =(a_max-a_min)/((double)(n_int-2));
  x_int[0]=0.;
  y_int[0]=0.;
  for(i_int=1,a=a_min;i_int<(n_int-1);i_int++,a+=da){
    x_int[i_int]=a;
    y_int[i_int]=1./(a*H_convert(H_z(z_of_a(a),(*cosmo))));
  }
  x_int[i_int]=a_max;
  y_int[i_int]=1./(a_max*H_convert(H_z(z_of_a(a_max),(*cosmo))));
  init_interpolate(x_int,y_int,(size_t)n_int,gsl_interp_cspline,&interp);
  SID_free(SID_FARG x_int);
  SID_free(SID_FARG y_int);
  ADaPS_store_interp(cosmo,
                     (void *)(interp),
                     "deltat_a_interp");
}

