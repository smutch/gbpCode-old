#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <gbpLib.h>
#include <gbpMisc.h>

void apply_rotation(double  x_hat,
                    double  y_hat,
                    double  z_hat,
                    double  theta,
                    double *x_i,
                    double *y_i,
                    double *z_i){
  double X;
  double Y;
  double Z;
  double W;
  double x_p;
  double y_p;
  double z_p;
  X     =x_hat*sin(theta/2);
  Y     =y_hat*sin(theta/2);
  Z     =z_hat*sin(theta/2);
  W     =      cos(theta/2);
  x_p   =(double)(*x_i);
  y_p   =(double)(*y_i);
  z_p   =(double)(*z_i);
  (*x_i)=(1.-2.*Y*Y-2.*Z*Z)*x_p+(   2.*X*Y-2.*Z*W)*y_p+(   2.*X*Z+2.*Y*W)*z_p;
  (*y_i)=(   2.*X*Y+2.*Z*W)*x_p+(1.-2.*X*X-2.*Z*Z)*y_p+(   2.*Y*Z-2.*X*W)*z_p;
  (*z_i)=(   2.*X*Z-2.*Y*W)*x_p+(   2.*Y*Z+2.*X*W)*y_p+(1.-2.*X*X-2.*Y*Y)*z_p;
}

