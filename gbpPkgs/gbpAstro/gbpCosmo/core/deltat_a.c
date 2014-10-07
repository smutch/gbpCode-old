#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double deltat_a(cosmo_info **cosmo,double a_1,double a_2){
   double       a_lo;
   double       a_hi;
   interp_info *interp;
   if(!ADaPS_exist(*cosmo,"deltat_a_interp"))
      init_deltat_a(cosmo);
   interp=(interp_info *)ADaPS_fetch(*cosmo,"deltat_a_interp");
   a_lo=MAX(MIN(a_1,a_2),DELTAT_A_MIN_A);
   a_hi=MAX(a_1,a_2);
   return(interpolate_integral(interp,a_lo,a_hi));
}

