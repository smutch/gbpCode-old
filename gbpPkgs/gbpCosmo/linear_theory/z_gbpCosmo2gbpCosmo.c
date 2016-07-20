#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double z_gbpCosmo2gbpCosmo(double z,gbpCosmo2gbpCosmo_info *gbpCosmo2gbpCosmo){
  double z_prime=z;

  // Perform scaling if it is given
  if(gbpCosmo2gbpCosmo!=NULL){
     // Sanity check
     if(z<gbpCosmo2gbpCosmo->z_min)
        SID_trap_error("A scaled redshift less than the given minimum has been requested (i.e. %le<%le).  Reinitialize gbpCosmo2gbpCosmo with a lower redshift.",ERROR_LOGIC,z,gbpCosmo2gbpCosmo->z_min);
     // Calculate scaling
     int    n_iter   =0;
     int    max_iter =100;
     double threshold=1e-3;
     double z_lo     =0.;
     double z_hi     =MAX(1e2,5.*z);
     double z_mid    =0.5*(z_lo+z_hi);
     double target   =linear_growth_factor(z,    gbpCosmo2gbpCosmo->cosmo_source)*gbpCosmo2gbpCosmo->D_ratio;
     double test_min =linear_growth_factor(z_hi, gbpCosmo2gbpCosmo->cosmo_target);
     double test_max =linear_growth_factor(z_lo, gbpCosmo2gbpCosmo->cosmo_target);
     double test     =linear_growth_factor(z_mid,gbpCosmo2gbpCosmo->cosmo_target);
     // Make sure the bracketing range covers the situation
     while(test<test_min){
        z_hi    *=2.;
        test_min =linear_growth_factor(z_hi,gbpCosmo2gbpCosmo->cosmo_target);
     }
     while(test>test_max){
        z_lo    *=0.5;
        test_min =linear_growth_factor(z_lo,gbpCosmo2gbpCosmo->cosmo_target);
     }
     // Perform bisection
     double dz=z_hi-z_lo;
     while(fabs(test/target)>threshold && (dz/z_mid)>threshold && n_iter<max_iter){
        if(test<target)
           z_hi=z_mid;
        else
           z_lo=z_mid;
        z_mid=0.5*(z_lo+z_hi);
        dz   =z_hi-z_lo;
        test =linear_growth_factor(z_mid,gbpCosmo2gbpCosmo->cosmo_target);
        n_iter++;
     }
     if(n_iter==max_iter)
        SID_trap_error("Maximum number of iterations (%d) reached in z_gbpCosmo2gbpCosmo.",ERROR_LOGIC,n_iter);
     z_prime=z_mid;
  }

  return(z_prime);
}

