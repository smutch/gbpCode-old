#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double R_NL_z(double z,cosmo_info **cosmo){
   int    n_iter   =0;
   int    max_iter =100;
   double threshold=1e-3;
   double R_lo     =1e-1*M_PER_MPC;
   double R_hi     =50.*M_PER_MPC;
   double R_mid    =0.5*(R_lo+R_hi);
   int    mode     =PSPEC_LINEAR_TF;
   int    component=PSPEC_ALL_MATTER;
   double target   =1.;
   double test_min =sigma_R(cosmo,R_hi, z,mode,component);
   double test_max =sigma_R(cosmo,R_lo, z,mode,component);
   double test     =sigma_R(cosmo,R_mid,z,mode,component);
   // Make sure the bracketing range covers the situation
   while(test_min>target){
      R_hi    *=2.;
      test_min =sigma_R(cosmo,R_hi, z,mode,component);
   }
   while(test_max<target){
      R_lo    *=0.5;
      test_max =sigma_R(cosmo,R_lo, z,mode,component);
   }
   // Perform bisection
   double dR=R_hi-R_lo;
   while(fabs(test/target)>threshold && (dR/R_mid)>threshold && n_iter<max_iter){
      if(test<target)
         R_hi=R_mid;
      else
         R_lo=R_mid;
      R_mid=0.5*(R_lo+R_hi);
      dR   =R_hi-R_lo;
      test =sigma_R(cosmo,R_mid,z,mode,component);
      n_iter++;
   }
   if(n_iter==max_iter)
      SID_trap_error("Maximum number of iterations (%d) reached in R_gbpCosmo2gbpCosmo.",ERROR_LOGIC,n_iter);

   return(R_mid);
}

