/********************************************************/
/* Generate a random field with a Gaussian distribution */
/********************************************************/
#include <math.h>
#include <gbpCommon.h>
#include <gbpRNG.h>

REAL random_gaussian(RNG_info *RNG){
  double rsq,fac,v1,v2,rval;
  double rand1,rand2;
  if((RNG->IGauss)<=0){
    rsq=0.;
    while(rsq>=1. || rsq<=0.){
      rand1=(double)random_number(RNG);
      rand2=(double)random_number(RNG);
      v1=2.*rand1-1.;
      v2=2.*rand2-1.;
      rsq=v1*v1+v2*v2;
    }
    fac=sqrt(-2.*log(rsq)/rsq);
    rval           =v1*fac;
    (RNG->GaussBak)=v2*fac;
    (RNG->IGauss)  =1;
  }
  else{
    rval       =RNG->GaussBak;
    RNG->IGauss=0;
  }
  return((REAL)rval);
}
