/*********************************************************/
/* Generate a random field with a lognormal distribution */
/*********************************************************/
#include <math.h>
#include <gbpLib.h>
#include <gbpMisc.h>
#include <gbpRNG.h>

double calc_pdf_local(double x,double mu,double sigma);
double calc_pdf_local(double x,double mu,double sigma){
  double term1,term2,term3;
  term1=1./(x*sqrt(TWO_PI)*sigma);
  term2=pow(take_ln(x/mu),2.);
  term3=2.*sigma*sigma;
  return(term1*exp(-term2/term3));
}

REAL random_lognormal(RNG_info *RNG,double mu,double sigma){
  double range1,range2;
  double rand1,rand2;
  int    flag_continue=TRUE;

  // Loop until we have sucessfully generated a random number
  range1=10.*mu*fabs(exp(sigma));    // Generate values out to 10-sigma
  range2=1./(sigma*sqrt(PI)*mu); 
  while(flag_continue){
    rand1=range1*(double)random_number(RNG); // Generate random x-value
    rand2=range2*(double)random_number(RNG); // Generate random y-value
    if(rand2<calc_pdf_local(rand1,mu,sigma)) flag_continue=FALSE;
  }
  return((REAL)rand1);
}
