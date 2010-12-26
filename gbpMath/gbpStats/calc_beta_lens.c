/* Needed for computing the lensing Sigma_crit */
#include <common.h>
#include <math.h>

double calc_beta_lens(double D_lens,
		      double D_source){
  return(1.0-D_lens/D_source);
}
