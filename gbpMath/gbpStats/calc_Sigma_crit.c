#include <common.h>
#include <math.h>

double calc_Sigma_crit(double Da,
		       double beta_lens){
  return(C_VACUUM*C_VACUUM/(FOUR_PI*G_NEWTON*Da*beta_lens));
}
