#include <common.h>
#include <math.h>

double calc_Bondi(double rho_gas,
                  double T_gas,
                  double M_9){  // M_9 is the black hole mass; units of 10^9M_sol
  // Formula from Nulsen's 2003 Conf. Proc. -- Output is in units of M_sol/yr
  return(0.012*rho_gas*1e-6*NE_PER_RHOGAS*M_9*M_9/pow(T_gas/K_PER_KEV,THREE_HALVES));
}
