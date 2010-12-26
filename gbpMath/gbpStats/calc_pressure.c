#include <common.h>
#include <math.h>

double calc_pressure(double rho_gas,
                     double T_gas){
  return(K_BOLTZMANN*T_gas*rho_gas/(MU_MMW*M_PROTON));
}
