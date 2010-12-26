#include <common.h>
#include <math.h>

double calc_tcool(double rho_gas,
                  double T_gas){
  return(MU_MMW*M_PROTON*THREE_HALVES*K_BOLTZMANN*T_gas/(rho_gas*Lambda_MB(T_gas,0.3)));
}
