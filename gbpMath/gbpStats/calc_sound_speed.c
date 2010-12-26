#include <common.h>
#include <math.h>

double calc_sound_speed(double gamma,
                        double T_gas){
  return(sqrt(gamma*K_BOLTZMANN*T_gas/(MU_MMW*M_PROTON)));
}
