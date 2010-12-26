#include <common.h>
#include <math.h>

double calc_Bondi_radius(double T_gas,
                         double M_9){  // M_9 is the black hole mass; units of 10^9M_sol
  // Formula from Allen et al 2006
  return(2.0*G_NEWTON*M_9*1e9*M_SOL/pow(calc_sound_speed(GAMMA_ICM,T_gas),2.0));
}
