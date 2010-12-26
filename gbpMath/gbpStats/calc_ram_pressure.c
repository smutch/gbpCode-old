#include <common.h>
#include <math.h>

double calc_ram_pressure(double rho_gas,
                         double v_gas){
  return(rho_gas*v_gas*v_gas);
}
