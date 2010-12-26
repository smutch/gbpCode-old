#include <common.h>
#include <math.h>

double calc_entropy(double rho_gas,
                    double T_gas){
  double rval;
  if(rho_gas>0.)
    rval=(T_gas/K_PER_KEV)/pow(rho_gas*NE_PER_RHOGAS*1e-6,TWO_THIRDS);
  else
    rval=0.;
  return(rval);
}
