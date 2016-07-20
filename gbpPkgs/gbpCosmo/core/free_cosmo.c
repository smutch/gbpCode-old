#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

void free_cosmo(cosmo_info **cosmo){
  ADaPS_free((void **) cosmo);
}

