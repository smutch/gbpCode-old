#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

void init_cosmo_std(cosmo_info **cosmo){
   read_gbpCosmo_file(cosmo,"default");
}

