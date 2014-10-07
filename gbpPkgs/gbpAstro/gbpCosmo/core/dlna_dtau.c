#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

/* Returns dln(a)_dtau=da/(Ho*dtau) where tau is conformal time */
double dlna_dtau(double      a,
                 cosmo_info *cosmo){
  return(a*a*Ha_Ho(a,cosmo));
}

