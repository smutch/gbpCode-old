#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

void init_cosmo(cosmo_info **cosmo,
      double       Omega_Lambda,
      double       Omega_M,
      double       Omega_k,
      double       Omega_b,
      double       f_gas,
      double       h_Hubble,
      double       sigma_8,
      double       n_spectral){
  ADaPS_init((ADaPS **)cosmo);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&Omega_Lambda),
             "Omega_Lambda",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&Omega_M),
             "Omega_M",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&Omega_k),
             "Omega_k",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&Omega_b),
             "Omega_b",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&f_gas),
             "f_gas",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&h_Hubble),
             "h_Hubble",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&sigma_8),
             "sigma_8",
             ADaPS_SCALAR_DOUBLE);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&n_spectral),
             "n_spectral",
             ADaPS_SCALAR_DOUBLE);
}

