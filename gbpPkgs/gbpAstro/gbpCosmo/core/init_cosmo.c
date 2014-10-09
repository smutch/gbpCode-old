#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

void init_cosmo(cosmo_info **cosmo,
                const char  *name,
                double       Omega_Lambda,
                double       Omega_M,
                double       Omega_k,
                double       Omega_b,
                double       f_gas,
                double       h_Hubble,
                double       sigma_8,
                double       n_spectral){

  // Print a log message
  char name_store[32];
  if(name!=NULL)
     strcpy(name_store,name);
  else
     sprintf(name_store,"unnamed");
  SID_log("Initializing %s cosmology...",SID_LOG_OPEN,name_store);

  // Store a name for the cosmology
  ADaPS_store(cosmo,
              name_store,
              "name",
              ADaPS_COPY,
              32);

  // Store the parameters
  SID_log("Omega_Lambda=%le",SID_LOG_COMMENT,Omega_Lambda);
  ADaPS_store(cosmo,
             &Omega_Lambda,
             "Omega_Lambda",
             ADaPS_SCALAR_DOUBLE);
  SID_log("Omega_M     =%le",SID_LOG_COMMENT,Omega_M);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&Omega_M),
             "Omega_M",
             ADaPS_SCALAR_DOUBLE);
  SID_log("Omega_k     =%le",SID_LOG_COMMENT,Omega_k);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&Omega_k),
             "Omega_k",
             ADaPS_SCALAR_DOUBLE);
  SID_log("Omega_b     =%le",SID_LOG_COMMENT,Omega_b);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&Omega_b),
             "Omega_b",
             ADaPS_SCALAR_DOUBLE);
  SID_log("f_gas       =%le",SID_LOG_COMMENT,f_gas);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&f_gas),
             "f_gas",
             ADaPS_SCALAR_DOUBLE);
  SID_log("h_Hubble    =%le",SID_LOG_COMMENT,h_Hubble);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&h_Hubble),
             "h_Hubble",
             ADaPS_SCALAR_DOUBLE);
  SID_log("sigma_8     =%le",SID_LOG_COMMENT,sigma_8);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&sigma_8),
             "sigma_8",
             ADaPS_SCALAR_DOUBLE);
  SID_log("n_spectral  =%le",SID_LOG_COMMENT,n_spectral);
  ADaPS_store((ADaPS **)cosmo,
             (void *)(&n_spectral),
             "n_spectral",
             ADaPS_SCALAR_DOUBLE);
  SID_log("Done.",SID_LOG_CLOSE);
}

