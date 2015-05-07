/**********************/
/* Cosmology routines */
/**********************/
#ifndef GBPCOSMO_MASS_FUNCTIONS_AWAKE
#define GBPCOSMO_MASS_FUNCTIONS_AWAKE

#include <gbpCosmo_core.h>
#include <gbpCosmo_linear_theory.h>

#define MF_JENKINS TTTP01
#define MF_PS      TTTP02
#define MF_ST      TTTP03
#define MF_WATSON  TTTP04

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif
double mass_function(double        M_interp,
                     double        z,
                     cosmo_info  **cosmo,
                     int           select_flag);
double mass_function_cumulative(double       M_interp,
                                double       z,
                                cosmo_info **cosmo,
                                int          select_flag);
double scaled_mass_function(double sigma,int select_flag);
#ifdef __cplusplus
}
#endif
#endif
