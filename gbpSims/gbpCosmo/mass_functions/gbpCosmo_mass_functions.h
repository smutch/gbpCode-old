/**********************/
/* Cosmology routines */
/**********************/
#ifndef GBPCOSMO_MASS_FUNCTIONS_AWAKE
#define GBPCOSMO_MASS_FUNCTIONS_AWAKE

#include <gbpCosmo_core.h>
#include <gbpCosmo_linear_theory.h>

#define MF_JENKINS        0
#define MF_PS             1
#define MF_ST             2

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif
double mass_function(double        M_interp,
                     double        z,
                     cosmo_info  **cosmo,
                     int           mode,
                     int           component,
                     int           select_flag);
double mass_function_cumulative(double       M_interp,
                                double       z,
                                cosmo_info  *cosmo,
                                int          mode,
                                int          component,
                                int          select_flag);
double scaled_mass_function(double sigma,int select_flag);
#ifdef __cplusplus
}
#endif
#endif
