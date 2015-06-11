/**********************/
/* Cosmology routines */
/**********************/
#ifndef GBPCOSMO_MASS_FUNCTIONS_AWAKE
#define GBPCOSMO_MASS_FUNCTIONS_AWAKE

#include <gbpCosmo_core.h>
#include <gbpCosmo_linear_theory.h>

#define MF_PASS_PARAMS TTTP01
#define MF_JENKINS     TTTP02
#define MF_PS          TTTP03
#define MF_ST          TTTP04
#define MF_WATSON      TTTP05

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif
double mass_function(double        M_interp,
                     double        z,
                     cosmo_info  **cosmo,
                     int           select_flag,...);
double mass_function_cumulative(double       M_interp,
                                double       z,
                                cosmo_info **cosmo,
                                int          select_flag,...);
double scaled_mass_function(double sigma,int mode,...);
#ifdef __cplusplus
}
#endif
#endif
