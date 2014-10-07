#ifndef GBPCOSMO_BIAS_AWAKE
#define GBPCOSMO_BIAS_AWAKE

#include <gbpCosmo_core.h>
#include <gbpCosmo_linear_theory.h>
#include <gbpCosmo_NFW_etc.h>

// Bias model stuff
#define BIAS_MODEL_TRK                  2
#define BIAS_MODEL_BPR                  4
#define BIAS_MODEL_POOLE_SUBSTRUCTURE   8
#define BIAS_MODEL_POOLE_HALO          16
#define BIAS_MODEL_POOLE_ZSPACE        32
#define BIAS_MODEL_POOLE_TOTAL         64
#define BIAS_MODEL_VMAX_ORDINATE      128
#define BIAS_MODEL_KAISER_BOOST       256
#define BIAS_MODEL_KAISER             512

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif
double bias_model(double       x_in,
                  double       delta_c,
                  double       z,
                  cosmo_info **cosmo,
                  int          mode);
#ifdef __cplusplus
}
#endif
#endif
