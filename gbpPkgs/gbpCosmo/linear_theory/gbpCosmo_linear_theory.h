/**********************/
/* Cosmology routines */
/**********************/
#ifndef GBPCOSMO_LINEAR_THEORY_AWAKE
#define GBPCOSMO_LINEAR_THEORY_AWAKE

#include <gbpCosmo_core.h>

#define PSPEC_LINEAR_TF   0
#define PSPEC_LINEAR_BBKS 1

#define PSPEC_ALL_MATTER  0
#define PSPEC_DARK_MATTER 1
#define PSPEC_BARYON      2

typedef struct sigma2_integrand_params_struct sigma2_integrand_params;
struct sigma2_integrand_params_struct {
  double      R;
  cosmo_info *cosmo;
  double      z;
  int         mode;
  int         component;
};

typedef struct gbpCosmo2gbpCosmo_info gbpCosmo2gbpCosmo_info;
struct gbpCosmo2gbpCosmo_info {
   double      s_L;
   double      s_M;
   double      M_min;
   double      M_max;
   double      z_min;
   double      z_min_scaled;
   double      D_prime_z_min;
   double      D_z_scaled;
   double      D_ratio;
   cosmo_info *cosmo_source;
   cosmo_info *cosmo_target;
};

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif
double R_of_k(double k);
double k_of_R(double R);
double M_of_R(double R,cosmo_info *cosmo);
double M_of_k(double k,cosmo_info *cosmo);
double R_of_M(double M,cosmo_info *cosmo);
double k_of_M(double M,cosmo_info *cosmo);
double R_NL_z(double z,cosmo_info **cosmo);
void init_gbpCosmo2gbpCosmo(cosmo_info      **cosmo_source,
                            cosmo_info      **cosmo_target,
                            double            z_min,
                            double            M_min,
                            double            M_max,
                            gbpCosmo2gbpCosmo_info *gbpCosmo2gbpCosmo);
double z_gbpCosmo2gbpCosmo(double z,gbpCosmo2gbpCosmo_info *gbpCosmo2gbpCosmo);
double L_gbpCosmo2gbpCosmo(double L,gbpCosmo2gbpCosmo_info *gbpCosmo2gbpCosmo);
double M_gbpCosmo2gbpCosmo(double M,gbpCosmo2gbpCosmo_info *gbpCosmo2gbpCosmo);
void   pspec_names(int   mode,
                   int   component,
                   char *mode_name,
                   char *component_name);
void   init_power_spectrum_TF(cosmo_info **cosmo);
void   init_transfer_function(cosmo_info **cosmo);
double power_spectrum(double k, double z, cosmo_info **cosmo, int mode,int component);
void   init_power_spectrum_variance(cosmo_info **cosmo,int mode,int component);
double power_spectrum_normalization(cosmo_info *cosmo,
                                    int         mode,
                                    int         component);
double power_spectrum_variance(double       k_interp,
                               double       redshift,
                               cosmo_info **cosmo,
                               int          mode,
                               int          component);
double dln_Inv_sigma_dlogM(cosmo_info **cosmo,
                           double       M_interp,
                           int          mode,
                           int          component);
void init_sigma_M(cosmo_info **cosmo,
                  int          mode,
                  int          component);
double sigma_M(cosmo_info **cosmo,
               double       M_interp,
               double       z,
               int          mode,
               int          component);
double sigma_R(cosmo_info **cosmo,
               double       R,
               double       z,
               int          mode,
               int          component);
double ln_sigma_M(cosmo_info *cosmo,
                  double      M_interp,
                  double      z,
                  int         mode,
                  int         component);
double ln_Inv_sigma_M(cosmo_info **cosmo,
                  double      M_interp,
                  double      z,
                  int         mode,
                  int         component);
double dln_sigma_dlnM(cosmo_info *cosmo,
                      double      M_interp,
                      double      z,
                      int         mode,
                      int         component);
double dlogDplus_dloga(double      a,
                       cosmo_info *cosmo);
double dDplus_da(double      a,
                 cosmo_info *cosmo);
double Dplus(double      a,
             cosmo_info *cosmo);
double linear_growth_factor(double       redshift,
                            cosmo_info  *cosmo);
double sigma2_integrand(double  k,
                        void   *params_in);
double W_k_tophat(double kR);
double M_sc(double       z,
            cosmo_info **cosmo,
            int          mode,
            int          component);
double lk_sc(double       z,
             cosmo_info **cosmo,
             int          mode,
             int          component);
#ifdef __cplusplus
}
#endif
#endif
