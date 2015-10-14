#ifndef GBPCOSMO_NFW_ETC_AWAKE
#define GBPCOSMO_NFW_ETC_AWAKE

#include <gbpCosmo_core.h>
#include <gbpCosmo_linear_theory.h>

#define NFW_MODE_DEFAULT 0

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif
void set_NFW_params(double       M,
                    double       z,
                    int          mode,
                    cosmo_info **cosmo,
                    double      *c_vir,
                    double      *R_vir);
double R_vir_NFW(double       M_vir,
                 double       z,
                 int          mode,
                 cosmo_info **cosmo);
double c_vir_NFW(double       M_vir,
                 double       z,
                 int          mode,
                 cosmo_info **cosmo);
double rho_NFW(double       r,
               double       M_vir,
               double       z,
               int          mode,
               cosmo_info **cosmo);
double rho_NFW_fft(double       k,
                   double       M_vir,
                   double       z,
                   int          mode,
                   cosmo_info **cosmo);
double M_r_NFW(double       r,
               double       M_vir,
               double       z,
               int          mode,
               cosmo_info **cosmo);
double V_circ_NFW(double       r,
                  double       M_vir,
                  double       z,
                  int          mode,
                  cosmo_info **cosmo);
double V_circ_vir_NFW(double       M_vir,
                      double       z,
                      int          mode,
                      cosmo_info **cosmo);
double V_max_NFW(double       M_vir,
                 double       z,
                 int          mode,
                 cosmo_info **cosmo);
void init_Vmax_to_Mvir_NFW(cosmo_info **cosmo,
                           int          mode,
                           double       z);
double Vmax_to_Mvir_NFW(double       V_max,
                        double       z,
                        int          mode,
                        cosmo_info **cosmo);
double R_half_V_max_NFW(double       M_vir,
                        double       z,
                        int          mode,
                        cosmo_info **cosmo);
double R_V_max_NFW(double       M_vir,
                   double       z,
                   int          mode,
                   cosmo_info **cosmo);
double Delta_half_V_max_NFW(double       M_vir,
                            double       z,
                            int          mode,
                            cosmo_info **cosmo);
double Delta_half_V_max_NFW(double       M_vir,
                            double       z,
                            int          mode,
                            cosmo_info **cosmo);
double V2_circ(double M, 
               double r);
#ifdef __cplusplus
}
#endif
#endif
