#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpCosmo.h>

#define PSPEC_ADD_VX     256   // Must start at 256 to allow for MAP2GRID flags
#define PSPEC_ADD_VY     512
#define PSPEC_ADD_VZ     1024
#define PSPEC_DEFAULT    2048

#define CFUNC_ADD_VX     32
#define CFUNC_ADD_VY     64
#define CFUNC_ADD_VZ     128
#define CFUNC_DEFAULT    256

void compute_cor_func(plist_info  *plist,
                      int          mode,
                      cosmo_info  *cosmo,
                      char        *species_name,
                      char        *random_name,
                      double       redshift,
                      double       box_size,
                      double       r_min_l1D,
                      double       r_max_1D,
                      double       r_min_2D,
                      double       r_max_2D,
                      int          n_1D,
                      int          n_2D,
                      int          n_jack,
                      double      *CFUNC_l1D,
                      double      *dCFUNC_l1D,
                      double      *COVMTX_l1D,
                      double      *CFUNC_1D,
                      double      *dCFUNC_1D,
                      double      *COVMTX_1D,
                      double      *CFUNC_2D,
                      double      *dCFUNC_2D,
                      double      *COVMTX_2D,
                      int         *flag_compute_RR,
                      long long  **DD_l1D,
                      long long  **DR_l1D,
                      long long  **RR_l1D,
                      long long  **DD_1D,
                      long long  **DR_1D,
                      long long  **RR_1D,
                      long long  **DD_2D,
                      long long  **DR_2D,
                      long long  **RR_2D);
void compute_power_spectrum(plist_info  *plist,
                            field_info  *FFT,
                            int          distribution_scheme,
                            int          mode,
                            cosmo_info  *cosmo,
                            char        *species_name,
                            double       redshift,
                            int          n_k_1D,
                            double       k_min_1D,
                            double       k_max_1D,
                            int          n_k_2D,
                            double       k_min_2D,
                            double       k_max_2D,
                            double      *k_1D,
                            double      *P_k_1D,
                            double      *dP_k_1D,
                            int         *n_modes_1D,
                            double      *P_k_2D,
                            double      *dP_k_2D,
                            int         *n_modes_2D);
