#ifndef GBPCOSMO_CORE_AWAKE
#define GBPCOSMO_CORE_AWAKE

// Define the default cosmology if it wasn't specified
//   when make was called.
#ifndef GBP_COSMOLOGY_DEFAULT
#define GBP_COSMOLOGY_DEFAULT "WMAP-7"
#endif

#define DELTAT_A_MIN_A 1e-5

typedef ADaPS cosmo_info;

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif
double a_of_z(double z);
double z_of_a(double a);
void   init_deltat_a(cosmo_info **cosmo);
double deltat_a(cosmo_info **cosmo,double a_1,double a_2);
double R_Delta_z(double      M_Delta,
                 double      Delta,
                 double      redshift,
                 cosmo_info *cosmo);

double Omega_z(double     redshift,cosmo_info *cosmo);
double Delta_vir(double redshift,cosmo_info *cosmo);
double Ha_Ho(double a,cosmo_info *cosmo);
double E_z(double Omega_M, double Omega_k, double Omega_Lambda, double z);
double H_z(double redshift,cosmo_info *cosmo);
double H_convert(double Hz);
double dlna_dtau(double      a,
       cosmo_info *cosmo);
double rho_crit_z(double redshift, cosmo_info *cosmo);
double rho_crit_z_strip(double redshift,
                        double h_Hubble,
                        double Omega_M, 
                        double Omega_Lambda);
double D_Hubble(double h_Hubble);
double E_z(double Omega_M, double Omega_k, double Omega_Lambda, double z);
double D_comove(double z,cosmo_info *cosmo);
double D_comove_transverse(double z,cosmo_info *cosmo);
double D_angular(double z,cosmo_info *cosmo);
double D_angular_1to2(double z_1,double z_2,cosmo_info *cosmo);
double D_luminosity(double z,cosmo_info *cosmo);
void   read_gbpCosmo_file(cosmo_info **cosmo,const char *filename_in);
void   init_cosmo_default(cosmo_info **cosmo);
void   init_cosmo(cosmo_info **cosmo,
                  const char  *name,
                  double       Omega_Lambda,
                  double       Omega_M,
                  double       Omega_k,
                  double       Omega_b,
                  double       f_gas,
                  double       h_Hubble,
                  double       sigma_8,
                  double       n_spectral);
void free_cosmo(cosmo_info **cosmo);

#ifdef __cplusplus
}
#endif
#endif
