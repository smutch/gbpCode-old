#ifndef GBPCLUSTERING_AWAKE
#define GBPCLUSTERING_AWAKE
#include <stdarg.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>

#define MAP2GRID_DIST_NGP   1
#define MAP2GRID_DIST_CIC   2
#define MAP2GRID_DIST_TSC   4
#define MAP2GRID_DIST_DWT12 8
#define MAP2GRID_DIST_DWT20 16

#define MAP2GRID_MODE_DEFAULT 0
#define MAP2GRID_MODE_NOCLEAN 2
#define MAP2GRID_MODE_NONORM  4

#define PSPEC_ADD_VX     256   // Must start at 256 to allow for MAP2GRID flags
#define PSPEC_ADD_VY     512
#define PSPEC_ADD_VZ     1024
#define PSPEC_DEFAULT    2048

#define CFUNC_ADD_VX     32
#define CFUNC_ADD_VY     64
#define CFUNC_ADD_VZ     128
#define CFUNC_DEFAULT    256

#define READ_GROUPING_SLAB    1
#define READ_GROUPING_ADD_VX  2
#define READ_GROUPING_ADD_VY  4
#define READ_GROUPING_ADD_VZ  8
#define READ_GROUPING_PHK     16
#define READ_GROUPING_DEFAULT READ_GROUPING_SLAB

// This structure stores everything pertaining to a
//   power spectrum calculation
typedef struct pspec_info pspec_info;
struct pspec_info {
   int         mass_assignment_scheme;
   int         initialized;
   double      redshift; 
   double      box_size; 
   int         grid_size;
   double      k_min_1D; 
   double      k_max_1D; 
   double      dk_1D;    
   double      k_min_2D; 
   double      k_max_2D; 
   double      dk_2D;    
   int         n_k_1D;
   int         n_k_2D;
   double     *k_1D;
   int        *n_modes_1D;
   int        *n_modes_2D;
   double    **P_k_1D;
   double    **dP_k_1D;
   double    **P_k_2D;
   double    **dP_k_2D;
   int         flag_processed[4];
   cosmo_info *cosmo;
   field_info  FFT;
};

// This structure stores everything pertaining to a
//   power spectrum calculation
typedef struct cfunc_info cfunc_info;
struct cfunc_info {
   int          initialized;
   int          flag_compute_RR;
   int          n_1D;
   int          n_2D;
   int          n_jack;
   int          n_2D_total;
   int          n_jack_total;
   int          F_random;
   int          n_data;
   int          n_random;
   double       r_max;
   double       dr_l1D;
   double       dr_1D;
   double       dr_2D;
   double       r_min_l1D;
   double       r_max_1D;
   double       r_min_2D;
   double       r_max_2D;
   double     **CFUNC_l1D;
   double     **dCFUNC_l1D;
   double     **COVMTX_l1D;
   double     **CFUNC_1D;
   double     **dCFUNC_1D;
   double     **COVMTX_1D;
   double     **CFUNC_2D;
   double     **dCFUNC_2D;
   double     **COVMTX_2D;
   long long  **DD_l1D;
   long long  **DR_l1D;
   long long  **RR_l1D;
   long long  **DD_1D;
   long long  **DR_1D;
   long long  **RR_1D;
   long long  **DD_2D;
   long long  **DR_2D;
   long long  **RR_2D;
   double       redshift;
   double       box_size;
   int          n_bits_PHK;
   int          PHK_width;
   cosmo_info  *cosmo;
};

// Function Definitions
#ifdef __cplusplus
extern "C" {
#endif

void read_groupings(char *filename_root,int grouping_number,plist_info *plist,int mode,...);
void read_atable(char *filename_in,plist_info *plist,int x_column,int y_column,int z_column,int vx_column,int vy_column,int vz_column,int mode,...);
void generate_randoms(cfunc_info *cfunc,plist_info *plist,char *species_name,char *filename_out_root,char *random_name);
void map_to_grid(size_t      n_particles_local,
                 GBPREAL    *x_particles_local,
                 GBPREAL    *y_particles_local,
                 GBPREAL    *z_particles_local,
                 GBPREAL    *m_particles_local,
                 cosmo_info *cosmo,
                 double      redshift,
                 int         distribution_scheme,
                 double      normalization_target,
                 field_info *field,
                 int         mode);

// Correlation function stuff
void init_cfunc(cfunc_info *cfunc,int   n_data,  int   F_random,int PHK_width,
                double redshift, double box_size,int n_jack,
                double r_min_l1D,double r_max_1D,double dr_1D,
                double r_min_2D, double r_max_2D,double dr_2D);
void free_cfunc(cfunc_info *cfunc);
void compute_cfunc(plist_info  *plist,
                   char        *species_name,
                   char        *random_name,
                   cfunc_info  *cfunc,
                   int          i_run);
void write_cfunc(cfunc_info *cfunc,char *filename_out_root,plist_info *plist,char *species_name,char *randoms_name,int i_run);

// Power spectrum stuff
void init_pspec(pspec_info *pspec,
                int mass_assignment_scheme,
                double redshift,double box_size,int grid_size,
                double k_min_1D,double k_max_1D,double dk_1D,
                double k_min_2D,double k_max_2D,double dk_2D);
void free_pspec(pspec_info *pspec);
void compute_pspec(plist_info  *plist,
                   char        *species_name,
                   pspec_info  *pspec,
                   int          i_run);
void write_grid(field_info *field,char *filename_out_root,int i_run,int n_run,int mass_assignment_scheme,double box_size);
void write_pspec(pspec_info *pspec,char *filename_out_root,plist_info *plist,char *species_name);

#ifdef __cplusplus
}
#endif
#endif

