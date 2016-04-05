#ifndef GBPHALOS_AWAKE
#define GBPHALOS_AWAKE

#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpCosmo.h>

//#define MATCH_SCORE_RANK_INDEX -1.5
#define MATCH_SCORE_RANK_INDEX -1.
//#define MATCH_SCORE_RANK_INDEX -TWO_THIRDS

#define MIN_PROFILE_BINS     5
#define MAX_PROFILE_BINS     1000
#define MAX_PROFILE_BINS_P1  1001 // This must be set to MAX_PROFILE_BINS+1

#define READ_GROUPS_DEFAULT         TTTP00
#define READ_GROUPS_ALL             TTTP01
#define READ_GROUPS_SUBGROUPS       TTTP02
#define READ_GROUPS_NOSUBGROUPS     TTTP03
#define READ_GROUPS_PROPERTIES      TTTP04
#define READ_GROUPS_NOPROPERTIES    TTTP05
#define READ_GROUPS_NOIDS           TTTP06
#define READ_GROUPS_MBP_IDS_ONLY    TTTP07
#define READ_GROUPS_PEANOHILBERT    TTTP08

#define READ_CATALOG_GROUPS      TTTP00
#define READ_CATALOG_SUBGROUPS   TTTP01
#define READ_CATALOG_PROPERTIES  TTTP02
#define READ_CATALOG_PROFILES    TTTP03
#define READ_CATALOG_DEFAULT     READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES

#define MATCH_SUBGROUPS     TTTP01 // Match subgroups (default)
#define MATCH_GROUPS        TTTP02 // Match groups
#define MATCH_BACK          TTTP03 // Switch the sence of matching between plists
#define MATCH_STORE_2       TTTP04 // Switch which plist results are stored in
#define MATCH_STORE_SCORE   TTTP05 // Store the matching score
#define MATCH_SUBSTRUCTURE  TTTP06 // Search for substructure; ie. ignore self matches
#define MATCH_APPLY_OFFSETS TTTP07 // When matching across multiple cores, apply
                                   //   rank offsets so that match_ids are global indices

typedef struct halo_properties_info halo_properties_info;
struct halo_properties_info{
  long long id_MBP;                    // ID of most bound particle in structure
  double    M_vir;                     // Bryan & Norman (ApJ 495, 80, 1998) virial mass [M_sol/h]
  int       n_particles;               // Number of particles in the structure
  float     position_COM[3];           // Centre-of-mass position      [Mpc/h]
  float     position_MBP[3];           // Most bound particle position [Mpc/h]
  float     velocity_COM[3];           // Centre-of-mass velocity      [km/s]
  float     velocity_MBP[3];           // Most bound particle velocity [km/s]
  float     R_vir;                     // Virial radius [Mpc/h]
  float     R_halo;                    // Distance of last halo particle from MBP [Mpc/h]
  float     R_max;                     // Radius of maximum circular velocity     [Mpc/h]
  float     V_max;                     // Maximum circular velocity               [km/s]
  float     sigma_v;                   // Total 3D velocity dispersion            [km/s]
  float     spin[3];                   // Specific angular momentum vector        [Mpc/h*km/s]
  float     q_triaxial;                // Triaxial shape parameter q=b/a
  float     s_triaxial;                // Triaxial shape parameter s=c/a
  float     shape_eigen_vectors[3][3]; // Normalized triaxial shape eigenvectors
  char      padding[8];                // Alignment padding
};

typedef struct halo_profile_bin_info halo_profile_bin_info;
struct halo_profile_bin_info{
  float   r_med;                     // The median  radius of particles in this bin [Mpc/h]
  float   r_max;                     // The maximum radius of particles in this bin [Mpc/h] 
  int     n_particles;               // Number of particles in this bin
  float   rho;                       // Density                        [M_sol*h^2/Mpc^3]
  double  M_r;                       // Mass of particles within r_max [M_sol/h]
  float   overdensity;               // Bryan & Norman overdensity within r_max
  float   position_COM[3];           // Position of the centre-of-mass of particles within r_max [Mpc/h]
  float   velocity_COM[3];           // Velocity of the centre-of-mass of particles within r_max [km/s]
  float   sigma_rad;                 // Radial velocity dispersion of particles in this bin      [km/s]
  float   sigma_tan;                 // Tangential velocity dispersion of particles in this bin  [km/s]
  float   sigma_tot;                 // Total 3D velocity dispersion of particles in this bin    [km/s]
  float   spin[3];                   // Specific angular momentum of particles within r_max      [Mpc/h*km/s]
  float   q_triaxial;                // Triaxial shape parameter q=b/a of particles within r_max
  float   s_triaxial;                // Triaxial shape parameter s=c/a of particles within r_max
  float   shape_eigen_vectors[3][3]; // Normalized triaxial shape eigenvectors of particles within r_max
};

typedef struct halo_profile_info halo_profile_info;
struct halo_profile_info{
  int                   n_bins;
  halo_profile_bin_info bins[MAX_PROFILE_BINS];
};

// This datastructure describes the halo catalog file-pointer
typedef struct fp_catalog_info fp_catalog_info;
struct fp_catalog_info{
   char  filename_properties_root[MAX_FILENAME_LENGTH];
   char  filename_properties_base[MAX_FILENAME_LENGTH];
   char  filename_profiles_root[MAX_FILENAME_LENGTH];
   char  filename_profiles_base[MAX_FILENAME_LENGTH];
   FILE *fp_properties;
   FILE *fp_profiles;
   int   snap_num;
   int   i_file;
   int   n_files;
   int   n_halos_total;
   int   i_halo;
   int   i_halo_start;
   int   i_halo_stop;
   int   n_halos_file;
   int   flag_read_properties;
   int   flag_read_profiles;
   int   flag_multifile;
};

// This is the format used as the SAGE structure
typedef struct halo_properties_SAGE_info halo_properties_SAGE_info;
struct halo_properties_SAGE_info{

  // merger tree pointers and match type
  int descendant;
  int progenitor_first;
  int progenitor_next;
  int group_halo_first;
  int group_halo_next;

  // properties of halo
  int       n_particles;
  float     M_Mean200; // Mean 200 values (Mvir=M_Crit200)
  float     M_vir;     // This is the FoF mass for the most massive substructure
  float     M_TopHat;
  float     pos[3];
  float     vel[3];
  float     sigma_v;
  float     v_max;
  float     spin[3];
  long long most_bound_id;

  // original position in halo-finder output
  int   snap_num;
  int   file_num;
  int   halo_index;
  float half_mass;
};

// This is the format used as the SHORT structure
typedef struct halo_properties_SHORT_info halo_properties_SHORT_info;
struct halo_properties_SHORT_info{
  float     pos[3];
  float     vel[3];
  double    M_vir;
  float     sigma_v;
  float     V_max; 
  float     R_vir; 
  int       n_particles;
};

typedef struct halo_MBP_info halo_MBP_info;
struct halo_MBP_info{
  // properties of halo
  GBPREAL   pos[3];
  GBPREAL   vel[3];
  long long most_bound_id;
  int       group_halo_first;

  // original position in halo-finder output
  int   snap_num;
  int   halo_index;
};

typedef struct halo_info halo_info;
struct halo_info{
   int                   flag_use_profiles;
   int                   flag_use_SO;
   int                   snapshot;
   halo_properties_info *properties_group;
   halo_properties_info *properties_subgroup;
   halo_profile_info    *profiles_group;
   halo_profile_info    *profiles_subgroup;
   float                 SO_data_group[6];
   int                   n_sub;
   int                   np_sub;
   int                   np_sub_largest;
};

typedef struct halo_trend_info halo_trend_info;
struct halo_trend_info{
   double *z_list;
   double  m_p;
   double  box_size;
   int     n_snaps;
};

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif

void read_groups(char        *filename_groups_root,
                 int          i_file,
                 int          mode,
                 plist_info  *plist,
                 const char  *catalog_name,...);
void read_groups_AHF(char        *filename_groups_root,
                     int          i_file,
                     int          mode,
                     plist_info  *plist,
                     char        *catalog_name);

int fopen_catalog(char            *filename_catalog_root,
                  int              snapshot_number,
                  int              mode,
                  fp_catalog_info *fp_out);

void init_halo_trend           (halo_trend_info *trend_data,trend_info **trend,const char *name);
void init_halo_trend_coordinate(halo_trend_info *trend_data,trend_info  *trend,const char *name);

void init_halo_trend_property_z(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs);
void free_halo_trend_property_z(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs);
int  calc_halo_trend_property_index_z(trend_property_info *property,hist_info *hist,void *halo);
void init_halo_trend_property_logM_FoF(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs);
void free_halo_trend_property_logM_FoF(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs);
int  calc_halo_trend_property_index_logM_FoF(trend_property_info *property,hist_info *hist,void *halo_in);
void init_halo_trend_property_SSFctn(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs);
void free_halo_trend_property_SSFctn(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs);
int  calc_halo_trend_property_index_SSFctn(trend_property_info *property,hist_info *hist,void *halo_in);

int  fopen_nth_catalog_file(fp_catalog_info *fp_in,int n);
int  fread_catalog_file(fp_catalog_info *fp_in,halo_properties_SHORT_info *properties_short_out,halo_properties_SAGE_info *properties_out,halo_properties_info *properties_all_out,halo_profile_info *profiles_out,int halo_index);
int  fread_catalog_raw(fp_catalog_info *fp_in,halo_properties_info *properties_out,halo_profile_info *profiles_out,int halo_index);
void fclose_catalog(fp_catalog_info *fp_in);
                 
int  compute_group_analysis(halo_properties_info *properties,
                            halo_profile_info    *profile,
                            double (*p_i_fctn) (void *,int,int), 
                            double (*v_i_fctn) (void *,int,int), 
                            size_t (*id_i_fctn)(void *,int), 
                            void   *params,
                            double       box_size,
                            double       particle_mass,
                            int          n_particles,
                            double       expansion_factor,
                            double      *x,
                            double      *y,
                            double      *z,
                            double      *vx,
                            double      *vy,
                            double      *vz,
                            double      *R,
                            size_t     **R_index,
                            int          flag_manual_centre,
                            int          flag_compute_shapes,
                            cosmo_info  *cosmo);
void write_group_properties(FILE *fp,halo_properties_info *properties,    double h_Hubble);
void write_group_profiles(FILE   *fp,halo_profile_info    *profiles,      double h_Hubble,double redshift,cosmo_info *cosmo);
void write_group_analysis(FILE                 *fp_properties,
                          FILE                 *fp_profiles,
                          FILE                 *fp_indices,
                          halo_properties_info *properties,
                          halo_profile_info    *profile,
                          size_t               *R_index,
                          int                   n_particles);
void match_halos(plist_info  *plist_1_in,
                 const char  *catalog_1,
                 int          i_file_1_in,
                 int         *mark_list_1,
                 int          n_mark_1,
                 plist_info  *plist_2_in,
                 const char  *catalog_2,
                 int          i_file_2_in,
                 int         *mark_list_2,
                 int          n_mark_2,
                 const char  *catalog_1to2,
                 int          mode,
                 float        match_weight_rank_index);

#ifdef __cplusplus
}
#endif

#endif

