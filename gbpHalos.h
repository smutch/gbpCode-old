#include <gbpSPH.h>
#include <gbpCosmo.h>

#define WRITE_GROUPS_ASCII

#define MIN_PROFILE_BINS     5
#define MAX_PROFILE_BINS     1000
#define MAX_PROFILE_BINS_P1  1001 // This must be set to MAX_PROFILE_BINS+1

#define READ_GROUPS_DEFAULT      1
#define READ_GROUPS_ALL          2
#define READ_GROUPS_SUBGROUPS    4
#define READ_GROUPS_NOSUBGROUPS  8
#define READ_GROUPS_PROPERTIES   16
#define READ_GROUPS_NOPROPERTIES 32
#define READ_GROUPS_NOIDS        64

#define MATCH_SUBGROUPS     2 // Match subgroups (default)
#define MATCH_GROUPS        4 // Match groups
#define MATCH_BACK          8 // Switch the sence of matching between plists
#define MATCH_STORE_2      16 // Switch which plist results are stored in
#define MATCH_STORE_SCORE  32 // Store the matching score
#define MATCH_SUBSTRUCTURE 64 // Search for substructure; ie. ignore self matches

#define TREE_PROGENITOR_ORDER_DEFAULT 0
#define TREE_PROGENITOR_ORDER_DELUCIA 2

typedef struct halo_properties_info halo_properties_info;
struct halo_properties_info{
  long long id_MBP;                    // ID of most bound particle in structure
  int       n_particles;               // Number of particles in the structure
  float     position_COM[3];           // Centre-of-mass position      [Mpc/h]
  float     position_MBP[3];           // Most bound particle position [Mpc/h]
  float     velocity_COM[3];           // Centre-of-mass velocity      [km/s]
  float     velocity_MBP[3];           // Most bound particle velocity [km/s]
  double    M_vir;                     // Bryan & Norman (ApJ 495, 80, 1998) virial mass [M_sol/h]
  float     R_vir;                     // Virial radius [Mpc/h]
  float     R_halo;                    // Distance of last halo particle from MBP [Mpc/h]
  float     R_max;                     // Radius of maximum circular velocity     [Mpc/h]
  float     V_max;                     // Maximum circular velocity               [km/s]
  float     sigma_v;                   // Total 3D velocity dispersion            [km/s]
  float     spin[3];                   // Specific angular momentum vector        [Mpc/h*km/s]
  float     q_triaxial;                // Triaxial shape parameter q=b/a
  float     s_triaxial;                // Triaxial shape parameter s=c/a
  float     shape_eigen_vectors[3][3]; // Normalized triaxial shape eigenvectors
};

typedef struct halo_profile_bin_info halo_profile_bin_info;
struct halo_profile_bin_info{
  float   r_med;                     // The median  radius of particles in this bin [Mpc/h]
  float   r_max;                     // The maximum radius of particles in this bin [Mpc/h] 
  int     n_particles;               // Number of particles in this bin
  double  M_r;                       // Mass of particles within r_max [M_sol/h]
  float   rho;                       // Density                        [M_sol*h^2/Mpc^3]
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

typedef struct halo_info halo_info;
struct halo_info{
  // merger tree pointers
  int descendant;
  int progenitor_first;
  int progenitor_next;
  int group_halo_first;
  int group_halo_next;

  // properties of halo
  int       n_particles;
  float     M_vir;
  float     R_vir;
  float     pos[3];
  float     vel[3];
  float     sigma_v;
  float     v_max;
  float     spin[3];
  long long most_bound_id;

  // original position in halo-finder output
  int   snap_num;
  int   halo_index;
  int   halo_id;
  int   group_id;
};

typedef struct halo_MBP_info halo_MBP_info;
struct halo_MBP_info{
  // properties of halo
  REAL      pos[3];
  REAL      vel[3];
  long long most_bound_id;
  int       group_halo_first;

  // original position in halo-finder output
  int   snap_num;
  int   halo_index;
};

typedef struct tree_node_info tree_node_info;
struct tree_node_info{
  halo_info       halo;
  int             depth_first_index;
  int             group_id;
  int             halo_id;
  int             n_progenitors;
  tree_node_info *descendant;
  tree_node_info *progenitor_first;
  tree_node_info *progenitor_next;
  tree_node_info *progenitor_last;
  tree_node_info *group_halo_first;
  tree_node_info *group_halo_next;
  tree_node_info *neighbour_halo_next;
  tree_node_info *next; // Points to the next halo, in the order they are added to the tree
};

typedef struct tree_info tree_info;
struct tree_info{
  tree_node_info  *root;
  tree_node_info  *last_leaf;
  int             *n_neighbours;
  tree_node_info **neighbour_halos;
  tree_node_info **neighbour_halo_last;
};

void read_groups(char        *filename_groups_root,
                 int          i_file,
                 int          mode,
                 plist_info  *plist,
                 char        *catalog_name);
void read_groups_AHF(char        *filename_groups_root,
                     int          i_file,
                     int          mode,
                     plist_info  *plist,
                     char        *catalog_name);
                 
void read_trees(char             *filename_trees_root,
                int               i_file_start,
                int               i_file_stop,
                int               mode,
                tree_node_info  **trees);
int  compute_group_analysis(halo_properties_info *properties,
                            halo_profile_info    *profile,
                            size_t               *id_array,
                            REAL                 *x_array,
                            REAL                 *y_array,
                            REAL                 *z_array,
                            REAL                 *vx_array,
                            REAL                 *vy_array,
                            REAL                 *vz_array,
                            size_t               *ids_sort_index,
                            double                box_size,
                            double                h_Hubble,
                            double                Omega_M,
                            double                particle_mass,
                            int                   n_particles,
                            double                redshift,
                            cosmo_info           *cosmo);
void read_group_properties(FILE  *fp,halo_info            *properties_out,int   halo_index,int i_read);
void write_group_properties(FILE *fp,halo_properties_info *properties,    double h_Hubble);
void write_group_profiles(FILE   *fp,halo_profile_info    *profiles,      double h_Hubble,double redshift,cosmo_info *cosmo);
void write_group_analysis(FILE                 *fp_properties,
                          FILE                 *fp_profiles,
                          halo_properties_info *properties,
                          halo_profile_info    *profile);
void match_halos(plist_info  *plist_1_in,
                 int          i_file_1_in,
                 int         *mark_list_1,
                 int          n_mark_1,
                 plist_info  *plist_2_in,
                 int          i_file_2_in,
                 int         *mark_list_2,
                 int          n_mark_2,
                 char        *catalog_1to2,
                 int          mode);

void compute_trees_horizontal(char *filename_halos_root_in,
                              char *filename_cat_root_in,
                              char *filename_root_out,
                              int   i_read_start,
                              int   i_read_stop,
                              int   i_read_step,
                              int   n_search,
                              int  *flag_clean);
void compute_trees_vertical(char *filename_root_out,
                            int   i_read_start,
                            int   i_read_stop,
                            int   i_read_step,
                            int   n_search,
                            int   n_files_groups,
                            int   n_files_subgroups,
                            int  *flag_clean);
void compute_trees_auxiliary(char *filename_root,
                             char *filename_snapshot_root,
                             char *filename_root_out,
                             int   i_read_start,
                             int   i_read_stop,
                             int   i_read_step,
                             int   n_search,
                             int   n_files_groups,
                             int   n_files_subgroups,
                             int  *flag_clean);

