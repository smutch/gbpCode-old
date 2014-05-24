#ifndef SPH_AWAKE
#define SPH_AWAKE
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>

#define MAKE_MAP_LOG           2
#define MAKE_MAP_APPLY_DIMMING 4
#define MAKE_MAP_COLOUR        8
#define MAKE_MAP_LUMINOSITY    16

#define GADGET_BUFFER_SIZE  SIZE_OF_MEGABYTE

#define GADGET_HEADER_SIZE  256
#define N_GADGET_TYPE       6
#define GADGET_TYPE_GAS     0
#define GADGET_TYPE_DARK    1
#define GADGET_TYPE_STAR    2
#define GADGET_TYPE_STAR2   3
#define GADGET_TYPE_STAR3   4
#define GADGET_TYPE_STAR4   5

#define GADGET_LENGTH       3.085678e22
#define GADGET_VELOCITY     1.e3
#define GADGET_MASS         1.989e40

#define SMOOTH_N_QUANTITIES 3
#define SMOOTH_DEFAULT      0
#define SMOOTH_DENSITY      2
#define SMOOTH_SIGMA_V      4
#define SMOOTH_NOT_PERIODIC 8
#define SMOOTH_LONGIDS      16

#define READ_SMOOTH_LOG_SIGMA 2
#define READ_SMOOTH_LOG_RHO   4

#define N_KERNEL_TABLE      20000
#define SPH_KERNEL_2D       2
#define SPH_KERNEL_GADGET   4
#define SPH_KERNEL_GASOLINE 8

#define READ_GADGET_DEFAULT      0
#define READ_GADGET_NO_HUBBLE    1
#define READ_GADGET_NONE         2
#define READ_GADGET_RANDOM       8
#define READ_GADGET_MARKED      16
#define READ_GADGET_X_MIN       32
#define READ_GADGET_X_MAX       64
#define READ_GADGET_Y_MIN      128
#define READ_GADGET_Y_MAX      256
#define READ_GADGET_Z_MIN      512
#define READ_GADGET_Z_MAX     1024

/*******************************/
/* Stuff for marking particles */
/*******************************/
#define MARK_DEFAULT    0
#define MARK_READ_ALL   1
#define MARK_LIST_ONLY  2
#define VOLUME_BOX      4
#define VOLUME_SPHERE   8
#define MARK_PROPERTY  16
#define MARK_INIT      32
#define MARK_AND       64
#define MARK_OR       128

typedef struct gadget_header_info gadget_header_info;
struct gadget_header_info{
  int          n_file[N_GADGET_TYPE];
  double       mass_array[N_GADGET_TYPE];
  double       time;
  double       redshift;
  int          flag_SFr;
  int          flag_feedback;
  unsigned int n_all[N_GADGET_TYPE];
  int          flag_cooling;
  int          n_files;
  double       box_size;
  double       Omega_M;
  double       Omega_Lambda;
  double       h_Hubble;
  int          flag_ages;
  int          flag_metals;
  unsigned int n_all_high_word[N_GADGET_TYPE];
  int          flag_entropyICs;
  char         unused[60];
};

typedef struct smooth_header_info smooth_header_info;
struct smooth_header_info{
   int       n_particles_file;
   int       offset;
   long long n_particles_total;
   int       n_files;
};

/**************************************************/
/* Structure to store Abstract Item Lists (AbILs) */
/**************************************************/
typedef struct AbIL_info AbIL_info;
struct AbIL_info{
  int          n_species;
  char       **species;
  big_int     *n_local;
  big_int     *n_global;
  big_int     *n_allocated;
  int         *flag_used;
  ADaPS      **data;
  ADaPS       *species_data;
  ADaPS       *constants;
  double       time;
  double       mass_unit;
  double       length_unit;
  double       time_unit;
  double       velocity_unit;
};

/************************************/
/* Structure to store particle info */
/************************************/
typedef struct plist_info plist_info;
struct plist_info{
  char      **species;
  int         n_species;
  double      time;
  double      mass_unit;
  double      length_unit;
  double      time_unit;
  double      velocity_unit;
  size_t      n_gas;
  size_t      n_dark;
  size_t      n_star;
  size_t      n_particle;
  double      d_x;
  double      d_y;
  double      d_z;
  double      d_vx;
  double      d_vy;
  double      d_vz;
  double      d_alpha;
  double      d_beta;
  double      d_gamma;
  ADaPS      *data;
};

typedef struct markfile_header_info markfile_header_info;
struct markfile_header_info {
  int n_type;
  int n_mark_species[N_GADGET_TYPE];
};

typedef struct smoothfile_header_info smoothfile_header_info;
struct smoothfile_header_info {
  int  n_particles;
  char species[256];
  int  n_quantities;
  int  n_quantities_used;
  int  used[SMOOTH_N_QUANTITIES];
};

typedef struct fp_gadget fp_gadget;
struct fp_gadget{
  FILE *fp;
  int   flag_filefound;
  int   flag_multifile;
  int   flag_file_type;
};

// Function Definitions
#ifdef __cplusplus
extern "C" {
#endif
void open_gadget_file(char      *filename_root_in,
                      int        snapshot_number,
                      fp_gadget *fp);

int  init_gadget_read(char *filename_root_in,int snapshot_number,int *flag_multifile,int *flag_file_type,gadget_header_info *header);
int  init_smooth_read(char *filename_root_in,int snapshot_number,int *flag_multifile,int *flag_file_type,smooth_header_info *header);
void set_gadget_filename(char *filename_root_in,int snapshot_number,int multifile_number,int flag_multifile,int flag_file_type,char *filename);
void set_smooth_filename(char *filename_root_in,int snapshot_number,int multifile_number,int flag_multifile,int flag_file_type,char *filename);
void change_gadget_filename(const char *filename_root_in,const char *filename_root,int snapshot_number,int multifile_number,int flag_multifile,int flag_file_type,char *filename);

void init_plist(plist_info *plist, slab_info *slab,double length_unit,double mass_unit, double velocity_unit); 
void free_plist(plist_info *plist);
void standard_to_system(plist_info *plist, char *species_name);
void close_plist(plist_info *plist);
void translate_sph(plist_info *plist,
                   double      x_cen,
                   double      y_cen,
                   double      z_cen,
                   double      vx_cen,
                   double      vy_cen,
                   double      vz_cen);
void fetch_particles(plist_info *plist,
                     char       *var_text,
                     char       *ptype_in_text,
                     char       *mark_name,
                     void      **rval,
                     size_t     *n_rval);
size_t mark_particles(plist_info *plist,
                      int         run_mode,
                      double     *input_vars,
                      const char *mark_name);
void compute_particle_quantity(plist_info *plist,
                               char       *var_text,
                               char       *ptype_in_text);
void read_mark_file(plist_info *plist,
                    char      *mark_name,
                    char      *filename_in,
                    int        mode);
void read_tipsy_mark_file(char        *filename,
                          plist_info  *plist,
                          char        *mark_name);
void read_tipsy_binary(char       *filename,
                       plist_info *plist);
void read_tipsy_ascii(char       *filename,
                      plist_info *plist);
void write_tipsy_binary(char       *filename,
                        plist_info *plist);
void write_ascii(char       *filename,
                 plist_info *plist);
void write_csv(char       *filename_out,
               plist_info *plist);
void write_mark_file(plist_info *plist,
                     const char *mark_name,
                     const char *filename_out);
void read_tipsy_ids(char       *filename,
                    plist_info *plist);
void read_smooth(plist_info *plist,
                 char       *filename_root,
                 int         snapshot_number,
                 int         mode);
int sph_project(plist_info  *plist,
                double       xmin,
                double       xmax,
                double       ymin,
                double       ymax,
                double       zmin,
                double       zmax,
                int          nx,
                int          ny,
                double     **image,
                int        **mask,
                int          proj_type,
                char        *parameter);
int sph_project2(plist_info  *plist,
                 double       x_c,
                 double       y_c,
                 double       z_c,
                 double       x_t,
                 double       y_t,
                 double       z_t,
                 double       FOV_x,
                 double       FOV_y,
                 double       alpha,
                 int          nx,
                 int          ny,
                 double     **image,
                 int        **mask,
                 int          mode,
                 char        *parameter);
void make_map(plist_info  *plist,
              double       x_c,
              double       y_c,
              double       z_c,
              double       x_t,
              double       y_t,
              double       z_t,
              double       FOV_x,
              double       FOV_y,
              double       box_size,
              int          nx,
              int          ny,
              double     **image,
              int        **mask,
              int          mode,
              char        *parameter);
void write_smooth(plist_info *plist,
                  char       *species_name,
                  char       *filename,
                  int         n_chunks,
                  int         mode);
void prep_types(char ***pname, int n_type, ... );
void read_gadget_binary_header(char        *filename_root_in,
                               int          snapshot_number,
                               plist_info  *plist);
void read_gadget_binary(char       *filename_root_in,
                        int        snapshot_number,
                        plist_info *plist,
                        int         mode);
void write_gadget_binary(char       *filename,
                         plist_info *plist);
void display_gadget_header(plist_info  *plist);
void display_gadget_info(plist_info  *plist);
void free_types(char ***pname,int n_species);

void system_to_standard(plist_info *plist, char *species_name);
void assign_particles_to_grid(plist_info *plist,
                              double     *L,
                              int        *n,
                              char       *species_name);
void smooth_grid(plist_info *plist,
                 int         grid_size,
                 char       *species_name,
                 int         mode);
#ifdef __cplusplus
}
#endif
#endif

