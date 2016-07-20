#ifndef SPH_AWAKE
#define SPH_AWAKE
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>

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
#define SMOOTH_DEFAULT           0
#define SMOOTH_DENSITY      TTTP01
#define SMOOTH_SIGMA_V      TTTP02
#define SMOOTH_NOT_PERIODIC TTTP03
#define SMOOTH_LONGIDS      TTTP04

#define READ_SMOOTH_LOG_SIGMA TTTP01
#define READ_SMOOTH_LOG_RHO   TTTP02

#define READ_GADGET_DEFAULT        0
#define READ_GADGET_NO_HUBBLE TTTP00
#define READ_GADGET_NONE      TTTP01 
#define READ_GADGET_RANDOM    TTTP02
#define READ_GADGET_MARKED    TTTP03
#define READ_GADGET_X_MIN     TTTP04
#define READ_GADGET_X_MAX     TTTP05
#define READ_GADGET_Y_MIN     TTTP06
#define READ_GADGET_Y_MAX     TTTP07
#define READ_GADGET_Z_MIN     TTTP08
#define READ_GADGET_Z_MAX     TTTP09

// Stuff for marking particles 
#define MARK_DEFAULT         0
#define MARK_READ_ALL   TTTP00
#define MARK_LIST_ONLY  TTTP01
#define VOLUME_BOX      TTTP02
#define VOLUME_SPHERE   TTTP03
#define MARK_PROPERTY   TTTP04
#define MARK_INIT       TTTP05
#define MARK_AND        TTTP06
#define MARK_OR         TTTP07

#define PROCESS_GADGET_BINARY_DEFAULT         0
#define PROCESS_GADGET_BINARY_ALL_TO_ALL TTTP00

typedef struct gadget_header_info gadget_header_info;
struct gadget_header_info{
  int          n_file[N_GADGET_TYPE];
  double       mass_array[N_GADGET_TYPE];
  double       time;
  double       redshift;
  int          flag_SFr;
  int          flag_feedback;
  unsigned int n_all_lo_word[N_GADGET_TYPE];
  int          flag_cooling;
  int          n_files;
  double       box_size;
  double       Omega_M;
  double       Omega_Lambda;
  double       h_Hubble;
  int          flag_ages;
  int          flag_metals;
  unsigned int n_all_hi_word[N_GADGET_TYPE];
  int          flag_entropyICs;
  char         unused[60];
};

typedef struct gadget_read_info gadget_read_info;
struct gadget_read_info{
   gadget_header_info header;
   size_t             n_all[N_GADGET_TYPE];
   size_t             n_particles;
   char filename_root[MAX_FILENAME_LENGTH];
   int  snapshot_number;
   int  flag_multifile;
   int  flag_file_type;
   int  first_select_call;
   int  first_action_call;
};

typedef struct smooth_header_info smooth_header_info;
struct smooth_header_info{
   int       n_particles_file;
   int       offset;
   long long n_particles_total;
   int       n_files;
};

// Structure to store Abstract Item Lists (AbILs) 
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

// Structure to store particle info 
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

#define PROCESS_GADGET_FILE_BUFFER_SIZE (1024*1024)
#define WRITE_GADGET_BINARY_DEFAULT 0
typedef struct select_gadget_volume_params_info select_gadget_volume_params_info;
struct select_gadget_volume_params_info{
   plist_info *plist;
   GBPREAL     cen[3];
   GBPREAL     size;
   GBPREAL     size2;
   GBPREAL     box_size;
};
typedef struct select_gadget_ids_params_info select_gadget_ids_params_info;
struct select_gadget_ids_params_info{
   plist_info *plist;
   int         n_ids;
   size_t     *id_list;
};

// Function Definitions
#ifdef __cplusplus
extern "C" {
#endif
void process_gadget_file(const char *status_message,
                         char   *filename_read_root,
                         int     snapshot_number,
                         int     select_function(gadget_read_info *fp_gadget,
                                                 void             *params,
                                                 size_t            i_particle,
                                                 size_t            i_particle_type,
                                                 int               i_type,
                                                 GBPREAL          *pos,
                                                 GBPREAL          *vel,
                                                 size_t            ID_i),
                         void    action_function(gadget_read_info *fp_gadget,
                                                 void             *params,
                                                 size_t            i_particle,
                                                 size_t            i_particle_type,
                                                 int               i_type,
                                                 GBPREAL          *pos,
                                                 GBPREAL          *vel,
                                                 size_t            ID_i),
                         void   *params,
                         size_t *n_particles_type_local_pass,
                         size_t *n_particles_type_pass,
                         int    *flag_long_IDs,
                         int     mode);
void allocate_gadget_particles(plist_info *plist,
                               size_t     *n_particles_type_local,
                               size_t     *n_particles_type,
                               int         flag_long_IDs);
int select_gadget_cube(gadget_read_info *fp_gadget,
                       void             *params,
                       size_t            i_particle,
                       size_t            i_particle_type,
                       int               i_type,
                       GBPREAL          *pos,
                       GBPREAL          *vel,
                       size_t            ID_i);
int select_gadget_sphere(gadget_read_info *fp_gadget,
                         void             *params,
                         size_t            i_particle,
                         size_t            i_particle_type,
                         int               i_type,
                         GBPREAL          *pos,
                         GBPREAL          *vel,
                         size_t            ID_i);
void process_gadget_file_fctn_null(gadget_read_info *fp_gadget,
                                   void             *params,
                                   size_t            i_particle,
                                   size_t            i_particle_type,
                                   int               i_type,
                                   GBPREAL          *pos,
                                   GBPREAL          *vel,
                                   size_t            ID_i);
void store_gadget_particles(gadget_read_info *fp_gadget,
                            void             *params,
                            size_t            i_particle,
                            size_t            i_particle_type,
                            int               i_type,
                            GBPREAL          *pos,
                            GBPREAL          *vel,
                            size_t            ID_i);
int select_gadget_ids(gadget_read_info *fp_gadget,
                      void             *params,
                      size_t            i_particle,
                      size_t            i_particle_type,
                      int               i_type,
                      GBPREAL          *pos,
                      GBPREAL          *vel,
                      size_t            ID_i);
int select_gadget_all(gadget_read_info *fp_gadget,
                      void             *params,
                      size_t            i_particle,
                      size_t            i_particle_type,
                      int               i_type,
                      GBPREAL          *pos,
                      GBPREAL          *vel,
                      size_t            ID_i);
void open_gadget_file(char      *filename_root_in,
                      int        snapshot_number,
                      fp_gadget *fp);

int  init_gadget_read(char *filename_root_in,int snapshot_number,gadget_read_info *fp);
int  init_smooth_read(char *filename_root_in,int snapshot_number,int *flag_multifile,int *flag_file_type,smooth_header_info *header);
void set_gadget_filename(gadget_read_info *fp,int i_file,char *filename);
void set_smooth_filename(char *filename_root_in,int snapshot_number,int multifile_number,int flag_multifile,int flag_file_type,char *filename);
void change_gadget_filename(const char *filename_root_in,const char *filename_root,int snapshot_number,int multifile_number,int flag_multifile,int flag_file_type,char *filename);
void allocate_gadget_particles(plist_info *plist,
                               size_t     *n_particles_type_local,
                               size_t     *n_particles_type,
                               int         flag_long_IDs);
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
                    const char *mark_name,
                    const char *filename_in,
                    int         mode);
void read_tipsy_mark_file(char        *filename,
                          plist_info  *plist,
                          char        *mark_name);
void read_tipsy_binary(char       *filename,
                       plist_info *plist);
void read_tipsy_ascii(char       *filename,
                      plist_info *plist);
void write_gadget_binary_new(plist_info  *plist,
                             char        *filename_out_root,
                             int          n_files,
                             int          mode);
void write_tipsy_binary(char       *filename,
                        plist_info *plist);
void write_gadget_ascii(char       *filename,
                        plist_info *plist,
                        size_t     *id_ordering);
void write_gadget_csv(char       *filename_out,
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

