#ifndef RENDER_AWAKE
#define RENDER_AWAKE
#include <stdint.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpTrees.h>
#include <gd.h>
#if USE_FFMPEG
  #include <avformat.h>
  #include <swscale.h>
#endif

#define MAKE_MAP_DEFAULT       0
#define MAKE_MAP_LOG           TTTP01
#define MAKE_MAP_APPLY_DIMMING TTTP02
#define MAKE_MAP_COLOUR        TTTP03
#define MAKE_MAP_LUMINOSITY    TTTP04
#define MAKE_MAP_NO_WEIGHTING  TTTP05
#define MAKE_MAP_MODE_RHO      TTTP06
#define MAKE_MAP_MODE_SIGMA    TTTP07
#define MAKE_MAP_INV_SIGMA     TTTP08

#define N_KERNEL_TABLE      20000
#define SPH_KERNEL_2D       TTTP01
#define SPH_KERNEL_GADGET   TTTP02
#define SPH_KERNEL_GASOLINE TTTP03
#define SPH_KERNEL_GAUSSIAN TTTP04

#define CAMERA_MONO            0
#define CAMERA_STEREO          2
#define CAMERA_PLANE_PARALLEL  4
#define CAMERA_DEFAULT         CAMERA_MONO

#define CAMERA_RGB_MODE_1CHANNEL TTTP01
#define CAMERA_RGB_MODE_3CHANNEL TTTP02
#define CAMERA_RGB_MODE_DEFAULT  CAMERA_RGB_MODE_1CHANNEL|CAMERA_RGB_MODE_3CHANNEL

#define RENDER_INIT_PERSPECTIVE 1
#define RENDER_INIT_EVOLVE      2
#define RENDER_INIT_DEFAULT     RENDER_INIT_PERSPECTIVE

#define SET_RENDER_DEFAULT 2
#define SET_RENDER_RESCALE 4
#define SET_RENDER_GADGET  8

#define WRITE_FRAME_DEFAULT  2
#define WRITE_FRAME_PNG_ONLY 4

#define RENDER_SWS_FLAGS      SWS_BICUBIC
#define RENDER_FORMAT_DEFAULT "mpeg"
#define RENDER_AV_LOG_LEVEL   -1
#define RENDER_GOP_SIZE        2
#define RENDER_QSCALE          1

#define WRITE_IMAGE_PNG      TTTP01
#define WRITE_IMAGE_RAW      TTTP02
#define WRITE_IMAGE_DEFAULT  (WRITE_IMAGE_PNG|WRITE_IMAGE_RAW)

#define READ_GADGET_RENDER_SCATTER    2
#define READ_GADGET_RENDER_ID_ORDERED 4
#define READ_GADGET_RENDER_DEFAULT    READ_GADGET_RENDER_SCATTER

#define RENDER_INVALID_SSIMPL_DIR ":%* invalid directory *%:"

// Data structure which holds all info about an image
typedef struct image_info image_info;
struct image_info{
  gdImagePtr       gd_ptr;
  int              width;
  int              height;
  int              n_pixels;
  double           range[2];
  int            **colour_table;
  double          *values;
  int              colourmapselect;
  int              n_colours;
};

// Data structure which holds all info about a movie
typedef struct movie_info movie_info;
struct movie_info{
  #if USE_FFMPEG
    AVStream        *video_stream;
    AVFormatContext *video_context;
    AVFrame         *picture;
    AVFrame         *temp_picture;
  #endif
  uint8_t         *video_outbuf;
  size_t           video_outbuf_size;
  int              width;
  int              height;
  int              n_pixels;
  int              n_frames;
  int              frame_count;
  int              frame_rate; // frames per second
  double           total_duration;
};

typedef struct perspective_info perspective_info;
struct perspective_info{
  double p_o[3];          // Object position
  double p_c[3];          // Camera position; computed from p_o, radius zeta and theta
  double d_o;             // Separation between object and camera
  double radius;          // Separation between object and camera modulo phi
  double theta;           // Rotation angle (azimuthal)   of camera position about x,y,z_o (after c-o transformation)
  double zeta;            // Rotation angle (altitudinal) of camera position about x,y,z_o (after c-o transformation)
  double phi;             // Zoom factor
  double FOV;             // Field of view at object position
  double time;            // Used for time-evolving movies
  double focus_shift_x;   // After all transformations are applied, this shifts everything in the image-frame.  Useful
  double focus_shift_y;   //    for cases where you want flag_comoving=FALSE but don't want (0,0,0) in the image-centre.
  perspective_info *next; // Used for perspective lists
};

// Perspectives are evenly spaced (in therms of frames) accross a scene
typedef struct perspective_interp_info perspective_interp_info;
struct perspective_interp_info{
  interp_info *p_o[3];  
  interp_info *radius;     
  interp_info *FOV;     
  interp_info *theta;   
  interp_info *zeta;   
  interp_info *phi;     
  interp_info *time;
};

typedef struct scene_info scene_info;
struct scene_info{
  int                      n_frames;
  int                      first_frame;
  int                      last_frame;
  int                      n_perspectives;
  perspective_info        *perspectives;
  perspective_info        *first_perspective;
  perspective_info        *last_perspective;
  perspective_info        *evolve;  // Sets magnitudes of linear gradients accross the scene
  perspective_interp_info *interp;
  int                      sealed;  // Set to TRUE if scene initialization has been finalized
  scene_info              *next;
};

typedef struct camera_info camera_info;
struct camera_info{
  int               camera_mode;
  int               colour_table;
  int               flag_velocity_space;
  int               width;
  int               height;
  double            stereo_ratio;
  double            f_near_field;
  double            f_taper_field;
  double            f_image_plane;
  perspective_info *perspective; // Present perspective state of camera
  image_info       *image_RGBY;
  image_info       *image_RGBY_left;
  image_info       *image_RGBY_right;
  char             *mask_RGBY;
  char             *mask_RGBY_left;
  char             *mask_RGBY_right;
  image_info       *image_RGBY_3CHANNEL;
  image_info       *image_RGBY_3CHANNEL_left;
  image_info       *image_RGBY_3CHANNEL_right;
  char             *mask_RGBY_3CHANNEL;
  char             *mask_RGBY_3CHANNEL_left;
  char             *mask_RGBY_3CHANNEL_right;
  int               RGB_mode;
  int               flag_calc_Z_image;
  char              RGB_param[64];
  double            RGB_range[2];
  interp_info      *RGB_gamma;
  ADaPS            *transfer_list;
  image_info       *image_RGB;
  image_info       *image_RGB_left;
  image_info       *image_RGB_right;
  char             *mask_RGB;
  char             *mask_RGB_left;
  char             *mask_RGB_right;
  int               Y_mode;
  char              Y_param[64];
  double            Y_range[2];
  interp_info      *Y_gamma;
  image_info       *image_Y;
  image_info       *image_Y_left;
  image_info       *image_Y_right;
  char             *mask_Y;
  char             *mask_Y_left;
  char             *mask_Y_right;
  int               Z_mode;
  double            Z_range[2];
  interp_info      *Z_gamma;
  image_info       *image_Z;
  image_info       *image_Z_left;
  image_info       *image_Z_right;
  // These images are used when colours can not be set with a colour table
  image_info       *image_RY;
  image_info       *image_RY_left;
  image_info       *image_RY_right;
  image_info       *image_GY;
  image_info       *image_GY_left;
  image_info       *image_GY_right;
  image_info       *image_BY;
  image_info       *image_BY_left;
  image_info       *image_BY_right;
};

typedef struct mark_arg_info mark_arg_info;
struct mark_arg_info{
   char            species[32];
   char            type[32];
   double          dval[8];
   int             ival[8];
   char            value;
   mark_arg_info *next;
};

typedef struct process_halo_info process_halo_info;
struct process_halo_info{
    int    n_particles;
    size_t offset;
    int    id;
    int    tree_case;
    int    descendant_id;
    int    tree_id;
    int    file_offset;
    int    file_index;
    int    n_subgroups;
    void  *ids;
};

typedef struct render_info render_info;
struct render_info{
  double         *kernel_radius;
  double         *kernel_table;
  double         *kernel_table_3d;
  double          kernel_table_avg;
  camera_info    *camera;
  scene_info     *scenes;
  scene_info     *first_scene;
  scene_info     *last_scene;
  double          f_interpolate;
  int             n_interpolate;
  int             n_frames;
  char            filename_SSimPL_root[256];
  char            filename_halos_version[256];
  char            filename_trees_version[256];
  char            filename_out_dir[256];
  char            snap_filename_root[256];
  char            mark_filename_root[256];
  char            smooth_filename_root[256];
  char            snap_a_list_filename[256];
  int             n_snap_a_list;
  int             snap_number;
  double         *snap_a_list;
  int            *snap_list;
  int             flag_comoving;
  int             flag_fade;
  double          alpha_fade;
  int             flag_force_periodic;
  int             flag_read_marked;
  int             flag_add_absorption;
  plist_info    **plist_list;
  tree_info      *trees;
  mark_arg_info  *mark_arg_first;
  mark_arg_info  *mark_arg_last;
  double          h_Hubble;
  double          f_absorption;
  int             w_mode;
  int             v_mode;
  int             sealed; // TRUE if the render is fully initialized
};

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif
void init_perspective(perspective_info **perspective,int mode);
void free_perspective(perspective_info **perspective);
void copy_perspective(perspective_info *from,perspective_info *to);
void init_perspective_interp(perspective_interp_info **perspective_interp);
void free_perspective_interp(perspective_interp_info **perspective_interp);
void add_scene_perspective(scene_info *scene);
void init_scene(scene_info **scene);
void seal_scenes(scene_info *scenes);
void free_scenes(scene_info **scene);
void init_camera(camera_info **camera, int mode);
void seal_render_camera(render_info *render);
void free_camera(camera_info **camera);
void init_render(render_info **render);
void free_render(render_info **render);
void add_render_scene(render_info *render);
void seal_render(render_info *render);
int  set_render_state(render_info *render,int frame,int mode);
void parse_render_file(render_info **render, char *filename);
void write_frame(render_info *render,int frame,int mode);
void read_frame(render_info *render,int frame);
void set_frame(camera_info *camera);
void set_render_scale(render_info *render,double RGB_min,double RGB_max,double Y_min,double Y_max,double Z_min,double Z_max);
int  set_transfer_function(char *line,int i_word,interp_info **return_interp);

void render_frame(render_info  *render);

void open_movie(char       *filename,
                int         width,
                int         height,
                int         n_frames,
                int         rate,
                movie_info *movie);
void write_image_to_movie(image_info *image, movie_info *movie);              
void close_movie(movie_info *movie);
void init_image(int          width,
                int          height,
                int          colourmapselect,
                image_info **image);
void free_image(image_info **image);
void create_colour_table(int     colourmapselect,
                         int     n_colours,
                         int  ***colour_table);

void set_image_RGB(image_info *image,
                   double      image_min,
                   double      image_max);
void set_image_RGBY(image_info *image_RGBY_in,
                    image_info *image_RGB_in,
                    image_info *image_Y_in,
                    double      RGB_min,
                    double      RGB_max,
                    double      Y_min,
                    double      Y_max);
void set_image_RGBY_3CHANNEL(image_info *image_RGBY_3CHANNEL_in,
                             image_info *image_RY_in,
                             image_info *image_GY_in,
                             image_info *image_BY_in,
                             image_info *image_Y_in,
                             double      Y_min,
                             double      Y_max);

void write_image(image_info *image,const char *filename_dir,const char *filename_root,int mode);
void read_image (image_info *image,const char *filename_dir,const char *filename_root);

void read_gadget_binary_render(char       *filename_root_in,
                               int         snapshot_number,
                               plist_info *plist,
                               int         mode);

void set_sph_kernel(double **kernel_radius,
                    double **kernel_table_3d,
                    double **kernel_table_2d,
                    double  *kernel_table_2d_average,
                    int      mode);

void add_mark_argument   (render_info *render,const char *species,int value,const char *type,...);
void create_mark_argument(render_info *render,mark_arg_info **new_arg);
void free_mark_arguments(mark_arg_info **argument);
void perform_marking     (render_info *render);
void pick_best_snap(double a_search,double *snap_a_list,int n_snap_a_list,int *snap_best,double *snap_diff_best);
void process_SSimPL_halos(render_info *render,
                          int          i_snap,
                          int          i_pass,
                          int          mode,
                          int          select_function(int                i_group,
                                                       int                j_subgroup,
                                                       int                i_subgroup,
                                                       int                flag_long_ids,
                                                       process_halo_info *group_i,
                                                       process_halo_info *subgroup_i,
                                                       void              *params),
                          int          action_function(int                i_group,
                                                       int                j_subgroup,
                                                       int                i_subgroup,
                                                       int                flag_long_ids,
                                                       process_halo_info *group_i,
                                                       process_halo_info *subgroup_i,
                                                       void              *params),
                          void        *params);
#ifdef __cplusplus
}
#endif

#endif

