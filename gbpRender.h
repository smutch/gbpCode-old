#ifndef RENDER_AWAKE
#define RENDER_AWAKE
#include <stdint.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gd.h>
#if USE_FFMPEG
  #include <avformat.h>
  #include <swscale.h>
#endif

#define CAMERA_MONO            0
#define CAMERA_STEREO          2
#define CAMERA_PLANE_PARALLEL  4
#define CAMERA_DEFAULT         CAMERA_MONO

#define RENDER_INIT_PERSPECTIVE 1
#define RENDER_INIT_EVOLVE      2
#define RENDER_INIT_DEFAULT     RENDER_INIT_PERSPECTIVE

#define SET_RENDER_DEFAULT 2
#define SET_RENDER_RESCALE 4

#define WRITE_FRAME_DEFAULT  2
#define WRITE_FRAME_PNG_ONLY 4

#define RENDER_SWS_FLAGS      SWS_BICUBIC
#define RENDER_FORMAT_DEFAULT "mpeg"
#define RENDER_AV_LOG_LEVEL   -1
#define RENDER_GOP_SIZE        2
#define RENDER_QSCALE          1

#define WRITE_IMAGE_DEFAULT  2
#define WRITE_IMAGE_PNG_ONLY 4

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
  double p_o[3];  // Object position
  double p_c[3];  // Camera position; computed from p_o, radius zeta and theta
  double d_o;     // Separation between object and camera
  double radius;  // Separation between object and camera modulo phi
  double theta;   // Rotation angle (azimuthal)   of camera position about x,y,z_o (after c-o transformation)
  double zeta;    // Rotation angle (altitudinal) of camera position about x,y,z_o (after c-o transformation)
  double phi;     // Zoom factor
  double FOV;     // Field of view at object position
  double time;    // Used for time-evolving movies
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
  int               width;
  int               height;
  double            stereo_ratio;
  perspective_info *perspective; // Present perspective state of camera
  image_info       *image_RGBY;
  image_info       *image_RGBY_left;
  image_info       *image_RGBY_right;
  char             *mask_RGBY;
  char             *mask_RGBY_left;
  char             *mask_RGBY_right;
  int               RGB_mode;
  char              RGB_param[64];
  double            RGB_range[2];
  interp_info      *RGB_gamma;
  ADaPS            *RGB_transfer;
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
  ADaPS            *Y_transfer;
  image_info       *image_Y;
  image_info       *image_Y_left;
  image_info       *image_Y_right;
  char             *mask_Y;
  char             *mask_Y_left;
  char             *mask_Y_right;
  int               Z_mode;
  double            Z_range[2];
  interp_info      *Z_gamma;
  ADaPS            *Z_transfer;
  image_info       *image_Z;
  image_info       *image_Z_left;
  image_info       *image_Z_right;
};

typedef struct render_info render_info;
struct render_info{
  camera_info *camera;
  scene_info  *scenes;
  scene_info  *first_scene;
  scene_info  *last_scene;
  int          n_frames;
  char         filename_out_dir[256];
  char         snap_filename_root[256];
  char         mark_filename_root[256];
  char         smooth_filename_root[256];
  char         snap_a_list_filename[256];
  int          n_snap_a_list;
  double      *snap_a_list;
  int          snap_number;
  int          snap_number_read;
  int          flag_comoving;
  int          flag_force_periodic;
  int          flag_read_marked;
  plist_info   plist;
  double       h_Hubble;
  double       near_field;
  int          mode;
  int          sealed; // TRUE if the render is fully initialized
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

void write_image(image_info *image,char *filename,int mode);
void read_image(image_info *image,char *filename_root);

#ifdef __cplusplus
}
#endif

#endif

