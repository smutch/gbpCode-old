#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void init_render(render_info **render){
  int         flag_scatter;
  int         flag_no_velocities;
  SID_log("Initializing render structure...",SID_LOG_OPEN);

  // Allocate memory for rendering
  (*render)=(render_info *)SID_malloc(sizeof(render_info));

  // Initialize the camera
  init_camera(&((*render)->camera),CAMERA_DEFAULT);

  // Create the first scene
  init_scene(&((*render)->scenes));
  (*render)->first_scene=(*render)->scenes;
  (*render)->last_scene =(*render)->scenes;

  (*render)->n_frames           = 0;
  (*render)->n_interpolate      = 1;
  (*render)->n_snap_a_list      = 0;
  (*render)->snap_a_list        = NULL;
  (*render)->h_Hubble           = 1.;
  (*render)->f_absorption       =-1.;
  (*render)->flag_read_marked   = FALSE;
  (*render)->flag_comoving      = TRUE;
  (*render)->flag_fade          = FALSE;
  (*render)->alpha_fade         = 2.;
  (*render)->flag_force_periodic= FALSE;
  (*render)->flag_add_absorption= FALSE;
  (*render)->sealed             = FALSE;
  (*render)->v_mode             = MAKE_MAP_DEFAULT;
  (*render)->w_mode             = MAKE_MAP_DEFAULT;
  (*render)->plist_list         = NULL;
  (*render)->trees              = NULL;
  (*render)->mark_arg_first     = NULL;
  (*render)->mark_arg_last      = NULL;
  (*render)->kernel_radius      = NULL;
  (*render)->kernel_table       = NULL;
  (*render)->kernel_table_3d    = NULL;
  (*render)->kernel_table_avg   = 0.;
  (*render)->f_interpolate      = 0.;

  // Initialize SSimPL directory to something invalid
  sprintf((*render)->filename_SSimPL_root,"%s",RENDER_INVALID_SSIMPL_DIR);

  SID_log("Done.",SID_LOG_CLOSE);
}

