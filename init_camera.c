#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void init_camera(camera_info **camera, int mode){
  SID_log("Initializing camera...",SID_LOG_OPEN);

  // Initialize camera
  (*camera)=(camera_info *)SID_malloc(sizeof(camera_info));
  (*camera)->camera_mode=mode;

  // Initialize the perspective information for this camera
  init_perspective(&((*camera)->perspective),RENDER_INIT_PERSPECTIVE);

  // Initialze image information
  (*camera)->flag_velocity_space=FALSE;
  (*camera)->flag_calc_Z_image  =FALSE;
  (*camera)->stereo_ratio       =0.;
  (*camera)->RGB_mode           =CAMERA_RGB_MODE_DEFAULT;
  (*camera)->RGB_gamma          =NULL;
  (*camera)->transfer_list      =NULL;
  (*camera)->Y_mode             =0;
  (*camera)->Y_gamma            =NULL;
  (*camera)->Z_mode             =0;
  (*camera)->Z_gamma            =NULL;
  (*camera)->f_near_field       =0.;
  (*camera)->f_taper_field      =0.;
  (*camera)->f_image_plane      =1.;
  strcpy((*camera)->RGB_param,"");
  strcpy((*camera)->Y_param,  "");

  // Initialize image buffers
  (*camera)->colour_table             =9;
  (*camera)->image_RGB                =NULL;
  (*camera)->image_Y                  =NULL;
  (*camera)->image_Z                  =NULL;
  (*camera)->image_RGBY               =NULL;
  (*camera)->image_RGBY_3CHANNEL      =NULL;
  (*camera)->image_RGB_left           =NULL;
  (*camera)->image_Y_left             =NULL;
  (*camera)->image_Z_left             =NULL;
  (*camera)->image_RGBY_left          =NULL;
  (*camera)->image_RGBY_3CHANNEL_left =NULL;
  (*camera)->image_RGB_right          =NULL;
  (*camera)->image_Y_right            =NULL;
  (*camera)->image_Z_right            =NULL;
  (*camera)->image_RGBY_right         =NULL;
  (*camera)->image_RGBY_3CHANNEL_right=NULL;
  (*camera)->mask_RGB                 =NULL;
  (*camera)->mask_Y                   =NULL;
  (*camera)->mask_RGBY                =NULL;
  (*camera)->mask_RGBY_3CHANNEL       =NULL;
  (*camera)->mask_RGB_left            =NULL;
  (*camera)->mask_Y_left              =NULL;
  (*camera)->mask_RGBY_left           =NULL;
  (*camera)->mask_RGBY_3CHANNEL_left  =NULL;
  (*camera)->mask_RGB_right           =NULL;
  (*camera)->mask_Y_right             =NULL;
  (*camera)->mask_RGBY_right          =NULL;
  (*camera)->mask_RGBY_3CHANNEL_right =NULL;
  (*camera)->image_RY                 =NULL;
  (*camera)->image_RY_left            =NULL;
  (*camera)->image_RY_right           =NULL;
  (*camera)->image_GY                 =NULL;
  (*camera)->image_GY_left            =NULL;
  (*camera)->image_GY_right           =NULL;
  (*camera)->image_BY                 =NULL;
  (*camera)->image_BY_left            =NULL;
  (*camera)->image_BY_right           =NULL;

  SID_log("Done.",SID_LOG_CLOSE);
}

