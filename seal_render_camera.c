#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void seal_render_camera(render_info *render){

  SID_log("Sealing camera...",SID_LOG_OPEN);

  // Decide if we're producing a stereo image based on stereo_factor
  if(render->camera->stereo_ratio>0.)
     render->camera->camera_mode|=CAMERA_STEREO;

  // Initialize the perspective information for this camera
  copy_perspective(render->first_scene->first_perspective,render->camera->perspective);

  // Initialize image buffers (use colour_table=1 ... ie B&W ... for Y and Z images)
  if(check_mode_for_flag(render->camera->camera_mode,CAMERA_STEREO)){
    // LEFT
    int flag_images_defined=FALSE;
    if(check_mode_for_flag(render->camera->RGB_mode,CAMERA_RGB_MODE_1CHANNEL)){
       init_image(render->camera->width,render->camera->height,render->camera->colour_table,&(render->camera->image_RGB_left));
       init_image(render->camera->width,render->camera->height,1,                           &(render->camera->image_Y_left));
       init_image(render->camera->width,render->camera->height,render->camera->colour_table,&(render->camera->image_RGBY_left));
       flag_images_defined=TRUE;
    }
    if(check_mode_for_flag(render->camera->RGB_mode,CAMERA_RGB_MODE_3CHANNEL)){
       init_image(render->camera->width,render->camera->height,2,&(render->camera->image_RY_left));
       init_image(render->camera->width,render->camera->height,3,&(render->camera->image_GY_left));
       init_image(render->camera->width,render->camera->height,4,&(render->camera->image_BY_left));
       init_image(render->camera->width,render->camera->height,1,&(render->camera->image_RGBY_3CHANNEL_left));
       flag_images_defined=TRUE;
    }
    if(!flag_images_defined)
       SID_trap_error("Invalid camera RGB mode (%d) specified in seal_render_camera().",ERROR_LOGIC,render->camera->RGB_mode);
    // RIGHT
    flag_images_defined=FALSE;
    if(check_mode_for_flag(render->camera->RGB_mode,CAMERA_RGB_MODE_1CHANNEL)){
       init_image(render->camera->width,render->camera->height,render->camera->colour_table,&(render->camera->image_RGB_right));
       init_image(render->camera->width,render->camera->height,1,                           &(render->camera->image_Y_right));
       init_image(render->camera->width,render->camera->height,render->camera->colour_table,&(render->camera->image_RGBY_right));
       flag_images_defined=TRUE;
    }
    if(check_mode_for_flag(render->camera->RGB_mode,CAMERA_RGB_MODE_3CHANNEL)){
       init_image(render->camera->width,render->camera->height,2,&(render->camera->image_RY_right));
       init_image(render->camera->width,render->camera->height,3,&(render->camera->image_GY_right));
       init_image(render->camera->width,render->camera->height,4,&(render->camera->image_BY_right));
       init_image(render->camera->width,render->camera->height,1,&(render->camera->image_RGBY_3CHANNEL_right));
       flag_images_defined=TRUE;
    }
    if(!flag_images_defined)
       SID_trap_error("Invalid camera RGB mode (%d) specified in seal_render_camera().",ERROR_LOGIC,render->camera->RGB_mode);
    // Z-image
    if(render->camera->flag_calc_Z_image){
       init_image(render->camera->width,render->camera->height,1,&(render->camera->image_Z_left));
       init_image(render->camera->width,render->camera->height,1,&(render->camera->image_Z_right));
    }
    /*
    render->camera->mask_RGB_left  =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
    render->camera->mask_Y_left    =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
    render->camera->mask_RGBY_left =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
    render->camera->mask_RGB_right =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
    render->camera->mask_Y_right   =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
    render->camera->mask_RGBY_right=(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
    */
  }
  else{
     int flag_images_defined=FALSE;
     if(check_mode_for_flag(render->camera->RGB_mode,CAMERA_RGB_MODE_1CHANNEL)){
        init_image(render->camera->width,render->camera->height,render->camera->colour_table,&(render->camera->image_RGB));
        init_image(render->camera->width,render->camera->height,1,                           &(render->camera->image_Y));
        init_image(render->camera->width,render->camera->height,render->camera->colour_table,&(render->camera->image_RGBY));
        flag_images_defined=TRUE;
     }
     if(check_mode_for_flag(render->camera->RGB_mode,CAMERA_RGB_MODE_3CHANNEL)){
        init_image(render->camera->width,render->camera->height,2,&(render->camera->image_RY));
        init_image(render->camera->width,render->camera->height,3,&(render->camera->image_GY));
        init_image(render->camera->width,render->camera->height,4,&(render->camera->image_BY));
        init_image(render->camera->width,render->camera->height,1,&(render->camera->image_RGBY_3CHANNEL));
        flag_images_defined=TRUE;
     }
     if(!flag_images_defined)
        SID_trap_error("Invalid camera RGB mode (%d) specified in seal_render_camera().",ERROR_LOGIC,render->camera->RGB_mode);
     if(render->camera->flag_calc_Z_image)
        init_image(render->camera->width,render->camera->height,1,&(render->camera->image_Z));
     /*
     render->camera->mask_RGB =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
     render->camera->mask_Y   =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
     render->camera->mask_RGBY=(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
     */
  }

  // Convert camera depth-range to Mpc/h
  if(render->camera->flag_calc_Z_image){
     render->camera->Z_range[0]*=M_PER_MPC/render->h_Hubble;
     render->camera->Z_range[1]*=M_PER_MPC/render->h_Hubble;
  }

  SID_log("Done.",SID_LOG_CLOSE);
}

