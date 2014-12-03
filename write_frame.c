#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void write_frame(render_info *render,int frame,int mode){
  char filename_RGB[256];
  char filename_Y[256];
  char filename_Z[256];
  char filename_RGBY[256];
  char filename_RY[256];
  char filename_GY[256];
  char filename_BY[256];

  SID_log("Writing rendered frame...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Create directory if needed
  mkdir(render->filename_out_dir,02755);

  // Write mono-images
  set_frame(render->camera);

  sprintf(filename_RGB, "%s/RGB_M_%05d", render->filename_out_dir,frame);
  sprintf(filename_Y,   "%s/Y_M_%05d",   render->filename_out_dir,frame);
  sprintf(filename_Z,   "%s/Z_M_%05d",   render->filename_out_dir,frame);
  sprintf(filename_RGBY,"%s/RGBY_M_%05d",render->filename_out_dir,frame);
  sprintf(filename_RY,  "%s/RY_M_%05d",  render->filename_out_dir,frame);
  sprintf(filename_GY,  "%s/GY_M_%05d",  render->filename_out_dir,frame);
  sprintf(filename_BY,  "%s/BY_M_%05d",  render->filename_out_dir,frame);
  if(check_mode_for_flag(render->camera->RGB_mode,CAMERA_RGB_MODE_TABLE)){
     write_image(render->camera->image_RGB, filename_RGB, mode);
     write_image(render->camera->image_Y,   filename_Y,   mode);
  }
  else if(check_mode_for_flag(render->camera->RGB_mode,CAMERA_RGB_MODE_NOTABLE)){
     write_image(render->camera->image_RY,filename_RY,mode);
     write_image(render->camera->image_GY,filename_GY,mode);
     write_image(render->camera->image_BY,filename_BY,mode);
  }
  else
     SID_trap_error("Invalid camera RGB mode (%d) specified in write_frame().",ERROR_LOGIC,render->camera->RGB_mode);
  write_image(render->camera->image_RGBY,filename_RGBY,mode);
  write_image(render->camera->image_Z,   filename_Z,   mode);

  // Write stereo-images
  if(check_mode_for_flag(render->camera->camera_mode,CAMERA_STEREO)){
    // Left
    sprintf(filename_RGB, "%s/RGB_L_%05d", render->filename_out_dir,frame);
    sprintf(filename_Y,   "%s/Y_L_%05d",   render->filename_out_dir,frame);
    sprintf(filename_Z,   "%s/Z_L_%05d",   render->filename_out_dir,frame);
    sprintf(filename_RGBY,"%s/RGBY_L_%05d",render->filename_out_dir,frame);
    sprintf(filename_RY,  "%s/RY_L_%05d",  render->filename_out_dir,frame);
    sprintf(filename_GY,  "%s/GY_L_%05d",  render->filename_out_dir,frame);
    sprintf(filename_BY,  "%s/BY_L_%05d",  render->filename_out_dir,frame);
    if(check_mode_for_flag(render->camera->RGB_mode,CAMERA_RGB_MODE_TABLE)){
       write_image(render->camera->image_RGB_left,  filename_RGB, mode);
       write_image(render->camera->image_Y_left,    filename_Y,   mode);
    }
    else if(check_mode_for_flag(render->camera->RGB_mode,CAMERA_RGB_MODE_NOTABLE)){
       write_image(render->camera->image_RY,filename_RY,mode);
       write_image(render->camera->image_GY,filename_GY,mode);
       write_image(render->camera->image_BY,filename_BY,mode);
    }
    else
       SID_trap_error("Invalid camera RGB mode (%d) specified in write_frame().",ERROR_LOGIC,render->camera->RGB_mode);
    write_image(render->camera->image_RGBY_left, filename_RGBY,mode);
    write_image(render->camera->image_Z_left,    filename_Z,   mode);

    // Right
    sprintf(filename_RGB, "%s/RGB_R_%05d", render->filename_out_dir,frame);
    sprintf(filename_Y,   "%s/Y_R_%05d",   render->filename_out_dir,frame);
    sprintf(filename_Z,   "%s/Z_R_%05d",   render->filename_out_dir,frame);
    sprintf(filename_RGBY,"%s/RGBY_R_%05d",render->filename_out_dir,frame);
    sprintf(filename_RY,  "%s/RY_R_%05d",  render->filename_out_dir,frame);
    sprintf(filename_GY,  "%s/GY_R_%05d",  render->filename_out_dir,frame);
    sprintf(filename_BY,  "%s/BY_R_%05d",  render->filename_out_dir,frame);
    if(check_mode_for_flag(render->camera->RGB_mode,CAMERA_RGB_MODE_TABLE)){
       write_image(render->camera->image_RGB_right, filename_RGB, mode);
       write_image(render->camera->image_Y_right,   filename_Y,   mode);
    }
    else if(check_mode_for_flag(render->camera->RGB_mode,CAMERA_RGB_MODE_NOTABLE)){
       write_image(render->camera->image_RY,filename_RY,mode);
       write_image(render->camera->image_GY,filename_GY,mode);
       write_image(render->camera->image_BY,filename_BY,mode);
    }
    else
       SID_trap_error("Invalid camera RGB mode (%d) specified in write_frame().",ERROR_LOGIC,render->camera->RGB_mode);
    write_image(render->camera->image_RGBY_right,filename_RGBY,mode);
    write_image(render->camera->image_Z_right,   filename_Z,   mode);
  }
  SID_log("Done.",SID_LOG_CLOSE);
}

