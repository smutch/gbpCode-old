#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void read_frame(render_info *render,int frame){
  char filename_RGB[256];
  char filename_Y[256];
  char filename_Z[256];
  char filename_RGBY[256];
  char filename_RY[256];
  char filename_GY[256];
  char filename_BY[256];

  // Read stereo-images
  SID_log("Reading rendered frame...",SID_LOG_OPEN|SID_LOG_TIMER);
  if(check_mode_for_flag(render->camera->camera_mode,CAMERA_STEREO)){
    // Left
    sprintf(filename_RGB, "%s/RGB_L_%05d", render->filename_out_dir,frame);
    sprintf(filename_Y,   "%s/Y_L_%05d",   render->filename_out_dir,frame);
    sprintf(filename_Z,   "%s/Z_L_%05d",   render->filename_out_dir,frame);
    sprintf(filename_RGBY,"%s/RGBY_L_%05d",render->filename_out_dir,frame);
    sprintf(filename_RY,  "%s/RY_L_%05d",  render->filename_out_dir,frame);
    sprintf(filename_GY,  "%s/GY_L_%05d",  render->filename_out_dir,frame);
    sprintf(filename_BY,  "%s/BY_L_%05d",  render->filename_out_dir,frame);
    if(render->camera->image_RGB_left!=NULL)
       read_image(render->camera->image_RGB_left,  filename_RGB);
    if(render->camera->image_Y_left!=NULL)
       read_image(render->camera->image_Y_left,    filename_Y);
    if(render->camera->image_Z_left!=NULL)
       read_image(render->camera->image_Z_left,    filename_Z);
    if(render->camera->image_RGBY_left!=NULL)
       read_image(render->camera->image_RGBY_left, filename_RGBY);
    if(render->camera->image_RY_left!=NULL)
       read_image(render->camera->image_RY_left, filename_RY);
    if(render->camera->image_GY_left!=NULL)
       read_image(render->camera->image_GY_left, filename_GY);
    if(render->camera->image_BY_left!=NULL)
       read_image(render->camera->image_BY_left, filename_BY);
    // Right
    sprintf(filename_RGB, "%s/RGB_R_%05d", render->filename_out_dir,frame);
    sprintf(filename_Y,   "%s/Y_R_%05d",   render->filename_out_dir,frame);
    sprintf(filename_Z,   "%s/Z_R_%05d",   render->filename_out_dir,frame);
    sprintf(filename_RGBY,"%s/RGBY_R_%05d",render->filename_out_dir,frame);
    sprintf(filename_RY,  "%s/RY_R_%05d",  render->filename_out_dir,frame);
    sprintf(filename_GY,  "%s/GY_R_%05d",  render->filename_out_dir,frame);
    sprintf(filename_BY,  "%s/BY_R_%05d",  render->filename_out_dir,frame);
    if(render->camera->image_RGB_right!=NULL)
       read_image(render->camera->image_RGB_right, filename_RGB);
    if(render->camera->image_Y_right!=NULL)
       read_image(render->camera->image_Y_right,   filename_Y);
    if(render->camera->image_Z_right!=NULL)
       read_image(render->camera->image_Z_right,   filename_Z);
    if(render->camera->image_RGBY_right!=NULL)
       read_image(render->camera->image_RGBY_right,filename_RGBY);
    if(render->camera->image_RY_right!=NULL)
       read_image(render->camera->image_RY_right,filename_RY);
    if(render->camera->image_GY_right!=NULL)
       read_image(render->camera->image_GY_right,filename_GY);
    if(render->camera->image_BY_right!=NULL)
       read_image(render->camera->image_BY_right,filename_BY);
  }
  else{
    // Mono
    sprintf(filename_RGB, "%s/RGB_M_%05d", render->filename_out_dir,frame);
    sprintf(filename_Y,   "%s/Y_M_%05d",   render->filename_out_dir,frame);
    sprintf(filename_Z,   "%s/Z_M_%05d",   render->filename_out_dir,frame);
    sprintf(filename_RGBY,"%s/RGBY_M_%05d",render->filename_out_dir,frame);
    sprintf(filename_RY,  "%s/RY_M_%05d",  render->filename_out_dir,frame);
    sprintf(filename_GY,  "%s/GY_M_%05d",  render->filename_out_dir,frame);
    sprintf(filename_BY,  "%s/BY_M_%05d",  render->filename_out_dir,frame);
    if(render->camera->image_RGB!=NULL)
       read_image(render->camera->image_RGB,filename_RGB);
    if(render->camera->image_Y!=NULL)
       read_image(render->camera->image_Y,filename_Y);
    if(render->camera->image_Z!=NULL)
       read_image(render->camera->image_Z,filename_Z);
    if(render->camera->image_RGBY!=NULL)
       read_image(render->camera->image_RGBY,filename_RGBY);
    if(render->camera->image_RY!=NULL)
       read_image(render->camera->image_RY,filename_RY);
    if(render->camera->image_GY!=NULL)
       read_image(render->camera->image_GY,filename_GY);
    if(render->camera->image_BY!=NULL)
       read_image(render->camera->image_BY,filename_BY);
  }

  set_frame(render->camera);

  SID_log("Done.",SID_LOG_CLOSE);
}

