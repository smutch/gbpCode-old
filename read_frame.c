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
  char filename_RGB[MAX_FILENAME_LENGTH];
  char filename_Y[MAX_FILENAME_LENGTH];
  char filename_Z[MAX_FILENAME_LENGTH];
  char filename_RGBY[MAX_FILENAME_LENGTH];
  char filename_RY[MAX_FILENAME_LENGTH];
  char filename_GY[MAX_FILENAME_LENGTH];
  char filename_BY[MAX_FILENAME_LENGTH];
  char filename_RGBY_3CHANNEL[MAX_FILENAME_LENGTH];

  // Read stereo-images
  SID_log("Reading rendered frame...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(int i_set=0;i_set<3;i_set++){
    // Set image set filesnames
    image_info *image_RGB          =NULL;
    image_info *image_RGBY         =NULL;
    image_info *image_RGBY_3CHANNEL=NULL;
    image_info *image_Y            =NULL;
    image_info *image_Z            =NULL;
    image_info *image_RY           =NULL;
    image_info *image_GY           =NULL;
    image_info *image_BY           =NULL;
    char set_label[8];
    if(i_set==0){
       sprintf(set_label,"L");
       image_RGB          =render->camera->image_RGB_left;
       image_RGBY         =render->camera->image_RGBY_left;
       image_RGBY_3CHANNEL=render->camera->image_RGBY_3CHANNEL_left;
       image_Y            =render->camera->image_Y_left;
       image_Z            =render->camera->image_Z_left;
       image_RY           =render->camera->image_RY_left;
       image_GY           =render->camera->image_GY_left;
       image_BY           =render->camera->image_BY_left;
    }
    else if(i_set==1){
       sprintf(set_label,"R");
       image_RGB          =render->camera->image_RGB_right;
       image_RGBY         =render->camera->image_RGBY_right;
       image_RGBY_3CHANNEL=render->camera->image_RGBY_3CHANNEL_right;
       image_Y            =render->camera->image_Y_right;
       image_Z            =render->camera->image_Z_right;
       image_RY           =render->camera->image_RY_right;
       image_GY           =render->camera->image_GY_right;
       image_BY           =render->camera->image_BY_right;
    }
    else if(i_set==2){
       sprintf(set_label,"M");
       image_RGB          =render->camera->image_RGB;
       image_RGBY         =render->camera->image_RGBY;
       image_RGBY_3CHANNEL=render->camera->image_RGBY_3CHANNEL;
       image_Y            =render->camera->image_Y;
       image_Z            =render->camera->image_Z;
       image_RY           =render->camera->image_RY;
       image_GY           =render->camera->image_GY;
       image_BY           =render->camera->image_BY;
    }
    else
       SID_trap_error("Undefined image set index in read_frame().",ERROR_LOGIC);
    sprintf(filename_RGB,          "%s/RGB_%s_%05d",          render->filename_out_dir,set_label,frame);
    sprintf(filename_Y,            "%s/Y_%s_%05d",            render->filename_out_dir,set_label,frame);
    sprintf(filename_Z,            "%s/Z_%s_%05d",            render->filename_out_dir,set_label,frame);
    sprintf(filename_RGBY,         "%s/RGBY_%s_%05d",         render->filename_out_dir,set_label,frame);
    sprintf(filename_RGBY_3CHANNEL,"%s/RGBY_3CHANNEL_%s_%05d",render->filename_out_dir,set_label,frame);
    sprintf(filename_RY,           "%s/RY_%s_%05d",           render->filename_out_dir,set_label,frame);
    sprintf(filename_GY,           "%s/GY_%s_%05d",           render->filename_out_dir,set_label,frame);
    sprintf(filename_BY,           "%s/BY_%s_%05d",           render->filename_out_dir,set_label,frame);
    if(image_RGB!=NULL)
       read_image(image_RGB,filename_RGB);
    if(image_Y!=NULL)
       read_image(image_Y,filename_Y);
    if(image_Z!=NULL)
       read_image(image_Z,filename_Z);
    if(image_RGBY!=NULL)
       read_image(image_RGBY,filename_RGBY);
    if(image_RY!=NULL)
       read_image(image_RY,filename_RY);
    if(image_GY!=NULL)
       read_image(image_GY,filename_GY);
    if(image_BY!=NULL)
       read_image(image_BY,filename_BY);
    if(image_RGBY_3CHANNEL!=NULL)
       read_image(image_RGBY_3CHANNEL,filename_RGBY_3CHANNEL);
  }

  // Finalize things
  set_frame(render->camera);

  SID_log("Done.",SID_LOG_CLOSE);
}

