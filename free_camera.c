#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void free_camera(camera_info **camera){
  SID_log("Freeing camera...",SID_LOG_OPEN);
  free_perspective(&((*camera)->perspective));
  free_image(&((*camera)->image_RGB));
  free_image(&((*camera)->image_Y));
  free_image(&((*camera)->image_Z));
  free_image(&((*camera)->image_RGBY));
  SID_free(SID_FARG (*camera)->mask_RGB);
  SID_free(SID_FARG (*camera)->mask_Y);
  SID_free(SID_FARG (*camera)->mask_RGBY);
  if(check_mode_for_flag((*camera)->camera_mode,CAMERA_STEREO)){
    free_image(&((*camera)->image_RGB_left));
    free_image(&((*camera)->image_Y_left));
    free_image(&((*camera)->image_Z_left));
    free_image(&((*camera)->image_RGBY_left));
    free_image(&((*camera)->image_RGB_right));
    free_image(&((*camera)->image_Y_right));
    free_image(&((*camera)->image_Z_right));
    free_image(&((*camera)->image_RGBY_right));
    free_image(&((*camera)->image_RY));
    free_image(&((*camera)->image_RY_left));
    free_image(&((*camera)->image_RY_right));
    free_image(&((*camera)->image_GY));
    free_image(&((*camera)->image_GY_left));
    free_image(&((*camera)->image_GY_right));
    free_image(&((*camera)->image_BY));
    free_image(&((*camera)->image_BY_left));
    free_image(&((*camera)->image_BY_right));
    SID_free(SID_FARG (*camera)->mask_RGB_left);
    SID_free(SID_FARG (*camera)->mask_Y_left);
    SID_free(SID_FARG (*camera)->mask_RGBY_left);
    SID_free(SID_FARG (*camera)->mask_RGB_right);
    SID_free(SID_FARG (*camera)->mask_Y_right);
    SID_free(SID_FARG (*camera)->mask_RGBY_right);
  }
  if((*camera)->RGB_gamma!=NULL)
    free_interpolate(SID_FARG (*camera)->RGB_gamma,NULL);
  if((*camera)->transfer_list!=NULL)
    ADaPS_free(SID_FARG (*camera)->transfer_list);
  if((*camera)->Y_gamma!=NULL)
    free_interpolate(SID_FARG (*camera)->Y_gamma,NULL);
  if((*camera)->Z_gamma!=NULL)
    free_interpolate(SID_FARG (*camera)->Z_gamma,NULL);
  SID_free((void **)camera);
  SID_log("Done.",SID_LOG_CLOSE);
}

