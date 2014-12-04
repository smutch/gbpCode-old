#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void set_frame(camera_info *camera){
  int         i_x,i_y,i_pixel;
  int         i_image;
  double      image_min,image_max;
  double      image_range;
  char        pixel_value;
  double      brightness;
  image_info *image_RGB;
  image_info *image_Y;
  image_info *image_Z;
  image_info *image_RGBY;
  image_info *image_RY;
  image_info *image_GY;
  image_info *image_BY;
  image_info *image;
  double     *values;

  // Loop over each set of images
  for(i_image=0;i_image<3;i_image++){
    switch(i_image){
    case 0:
      image_RGB =camera->image_RGB;
      image_Y   =camera->image_Y;
      image_Z   =camera->image_Z;
      image_RGBY=camera->image_RGBY;
      image_RY  =camera->image_RY;
      image_GY  =camera->image_GY;
      image_BY  =camera->image_BY;
      break;
    case 1:
      image_RGB =camera->image_RGB_left;
      image_Y   =camera->image_Y_left;
      image_Z   =camera->image_Z_left;
      image_RGBY=camera->image_RGBY_left;
      image_RY  =camera->image_RY_left;
      image_GY  =camera->image_GY_left;
      image_BY  =camera->image_BY_left;
      break;
    case 2:
      image_RGB =camera->image_RGB_right;
      image_Y   =camera->image_Y_right;
      image_Z   =camera->image_Z_right;
      image_RGBY=camera->image_RGBY_right;
      image_RY  =camera->image_RY_right;
      image_GY  =camera->image_GY_right;
      image_BY  =camera->image_BY_right;
      break;
    }

    // Set RGB Image
    if(image_RGB!=NULL)
       set_image_RGB(image_RGB,camera->RGB_range[0],camera->RGB_range[1]);

    // Set Y Image
    if(image_Y!=NULL)
       set_image_RGB(image_Y,camera->Y_range[0],camera->Y_range[1]);

    // Set Z Image
    if(image_Z!=NULL)
       set_image_RGB(image_Z,camera->Z_range[0],camera->Z_range[1]);

    // Set RGBY Image
    if(check_mode_for_flag(camera->RGB_mode,CAMERA_RGB_MODE_TABLE)){
       // n.b.: This check may fail deliberately (generally, each case has only a few sets of images defined).  Don't throw an error when it does.
       if(image_RGBY!=NULL && image_RGB!=NULL && image_Y!=NULL)
          set_image_RGBY(image_RGBY,
                         image_RGB,
                         image_Y,
                         camera->RGB_range[0],
                         camera->RGB_range[1],
                         camera->Y_range[0],
                         camera->Y_range[1]);
    }
    else if(check_mode_for_flag(camera->RGB_mode,CAMERA_RGB_MODE_NOTABLE)){
       // n.b.: This check may fail deliberately (generally, each case has only a few sets of images defined).  Don't throw an error when it does.
       if(image_RGBY!=NULL && image_RY!=NULL && image_GY!=NULL && image_BY!=NULL)
          set_image_RGBY_no_table(image_RGBY,
                                  image_RY,
                                  image_GY,
                                  image_BY,
                                  camera->Y_range[0],
                                  camera->Y_range[1]);
    }
    else
       SID_trap_error("Invalid camera RGB mode (%d) specified in set_frame().",ERROR_LOGIC,camera->RGB_mode);
  }
}

