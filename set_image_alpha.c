#include <gbpLib.h>
#include <gbpRender.h>
void set_image_alpha(image_info *image,
                     double     *image_in,
                     double      image_min,
                     double      image_max){
  int    i_x,i_y,i_pixel;
  double image_range;
  double brightness;
  if(SID.I_am_Master){
    SID_log("Scaling alpha to the range %le->%le...",SID_LOG_OPEN,image_min,image_max);
    image_range=image_max-image_min;
    for(i_x=0,i_pixel=0;i_x<image->width;i_x++){
      for(i_y=0;i_y<image->height;i_y++,i_pixel++){
        // Compute the brightness of the pixel [0.->1.]
        brightness            =(image_in[i_pixel]-image_min)/image_range;
        brightness            =MAX(0.,MIN(brightness,1.));
        // Apply to colours
        image->red[i_pixel]   =(int)(brightness*(double)image->red[i_pixel]);
        image->green[i_pixel] =(int)(brightness*(double)image->green[i_pixel]);
        image->blue[i_pixel]  =(int)(brightness*(double)image->blue[i_pixel]);
        // Make sure colours are still in the range [0->255]
        image->red[i_pixel]   =MAX(0,MIN(image->red[i_pixel],  255));
        image->blue[i_pixel]  =MAX(0,MIN(image->blue[i_pixel], 255));
        image->green[i_pixel] =MAX(0,MIN(image->green[i_pixel],255));
        // Apply changes
        gdImageSetPixel(image->gd_ptr,i_x,i_y,gdTrueColorAlpha(image->red[i_pixel],image->green[i_pixel],image->blue[i_pixel],image->alpha[i_pixel]));
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
}
