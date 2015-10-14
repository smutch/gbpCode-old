#include <gbpLib.h>
#include <gbpRender.h>
#include <gd.h>

void set_image_RGBY_3CHANNEL(image_info *image_RGBY_in,
                             image_info *image_RY_in,
                             image_info *image_GY_in,
                             image_info *image_BY_in,
                             image_info *image_Y_in,
                             double      Y_min,
                             double      Y_max){
  int     i_x,i_y,i_pixel;
  int     n_colours;
  double  Y_range;
  int     red,green,blue,alpha;
  
  double *image_RY =image_RY_in->values;
  double *image_GY =image_GY_in->values;
  double *image_BY =image_BY_in->values;
  double *image_Y  =image_Y_in->values;
  Y_range  =Y_max-Y_min;
  n_colours=image_RGBY_in->n_colours-1;
  for(i_x=0,i_pixel=0;i_x<image_RGBY_in->width;i_x++){
    for(i_y=0;i_y<image_RGBY_in->height;i_y++,i_pixel++){
       // Draw colours from the RY, GY and BY images
       red  =MAX(0,MIN((int)(((double)n_colours)*image_RY[i_pixel]),n_colours));
       green=MAX(0,MIN((int)(((double)n_colours)*image_GY[i_pixel]),n_colours));
       blue =MAX(0,MIN((int)(((double)n_colours)*image_BY[i_pixel]),n_colours));
       alpha=gdAlphaOpaque;

       // Compute the brightness of the pixel [0.->1.]
       double brightness;
       brightness =(image_Y[i_pixel]-Y_min)/Y_range;
       brightness =MAX(0.,MIN(brightness,1.));

       // Scale the pixel colours by the brightness
       red  =(int)(brightness*(double)red);
       green=(int)(brightness*(double)green);
       blue =(int)(brightness*(double)blue);

       gdImageSetPixel(image_RGBY_in->gd_ptr,i_x,i_y,gdTrueColorAlpha(red,green,blue,alpha));
    }
  }
}
