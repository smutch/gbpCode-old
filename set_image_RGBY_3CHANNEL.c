#include <gbpLib.h>
#include <gbpRender.h>
#include <gd.h>

void set_image_RGBY_3CHANNEL(image_info *image_RGBY_in,
                             image_info *image_RY_in,
                             image_info *image_GY_in,
                             image_info *image_BY_in,
                             double      Y_min,
                             double      Y_max){
  double *image_RY;
  double *image_GY;
  double *image_BY;
  int     i_x,i_y,i_pixel;
  int     n_colours;
  double  Y_range;
  int     red,green,blue,alpha;
  
  image_RY =image_RY_in->values;
  image_GY =image_GY_in->values;
  image_BY =image_BY_in->values;
  Y_range  =Y_max-Y_min;
  n_colours=image_RGBY_in->n_colours-1;
  for(i_x=0,i_pixel=0;i_x<image_RGBY_in->width;i_x++){
    for(i_y=0;i_y<image_RGBY_in->height;i_y++,i_pixel++){
       red  =MAX(0,MIN((int)((double)n_colours*(image_RY[i_pixel]-Y_min)/Y_range),n_colours));
       green=MAX(0,MIN((int)((double)n_colours*(image_GY[i_pixel]-Y_min)/Y_range),n_colours));
       blue =MAX(0,MIN((int)((double)n_colours*(image_BY[i_pixel]-Y_min)/Y_range),n_colours));
       alpha=gdAlphaOpaque;
       gdImageSetPixel(image_RGBY_in->gd_ptr,i_x,i_y,gdTrueColorAlpha(red,green,blue,alpha));
    }
  }
}
