#include <gbpLib.h>
#include <gbpRender.h>
#include <gd.h>

void set_image_RGB(image_info *image,
                   double     *image_in,
                   double      image_min,
                   double      image_max){
  int    i_x,i_y,i_pixel;
  double image_range;
  int    pixel_value;
  int    red,green,blue,colour,alpha;
  if(SID.I_am_Master){
    SID_log("Scaling RGB to the range %le->%le...",SID_LOG_OPEN,image_min,image_max);
    image_range=image_max-image_min;
    for(i_x=0,i_pixel=0;i_x<image->width;i_x++){
      for(i_y=0;i_y<image->height;i_y++,i_pixel++){
        pixel_value          =(int)(255.*(image_in[i_pixel]-image_min)/image_range);
        pixel_value          =MAX(0,MIN(pixel_value,255));
        image->red[i_pixel]  =image->colour_table[0][pixel_value];
        image->green[i_pixel]=image->colour_table[1][pixel_value];
        image->blue[i_pixel] =image->colour_table[2][pixel_value];
        gdImageSetPixel(image->gd_ptr,i_x,i_y,gdTrueColorAlpha(image->red[i_pixel],image->green[i_pixel],image->blue[i_pixel],image->alpha[i_pixel]));
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
}
