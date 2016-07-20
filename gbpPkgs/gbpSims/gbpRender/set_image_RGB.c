#include <gbpLib.h>
#include <gbpRender.h>
#include <gd.h>

void set_image_RGB(image_info *image,
                   double      image_min,
                   double      image_max){
  double *image_in;
  int     i_x,i_y,i_pixel;
  double  image_range;
  int     pixel_value;
  int     red,green,blue,alpha;
  int     n_colours;

  image_range=image_max-image_min;
  image_in   =image->values;
  n_colours  =image->n_colours-1;
  for(i_x=0,i_pixel=0;i_x<image->width;i_x++){
    for(i_y=0;i_y<image->height;i_y++,i_pixel++){
      pixel_value=(int)((double)n_colours*(image_in[i_pixel]-image_min)/image_range);
      pixel_value=MAX(0,MIN(pixel_value,n_colours));
      red        =(int)image->colour_table[0][pixel_value];
      green      =(int)image->colour_table[1][pixel_value];
      blue       =(int)image->colour_table[2][pixel_value];
      alpha      =gdAlphaOpaque;
      gdImageSetPixel(image->gd_ptr,i_x,i_y,gdTrueColorAlpha(red,green,blue,alpha));
    }
  }
}
