#include <gbpLib.h>
#include <gbpRender.h>
#include <gd.h>

void set_image_RGBY(image_info *image_RGBY_in,
                    image_info *image_RGB_in,
                    image_info *image_Y_in,
                    double      RGB_min,
                    double      RGB_max,
                    double      Y_min,
                    double      Y_max){
  double *image_RGBY;
  double *image_RGB;
  double *image_Y;
  int     i_x,i_y,i_pixel;
  double  RGB_range;
  double  Y_range;
  double  brightness;
  int     pixel_value;
  int     red,green,blue,alpha;
  int     n_colours;
  
  image_RGBY=image_RGBY_in->values;
  image_RGB =image_RGB_in->values;
  image_Y   =image_Y_in->values;
  RGB_range =RGB_max-RGB_min;
  Y_range   =Y_max  -Y_min;
  n_colours =image_RGB_in->n_colours-1;
  for(i_x=0,i_pixel=0;i_x<image_RGBY_in->width;i_x++){
    for(i_y=0;i_y<image_RGBY_in->height;i_y++,i_pixel++){

       // Set the pixel colour
       pixel_value=(int)(n_colours*(image_RGB[i_pixel]-RGB_min)/RGB_range);
       pixel_value=MAX(0,MIN(pixel_value,n_colours));
       red        =(int)image_RGB_in->colour_table[0][pixel_value];
       green      =(int)image_RGB_in->colour_table[1][pixel_value];
       blue       =(int)image_RGB_in->colour_table[2][pixel_value];
       alpha      =gdAlphaOpaque;

       // Compute the brightness of the pixel [0.->1.]
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
