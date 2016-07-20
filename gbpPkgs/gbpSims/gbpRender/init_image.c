#include <string.h>
#include <gbpLib.h>
#include <gbpRender.h>

void init_image(int           width,
                int           height,
                int           colourmapselect,
                image_info **image){
  int i_pixel;

  // Allocate image
  (*image)=(image_info *)SID_malloc(sizeof(image_info));

  // Set image specs
  (*image)->width    =width;
  (*image)->height   =height;
  (*image)->n_pixels =width*height;

  // Init gdlib stuff
  (*image)->gd_ptr=gdImageCreateTrueColor((*image)->width,(*image)->height);

  // Allocate image buffers
  (*image)->values=(double *)malloc(sizeof(double)*(*image)->n_pixels);
  for(i_pixel=0;i_pixel<(*image)->n_pixels;i_pixel++)
    (*image)->values[i_pixel]=0.;

  // Set colourtable
  (*image)->n_colours=256;
  (*image)->colourmapselect=colourmapselect;
  create_colour_table(colourmapselect,(*image)->n_colours,&((*image)->colour_table));
}
