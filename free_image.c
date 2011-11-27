#include <string.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpRender.h>

void free_image(image_info **image){
  if((*image)!=NULL){
     (*image)->width   =0;
     (*image)->height  =0;
     (*image)->n_pixels=0;
     SID_free(SID_FARG (*image)->values);
     gdImageDestroy((*image)->gd_ptr);
  
     // Free colour table
     SID_free(SID_FARG (*image)->colour_table[0]);
     SID_free(SID_FARG (*image)->colour_table[1]);
     SID_free(SID_FARG (*image)->colour_table[2]);
     SID_free(SID_FARG (*image)->colour_table);

     // Free the image structure
     SID_free((void **)image);
  }
}
