#include <string.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpRender.h>

void free_image(image_info **image){
  (*image)->width   =0;
  (*image)->height  =0;
  (*image)->n_pixels=0;
  SID_free(SID_FARG (*image)->red);
  SID_free(SID_FARG (*image)->green);
  SID_free(SID_FARG (*image)->blue);
  SID_free(SID_FARG (*image)->alpha);
  SID_free(SID_FARG (*image)->values);
  gdImageDestroy((*image)->gd_ptr);
  
  // Free colour table
  SID_free(SID_FARG (*image)->colour_table[0]);
  SID_free(SID_FARG (*image)->colour_table[1]);
  SID_free(SID_FARG (*image)->colour_table[2]);
  SID_free(SID_FARG (*image)->colour_table);

  SID_free((void **)image);
}
