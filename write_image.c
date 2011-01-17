#include <stdio.h>
#include <stdlib.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpRender.h>

void write_image(image_info *image,char *filename_root,int mode){
  char  filename[256];
  int   i_pixel;
  FILE *fp;

  // Write raw image
  if(!check_mode_for_flag(mode,WRITE_IMAGE_PNG_ONLY)){
    sprintf(filename,"%s.raw",filename_root);
    fp=fopen(filename,"w");
    fwrite(&(image->width), sizeof(int),1,fp);
    fwrite(&(image->height),sizeof(int),1,fp);
    fwrite(image->values,sizeof(double),image->n_pixels,fp);
    fclose(fp);
  }

  // Write png
  sprintf(filename,"%s.png",filename_root);
  fp=fopen(filename,"w");
  gdImagePng(image->gd_ptr,fp);
  fclose(fp);

}
