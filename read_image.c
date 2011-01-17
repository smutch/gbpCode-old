#include <stdio.h>
#include <stdlib.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpRender.h>

void read_image(image_info *image,char *filename_root){
  char  filename[256];
  int   i_pixel;
  FILE *fp;

  // Read raw image
  sprintf(filename,"%s.raw",filename_root);
  fp=fopen(filename,"r");
  fread(&(image->width), sizeof(int),1,fp);
  fread(&(image->height),sizeof(int),1,fp);
  image->n_pixels=image->width*image->height;
  fread(image->values,sizeof(double),image->n_pixels,fp);
  fclose(fp);

  // To Do: Add .png image read here
  sprintf(filename,"%s.png",filename_root);
  fp=fopen(filename,"r");
  fclose(fp);

}
