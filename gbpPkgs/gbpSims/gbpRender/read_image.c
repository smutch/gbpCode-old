#include <stdio.h>
#include <stdlib.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpRender.h>

void read_image(image_info *image,const char *filename_dir,const char *filename_root){
  char  filename[256];
  int   i_pixel;
  FILE *fp;

  // Read raw image
  sprintf(filename,"%s/raw/%s.raw",filename_dir,filename_root);
  if((fp=fopen(filename,"r"))==NULL)
     SID_trap_error("Could not open image {%s}.",ERROR_IO_OPEN,filename);
  fread_verify(&(image->width), sizeof(int),1,fp);
  fread_verify(&(image->height),sizeof(int),1,fp);
  image->n_pixels=image->width*image->height;
  fread_verify(image->values,sizeof(double),image->n_pixels,fp);
  fclose(fp);

  // To Do: Add .png image read here
  //sprintf(filename,"%s/%s.png",filename_dir,filename_root);
  //fp=fopen(filename,"r");
  //fclose(fp);

}
