#include <stdio.h>
#include <stdlib.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpRender.h>

void write_image(image_info *image,const char *filename_dir,const char *filename_root,int mode){
  char  filename[256];
  int   i_pixel;
  FILE *fp;

  if(image!=NULL){

     // Compute image statistics
     double image_min= DBL_MAX;
     double image_max=-DBL_MAX;
     for(int i_w=0;i_w<image->width;i_w++){
        for(int i_h=0;i_h<image->height;i_h++){
           double image_i=image->values[i_h*image->width+i_w];
           if(image_i>image_max) image_max=image_i;
           if(image_i<image_min and image_i!=LOG_ZERO) image_min=image_i;
        }
     }
     //SID_log("Image range =%le -> %le",SID_LOG_COMMENT,image_min,image_max);

     // Write raw image
     if(check_mode_for_flag(mode,WRITE_IMAGE_PNG)){
       sprintf(filename,"%s/raw/%s.raw",filename_dir,filename_root);
       fp=fopen(filename,"w");
       fwrite(&(image->width), sizeof(int),1,fp);
       fwrite(&(image->height),sizeof(int),1,fp);
       fwrite(image->values,sizeof(double),image->n_pixels,fp);
       fclose(fp);
     }

     // Write png
     if(check_mode_for_flag(mode,WRITE_IMAGE_RAW)){
        sprintf(filename,"%s/%s.png",filename_dir,filename_root);
        fp=fopen(filename,"w");
        gdImagePng(image->gd_ptr,fp);
        fclose(fp);
     }
  }

}
