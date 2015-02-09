#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpRender.h>

int main(int argc, char *argv[]){
  int         i_frame;
  int         width;
  int         height;
  int         colour_table;
  image_info *image;
  char        file_out[256];
  double     *image_array;
  FILE       *fp;
  int         i_x,i_y,i_pixel;

  SID_init(&argc,&argv,NULL,NULL);

  width       =atoi(argv[1]);
  height      =atoi(argv[2]);
  colour_table=atoi(argv[3]);
  strcpy(file_out,argv[4]);

  SID_log("Creating image of colourtable {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,file_out);

  // Initialize movie
  init_image(width,height,colour_table,&image);
  image_array=image->values;
  for(i_x=0,i_pixel=0;i_x<width;i_x++){
    for(i_y=0;i_y<height;i_y++,i_pixel++){
      image_array[i_pixel]=(double)(i_x+1)/(double)width;
    }
  }
  set_image_RGB(image,
                0.,
                1.);
  write_image(image,file_out,WRITE_IMAGE_DEFAULT);
  free_image(&image);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}
