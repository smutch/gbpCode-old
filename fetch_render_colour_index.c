#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

// WARNING: this code is very inefficient
char fetch_render_colour_index(render_info *render,const char *name){
  // Default in the case of a new entry
  char index_return=render->n_colour_list;

  // Check if the colour is already in the list
  int flag_found=FALSE;
  for(char i_list=0;i_list<render->n_colour_list;i_list++){
     if(!strcmp(name,render->colour_name[i_list])){
        index_return=i_list;
        flag_found=TRUE;
     }
  }

  // If we didn't find a colour already loaded, add it (if we can)
  if(!flag_found){
     // Don't allow n_list to exceed 256
     if(render->n_colour_list==256)
        SID_trap_error("The number of loaded colours in the render colour list can not exceed 256.",ERROR_LOGIC);

     // Open the file and look for the colour
     char filename_rgb[MAX_FILENAME_LENGTH];
     sprintf(filename_rgb,"%s/rgb.txt",GBP_DATA_DIR);
     FILE *fp=fopen(filename_rgb,"r");
     int   n_colours_file=count_lines_data(fp);
     size_t line_length=0;
     char  *line       =NULL;
     int   flag_success=FALSE;
     int   R_i,G_i,B_i;
     for(int i_file=0;i_file<n_colours_file && !flag_success;i_file++){
        char colour_name_i[32];
        grab_next_line_data(fp,&line,&line_length);
        grab_tail(line,4,colour_name_i);
        if(!strcmp(colour_name_i,name)){
           grab_int(line,1,&R_i);
           grab_int(line,2,&G_i);
           grab_int(line,3,&B_i);
           flag_success=TRUE;
        }
     }
     fclose(fp);
     if(!flag_success)
        SID_trap_error("Could not initialize colour {%s} in fetch_render_colour_index.",ERROR_LOGIC,name);
     else{
        // This is the bit that's poor: realloc the arrays and add the new entry
        char  **colour_name_old =render->colour_name;
        int   **colour_RGB_old  =render->colour_RGB;
        float **colour_f_RGB_old=render->colour_f_RGB;
        render->n_colour_list++;
        render->colour_name =(char  **)SID_malloc(sizeof(char  *)*render->n_colour_list);
        render->colour_RGB  =(int   **)SID_malloc(sizeof(int   *)*render->n_colour_list);
        render->colour_f_RGB=(float **)SID_malloc(sizeof(float *)*render->n_colour_list);
        for(int i_list=0;i_list<render->n_colour_list;i_list++){
           render->colour_name[i_list] =(char  *)SID_malloc(sizeof(char )*32);
           render->colour_RGB[i_list]  =(int   *)SID_malloc(sizeof(int  )*3);
           render->colour_f_RGB[i_list]=(float *)SID_malloc(sizeof(float)*3);
           if(i_list==(render->n_colour_list-1)){
              strcpy(render->colour_name[i_list],name);
              render->colour_RGB[i_list][0]  =R_i;
              render->colour_RGB[i_list][1]  =G_i;
              render->colour_RGB[i_list][2]  =B_i;
              render->colour_f_RGB[i_list][0]=((float)R_i)/256.;
              render->colour_f_RGB[i_list][1]=((float)G_i)/256.;
              render->colour_f_RGB[i_list][2]=((float)B_i)/256.;
           }
           else{
             memcpy(render->colour_name[i_list], colour_name_old[i_list], sizeof(char )*32); 
             memcpy(render->colour_RGB[i_list],  colour_RGB_old[i_list],  sizeof(int  )*3); 
             memcpy(render->colour_f_RGB[i_list],colour_f_RGB_old[i_list],sizeof(float)*3);
             SID_free(SID_FARG colour_name_old[i_list]);
             SID_free(SID_FARG colour_RGB_old[i_list]);
             SID_free(SID_FARG colour_f_RGB_old[i_list]);
           }
        }
        SID_free(SID_FARG colour_name_old);
        SID_free(SID_FARG colour_RGB_old);
        SID_free(SID_FARG colour_f_RGB_old);
     }
  }
  return(index_return);
}

