#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void parse_render_file(render_info **render, char *filename){
  FILE   *fp;
  int     n_lines;
  int     i_line;
  char   *line=NULL;
  size_t  line_length=0;
  int     i_word,j_word;
  char    temp_word[1024];
  char    line_print[256];
  char    command[64];
  char    parameter[64];
  char    variable[64];
  char    image[64];
  char    c_value[64];
  double  d_value;
  int     i_value;
  int     i;
  int     flag;
  int     n_transfer;
  double      *transfer_array_x;
  double      *transfer_array_y;
  interp_info *temp_interp;
  FILE        *fp_list;
  
  SID_log("Parsing render file {%s}...",SID_LOG_OPEN,filename);
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
  if((fp=fopen(filename,"r"))==NULL)
    SID_trap_error("Could not open render file {%s}",ERROR_IO_OPEN,filename);
  init_render(render);
  sprintf((*render)->filename_out_dir,"%s.images/",filename);
  n_lines=count_lines_data(fp);
  for(i_line=0;i_line<n_lines;i_line++){
    SID_log("Processing line %d of %d...",SID_LOG_OPEN,i_line+1,n_lines);SID_Barrier(SID.COMM_WORLD);
    grab_next_line_data(fp,&line,&line_length);
    i_word=1;
    grab_word(line,i_word++,command);
    if(strlen(line)>0){
      sprintf(line_print,"%s",line);
      // Interpret command
      if(!strcmp(command,"set")){
        grab_word(line,i_word++,parameter);
        if(!strcmp(parameter,"mark")){
           char mark_species[32];
           int  mark_value;
           grab_int(line,i_word++, &mark_value);
           grab_word(line,i_word++,mark_species);
           if(i_word>count_words(line))
              add_mark_argument((*render),mark_species,mark_value,"");
           else{
              char mark_type[32];
              grab_word(line,i_word++,mark_type);
              if(!strcmp(mark_type,"sphere")){
                 double x,y,z,r;
                 grab_double(line,i_word++,&x);
                 grab_double(line,i_word++,&y);
                 grab_double(line,i_word++,&z);
                 grab_double(line,i_word++,&r);
                 add_mark_argument((*render),mark_species,mark_value,mark_type,x,y,z,r);
              }
              else if(!strcmp(mark_type,"group_index")){
                 int halo_index;
                 grab_int(line,i_word++,&halo_index);
                 add_mark_argument((*render),mark_species,mark_value,mark_type,halo_index);
              }
              else
                 SID_trap_error("Invalid mark type {%s} in line {%s}.",ERROR_LOGIC,mark_type,line);
           }
        }
        else if(!strcmp(parameter,"scene")){
          grab_word(line,i_word++,variable);
          if(!strcmp(variable,"n_frames")){
            grab_int(line,i_word++,&i_value);
            (*render)->last_scene->n_frames=i_value;
          }
          else if(!strcmp(variable,"FOV")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->FOV=d_value;
          }
          else if(!strcmp(variable,"p_o")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->p_o[0]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->p_o[1]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->p_o[2]=d_value;
          }
          else if(!strcmp(variable,"radius")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->radius=d_value;
          }
          else if(!strcmp(variable,"theta")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->theta=d_value;
          }
          else if(!strcmp(variable,"phi")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->phi=d_value;
          }
          else if(!strcmp(variable,"zeta")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->zeta=d_value;
          }
          else if(!strcmp(variable,"time")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->time=d_value;
          }
          else if(!strcmp(variable,"focus_shift_x")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->focus_shift_x=d_value;
          }
          else if(!strcmp(variable,"focus_shift_y")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->focus_shift_y=d_value;
          }
          else
            SID_trap_error("Unknown variable {%s} for parameter {%s} for command {%s} on line %d",ERROR_LOGIC,variable,parameter,command,i_line);
        }
        else if(!strcmp(parameter,"n_interpolate"))
          grab_int(line,i_word++,&((*render)->n_interpolate));
        else if(!strcmp(parameter,"SSimPL")){
          grab_word(line,i_word++,(*render)->filename_SSimPL_dir);
          grab_word(line,i_word++,(*render)->filename_halo_type);
          if(i_word>count_words(line))
             sprintf((*render)->filename_tree_version,"version_nominal_res");
          else
             grab_word(line,i_word++,(*render)->filename_tree_version);
        }
        else if(!strcmp(parameter,"snap_file"))
          grab_word(line,i_word++,(*render)->snap_filename_root);
        else if(!strcmp(parameter,"mark_file")){
          grab_word(line,i_word++,(*render)->mark_filename_root);
          (*render)->flag_read_marked=TRUE;
        }
        else if(!strcmp(parameter,"smooth_file"))
          grab_word(line,i_word++,(*render)->smooth_filename_root);
        else if(!strcmp(parameter,"snap_a_list_file")){
          SID_log("Reading snapshot a_list...",SID_LOG_OPEN);
          grab_word(line,i_word++,(*render)->snap_a_list_filename);
          if((fp_list=fopen((*render)->snap_a_list_filename,"r"))==NULL)
            SID_trap_error("Could not open snapshot a_list {%s}",ERROR_IO_OPEN,(*render)->snap_a_list_filename);
          (*render)->n_snap_a_list=count_lines_data(fp_list);
          (*render)->snap_a_list  =(double *)SID_malloc(sizeof(double)*(*render)->n_snap_a_list);
          for(i=0;i<(*render)->n_snap_a_list;i++){
            grab_next_line_data(fp_list,&line,&line_length);
            grab_double(line,1,&((*render)->snap_a_list[i]));
          }
          fclose(fp_list);
          for(i=0;i<(*render)->n_snap_a_list;i++)
            SID_log("a[%3d]=%lf",SID_LOG_COMMENT,i,(*render)->snap_a_list[i]);
          SID_log("Done.",SID_LOG_CLOSE);
        }
        else if(!strcmp(parameter,"h_Hubble"))
          grab_double(line,i_word++,&((*render)->h_Hubble));
        else if(!strcmp(parameter,"f_absorption")){
          if((*render)->flag_add_absorption)
             SID_trap_error("There are conflicting absorption criteria.",ERROR_LOGIC);
          grab_double(line,i_word++,&((*render)->f_absorption));
          if((*render)->f_absorption<0.)
             SID_trap_error("f_absorption has been set to %le but must be >=0.",ERROR_LOGIC,(*render)->f_absorption);
          if((*render)->f_absorption>0.)
             (*render)->flag_add_absorption=TRUE;
        }
        else if(!strcmp(parameter,"snap_number"))
          grab_int(line,i_word++,&((*render)->snap_number));
        else if(!strcmp(parameter,"flag_comoving"))
          grab_int(line,i_word++,&((*render)->flag_comoving));
        else if(!strcmp(parameter,"flag_fade"))
          (*render)->flag_fade=TRUE;
        else if(!strcmp(parameter,"force_periodic"))
          (*render)->flag_force_periodic=TRUE;
        else if(!strcmp(parameter,"camera")){
          grab_word(line,i_word++,variable);
          if(!strcmp(variable,"size")){
            grab_int(line,i_word++,&i_value);
            (*render)->camera->width =i_value;
            grab_int(line,i_word++,&i_value);
            (*render)->camera->height=i_value;
          }
          else if(!strcmp(variable,"stereo_factor")){
            grab_double(line,i_word++,&d_value);
            (*render)->camera->stereo_ratio=d_value;
            (*render)->camera->camera_mode|=CAMERA_STEREO;
          }
          else if(!strcmp(variable,"colour_table")){
            grab_int(line,i_word++,&i_value);
            (*render)->camera->colour_table=i_value;
          }
          else if(!strcmp(variable,"plane_parallel")){
            (*render)->camera->camera_mode|=CAMERA_PLANE_PARALLEL|CAMERA_MONO;
          }
          else if(!strcmp(variable,"RGB_param")){
            grab_word(line,i_word++,c_value);
            strcpy((*render)->camera->RGB_param,c_value);
          }
          else if(!strcmp(variable,"RGB_range")){
            grab_double(line,i_word++,&d_value);
            (*render)->camera->RGB_range[0]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->camera->RGB_range[1]=d_value;
            while(count_words(line)>=i_word){
               char next_word[32];
               grab_word(line,i_word,next_word);
               if(!strcmp(next_word,"linear"))
                  (*render)->v_mode&=(~MAKE_MAP_LOG);
               else if(!strcmp(next_word,"log"))
                  (*render)->v_mode|=MAKE_MAP_LOG;
               else
                  SID_log("Ignoring unknown RGB_range modifier {%s}.",SID_LOG_COMMENT,next_word);
               i_word++;
            }
          }
          else if(!strcmp(variable,"RGB_gamma"))
            i_word=set_transfer_function(line,i_word,&((*render)->camera->RGB_gamma));
          else if(!strcmp(variable,"transfer")){
            grab_word(line,i_word++,parameter);
            grab_word(line,i_word++,temp_word);
            if(!strcmp(temp_word,"log") || !strcmp(temp_word,"LOG"))
              flag=TRUE;
            else if(!strcmp(temp_word,"linear") || !strcmp(temp_word,"LINEAR"))
              flag=FALSE;
            else
              SID_trap_error("log/linear flag not set to 'log' or 'linear' {%s}",ERROR_LOGIC,temp_word);
            i_word=set_transfer_function(line,i_word,&temp_interp);
            // Add to the transfer function list
            ADaPS_store(&((*render)->camera->transfer_list),(void *)temp_interp,parameter,ADaPS_DEFAULT);
            // If this is a log-defined transfer function, store that fact
            if(flag)
              ADaPS_store(&((*render)->camera->transfer_list),(void *)&flag,"%s_log",ADaPS_SCALAR_INT,parameter);
          }
          else if(!strcmp(variable,"Y_param")){
            grab_word(line,i_word++,c_value);
            strcpy((*render)->camera->Y_param,c_value);
          }
          else if(!strcmp(variable,"Y_range")){
            grab_double(line,i_word++,&d_value);
            (*render)->camera->Y_range[0]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->camera->Y_range[1]=d_value;
            while(count_words(line)>=i_word){
               char next_word[32];
               grab_word(line,i_word,next_word);
               if(!strcmp(next_word,"linear"))
                  (*render)->w_mode&=(~MAKE_MAP_LOG);
               else if(!strcmp(next_word,"log"))
                  (*render)->w_mode|=MAKE_MAP_LOG;
               else
                  SID_log("Ignoring unknown RGB_range modifier {%s}.",SID_LOG_COMMENT,next_word);
               i_word++;
            }
          }
          else if(!strcmp(variable,"Y_gamma")){
            n_transfer=count_words(line)-i_word+1;
            if(n_transfer>2){
              transfer_array_x=(double *)SID_malloc(sizeof(double)*n_transfer);
              transfer_array_y=(double *)SID_malloc(sizeof(double)*n_transfer);
              if((*render)->camera->Y_gamma!=NULL)
                free_interpolate(SID_FARG (*render)->camera->Y_gamma,NULL);
              for(j_word=0;j_word<n_transfer;j_word++){
                grab_double(line,i_word++,&d_value);
                transfer_array_x[j_word]=255.*(double)j_word/(double)(n_transfer-1);
                transfer_array_y[j_word]=d_value;
              }
              init_interpolate(transfer_array_x,transfer_array_y,n_transfer,gsl_interp_cspline,&((*render)->camera->Y_gamma));
              SID_free(SID_FARG transfer_array_x);
              SID_free(SID_FARG transfer_array_y);
            }
            else
              SID_log_warning("Gamma arrays bust be >2 elements long.",ERROR_LOGIC);
          }
          else if(!strcmp(variable,"Z_range")){
            grab_double(line,i_word++,&d_value);
            (*render)->camera->Z_range[0]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->camera->Z_range[1]=d_value;
          }
          else if(!strcmp(variable,"Z_gamma")){
            n_transfer=count_words(line)-i_word+1;
            if(n_transfer>2){
              transfer_array_x=(double *)SID_malloc(sizeof(double)*n_transfer);
              transfer_array_y=(double *)SID_malloc(sizeof(double)*n_transfer);
              if((*render)->camera->Z_gamma!=NULL)
                free_interpolate(SID_FARG (*render)->camera->Z_gamma,NULL);
              for(j_word=0;j_word<n_transfer;j_word++){
                grab_double(line,i_word++,&d_value);
                transfer_array_x[j_word]=255.*(double)j_word/(double)(n_transfer-1);
                transfer_array_y[j_word]=d_value;
              }
              init_interpolate(transfer_array_x,transfer_array_y,n_transfer,gsl_interp_cspline,&((*render)->camera->Z_gamma));
              SID_free(SID_FARG transfer_array_x);
              SID_free(SID_FARG transfer_array_y);
            }
            else
              SID_log_warning("Gamma arrays bust be >2 elements long.",ERROR_LOGIC);
          }
          else if(!strcmp(variable,"f_near_field"))
            grab_double(line,i_word++,&((*render)->camera->f_near_field));
          else if(!strcmp(variable,"f_taper_field"))
            grab_double(line,i_word++,&((*render)->camera->f_taper_field));
          else if(!strcmp(variable,"f_image_plane")){
            grab_double(line,i_word++,&((*render)->camera->f_image_plane));
            if((*render)->camera->f_image_plane<=0.)
               SID_trap_error("f_image_plane (%le) must be >0.  Consider using a plane-parallel projection.",ERROR_LOGIC,(*render)->camera->f_image_plane);
          }
          else
            SID_trap_error("Unknown variable {%s} for parameter {%s} for command {%s} on line %d",ERROR_LOGIC,variable,parameter,command,i_line);
        }
        else if(!strcmp(parameter,"evolve")){
          grab_word(line,i_word++,variable);
          if(!strcmp(variable,"FOV")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->FOV=d_value;
          }
          else if(!strcmp(variable,"p_o")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->p_o[0]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->p_o[1]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->p_o[2]=d_value;
          }
          else if(!strcmp(variable,"radius")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->radius=d_value;
          }
          else if(!strcmp(variable,"theta")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->theta=d_value;
          }
          else if(!strcmp(variable,"phi")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->phi=d_value;
          }
          else if(!strcmp(variable,"zeta")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->zeta=d_value;
          }
          else if(!strcmp(variable,"time")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->time=d_value;
          }
          else
            SID_trap_error("Unknown variable {%s} for parameter {%s} for command {%s} on line %d",ERROR_LOGIC,variable,parameter,command,i_line);
        }      
        else
          SID_trap_error("Unknown parameter {%s} for command {%s} on line %d",ERROR_LOGIC,parameter,command,i_line);
      }
      else if(!strcmp(command,"add")){
        grab_word(line,i_word++,parameter);
        if(!strcmp(parameter,"scene"))
          add_render_scene((*render));
        else if(!strcmp(parameter,"perspective"))
          add_scene_perspective((*render)->last_scene);
        else
          SID_trap_error("Unknown parameter {%s} for command {%s} on line %d",ERROR_LOGIC,parameter,command,i_line);
      }
      else
        SID_trap_error("Unknown render command {%s} on line %d",ERROR_LOGIC,command,i_line);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
  fclose(fp);
  if(!(*render)->sealed)
    seal_render(*render);
  if(line!=NULL)
     SID_free(SID_FARG line);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  SID_log("Done.",SID_LOG_CLOSE);
}

