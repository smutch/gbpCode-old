#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>
#include <assert.h>

void init_treenode_hist(tree_info           *trees,
                        const char          *hist_name,
                        const char          *x_name_in,
                        const char          *y_name_in,
                        int                  mode,
                        treenode_hist_info **hist, ...){
  va_list vargs;
  va_start(vargs,hist);

  (*hist)=(treenode_hist_info *)SID_malloc(sizeof(treenode_hist_info));

  // Set histogram name
  sprintf((*hist)->name,"%s",hist_name);

  // Parse the mode
  if(check_mode_for_flag(mode,TREENODE_HIST_LOG_X))
    (*hist)->flag_log_x=TRUE;
  else
    (*hist)->flag_log_x=FALSE;
  if(check_mode_for_flag(mode,TREENODE_HIST_LOG_Y))
    (*hist)->flag_log_y=TRUE;
  else
    (*hist)->flag_log_y=FALSE;
  
  // Parse the axis information
  int flag_x_found=FALSE;
  int flag_y_found=FALSE;
  int n_args_x    =0;
  int n_args_y    =0;
  int i_prop_x    =0;
  int i_prop_y    =0;
  int i_prop      =0;
  for(i_prop=0;i_prop<TREENODE_HIST_N_PROPS;i_prop++){
     if(!strcmp(x_name_in,treenode_hist_props.name[i_prop])){
        i_prop_x       =i_prop;
        (*hist)->x_prop=i_prop_x;
        n_args_x       =treenode_hist_props.n_args[i_prop_x];
        flag_x_found   =TRUE;
        if(n_args_x>TREENODE_HIST_N_ARGS_MAX)
           SID_trap_error("Increase TREENODE_HIST_N_ARGS_MAX to at least %d.",ERROR_LOGIC,n_args_x);
     }
     if(!strcmp(y_name_in,treenode_hist_props.name[i_prop])){
        i_prop_y       =i_prop;
        (*hist)->y_prop=i_prop_y;
        n_args_y       =treenode_hist_props.n_args[i_prop_y];
        flag_y_found   =TRUE;
        if(n_args_y>TREENODE_HIST_N_ARGS_MAX)
           SID_trap_error("Increase TREENODE_HIST_N_ARGS_MAX to at least %d.",ERROR_LOGIC,n_args_y);
     }
  }

  //Sanity check
  if(!flag_x_found)
     SID_trap_error("Could not find property definition for {%s} in init_treenode_hist.",ERROR_LOGIC,x_name_in);
  if(!flag_y_found)
     SID_trap_error("Could not find property definition for {%s} in init_treenode_hist.",ERROR_LOGIC,y_name_in);

  // Parse the variable arguments depending on what we've specified for axes
  for(int i_arg=0;i_arg<n_args_x;i_arg++){
     if(treenode_hist_props.arg_type[i_prop_x][i_arg]==SID_DOUBLE)
        (*hist)->x_args_d[i_arg]=(double)va_arg(vargs,double);
     else if(treenode_hist_props.arg_type[i_prop_x][i_arg]==SID_INT)
        (*hist)->x_args_i[i_arg]=(int)va_arg(vargs,int);
     else
        SID_trap_error("Invalid datatype specified in init_treenode_hist for x-axis.",ERROR_LOGIC);
  }
  for(int i_arg=0;i_arg<n_args_y;i_arg++){
     if(treenode_hist_props.arg_type[i_prop_y][i_arg]==SID_DOUBLE)
        (*hist)->y_args_d[i_arg]=(double)va_arg(vargs,double);
     else if(treenode_hist_props.arg_type[i_prop_y][i_arg]==SID_INT)
        (*hist)->y_args_i[i_arg]=(int)va_arg(vargs,int);
     else
        SID_trap_error("Invalid datatype specified in init_treenode_hist for y-axis.",ERROR_LOGIC);
  }

  // Specify the size of the histogram in each dimension
  for(int i_axis=0;i_axis<2;i_axis++){
     int    *args_i;
     double *args_d;
     int    *n_d;
     switch(i_axis){
        case 0:
          i_prop=i_prop_x;
          args_i=(*hist)->x_args_i;
          args_d=(*hist)->x_args_d;
          n_d   =&((*hist)->n_x);
          break;
        case 1:
          i_prop=i_prop_y;
          args_i=(*hist)->y_args_i;
          args_d=(*hist)->y_args_d;
          n_d   =&((*hist)->n_y);
          break;
     }
     switch(i_prop){
       case 0:{
          int n_z_bin=args_i[0];
          (*n_d)=trees->n_snaps/n_z_bin;
          break;
       }
       case 1:
       case 2:{
          (*n_d)=args_i[2];
          break;
       }
     }
  }

  // Allocate the array
  (*hist)->array=(int *)SID_calloc(sizeof(int)*(*hist)->n_x*(*hist)->n_y);

  va_end(vargs);
}

