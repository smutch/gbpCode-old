#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>
#include <assert.h>

void write_treenode_hist(tree_info *trees,const char *filename_out_root,treenode_hist_info *hist){

  SID_log("Writing treenode histograms to",SID_LOG_OPEN);

  // Set filename roots
  char filename_out[MAX_FILENAME_LENGTH];
  char x_name[TREENODE_HIST_NAME_LENGTH];
  char y_name[TREENODE_HIST_NAME_LENGTH];
  sprintf(x_name,"%s",treenode_hist_props.name[hist->x_prop]);
  sprintf(y_name,"%s",treenode_hist_props.name[hist->y_prop]);
  sprintf(filename_out,"%s_%s_hist_%s_%s.txt",filename_out_root,hist->name,x_name,y_name);

  SID_log(" {%s}...",SID_LOG_CONTINUE,filename_out);

  // Open file
  FILE *fp_out=fopen(filename_out,"w");

  // Write axis information
  int n_x=0;
  int n_y=0;
  fprintf(fp_out,"%s\n",hist->name);
  for(int i_axis=0;i_axis<2;i_axis++){
     int     n_bins;
     int     i_prop;
     char   *axis_name;
     int    *args_i;
     double *args_d;
     int     flag_log_axis;
     if(i_axis==0){
        n_bins=hist->n_x;
        i_prop=hist->x_prop;
        args_i=hist->x_args_i;
        args_d=hist->x_args_d;
        n_x   =n_bins;
        flag_log_axis=hist->flag_log_x;
     }
     else{
        n_bins=hist->n_y;
        i_prop=hist->y_prop;
        args_i=hist->y_args_i;
        args_d=hist->y_args_d;
        n_y   =n_bins;
        flag_log_axis=hist->flag_log_y;
     }
     if(flag_log_axis)
        axis_name=treenode_hist_props.log_axis_text[i_prop];
     else
        axis_name=treenode_hist_props.axis_text[i_prop];
     fprintf(fp_out,"%d %s\n",n_bins,axis_name);
     for(int i_bin=0;i_bin<=n_bins;i_bin++){
        double bin_i;
        switch(i_prop){
           case 0:{ // z
             int di   =args_i[0];
             int i_min=trees->n_snaps%di;
             if(i_bin==0)
                bin_i=trees->z_list[i_min];
             else
                bin_i=trees->z_list[MIN(i_min+i_bin*di-1,trees->n_snaps-1)];
             break;
           }
           case 1: // M
           case 2: // M
           case 3:{ // N
             double d_min=args_d[0];
             double dd   =args_d[1];
             bin_i       =d_min+((double)i_bin)*dd;
             break;
           }
           default:
             SID_trap_error("Invalid property passed to write_treenode_hist().",ERROR_LOGIC);
             break;
        }
        fprintf(fp_out,"%le\n",bin_i);
     }
  }

  for(int i_y=0;i_y<n_y;i_y++){
     for(int i_x=0;i_x<n_x;i_x++) fprintf(fp_out,"%le ",(double)(hist->array[i_y*n_x+i_x]));fprintf(fp_out,"\n");
  }
  fclose(fp_out);

  SID_log("Done.",SID_LOG_CLOSE);
}

