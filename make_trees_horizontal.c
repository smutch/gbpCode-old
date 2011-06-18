#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){
  char        filename_halo_root_in[256];
  char        filename_cat_root_in[256];
  char        filename_root_out[256];
  char        filename_snap_list_in[256];
  char        filename_snap_list_out[256];
  int         i_read_start;
  int         i_read_stop;
  int         i_read_step;
  int         n_search;
  int         n_files_groups;
  int         n_files_subgroups;
  int         n_k_match=2;
  int         flag_clean=FALSE;
  FILE       *fp_in;
  FILE       *fp_out;
  char       *line=NULL;
  size_t      line_length=0;
  int         n_snaps,i_read,i_next,i_write,n_keep;
  double     *a_list_in;
  double     *a_list_out;

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  strcpy(filename_halo_root_in,argv[1]);
  strcpy(filename_cat_root_in, argv[2]);
  strcpy(filename_root_out,    argv[3]);
  strcpy(filename_snap_list_in,argv[4]);
  i_read_start     =atoi(argv[5]);
  i_read_stop      =atoi(argv[6]);
  i_read_step      =atoi(argv[7]);
  n_search         =atoi(argv[8]);
  n_files_groups   =atoi(argv[9]);
  n_files_subgroups=atoi(argv[10]);

  // Create snapshot expansion factor list
  if(SID.I_am_Master){
    sprintf(filename_snap_list_out,"%s.a_list",filename_root_out);
    SID_log("Creating snapshot list {%s}...",SID_LOG_OPEN,filename_snap_list_out);

    // Read the original list
    fp_in     =fopen(filename_snap_list_in, "r");
    n_snaps   =count_lines_data(fp_in);
    a_list_in =(double *)SID_malloc(sizeof(double)*n_snaps);
    a_list_out=(double *)SID_malloc(sizeof(double)*n_snaps);
    for(i_read=0;i_read<n_snaps;i_read++){
      grab_next_line_data(fp_in,&line,&line_length);
      grab_double(line,1,&(a_list_in[i_read]));
    }
    fclose(fp_in);

    // Select the snapshots we've used
    for(i_read=i_read_stop,i_next=i_read_stop,n_keep=0;i_read>=i_read_start;i_read--){
      if(i_read==i_next){
        a_list_out[n_keep++]=a_list_in[i_read];
        i_next-=i_read_step;
      }
    }

    // Write them to a file
    fp_out=fopen(filename_snap_list_out,"w");
    for(i_write=n_keep-1;i_write>=0;i_write--){
      fprintf(fp_out,"%le\n",a_list_out[i_write]);
    }
    fclose(fp_out);

    SID_free(SID_FARG a_list_in);
    SID_log("Done.",SID_LOG_CLOSE);
  }

  SID_log("Constructing merger trees for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_start,i_read_stop,i_read_step);
  compute_trees_horizontal(filename_halo_root_in,
                           filename_cat_root_in,
                           filename_root_out,
                           a_list_out,
                           i_read_start,
                           i_read_stop,
                           i_read_step,
                           n_search,
                           &flag_clean);
  SID_log("Done.",SID_LOG_CLOSE);

  if(SID.I_am_Master)
     SID_free(SID_FARG a_list_out);
  
  SID_exit(ERROR_NONE);
}
