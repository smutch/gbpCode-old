#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){
  char        filename_halo_root_in[256];
  char        filename_cat_root_in[256];
  char        filename_root_out[256];
  char        filename_output_dir_horizontal_trees[256];
  char        filename_output_dir[256];
  char        filename_output_dir_horizontal[256];
  char        filename_root_matches[256];
  char        filename_snap_list_in[256];
  char        filename_snap_list_out[256];
  char        filename_output_file_root[256];
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
  cosmo_info *cosmo;
  int         flag_fix_bridges;

  SID_init(&argc,&argv,NULL);

  // Initialize cosmology
  init_cosmo_std(&cosmo);

  // Fetch user inputs
  if(argc!=3)
    SID_trap_error("Incorrect syntax",ERROR_SYNTAX);
  strcpy(filename_root_out,    argv[1]);
  strcpy(filename_snap_list_in,argv[2]);

  SID_log("Constructing snapshot list for trees in {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_root_out);

  int flag_compute_fragmented;
  int flag_compute_ghosts;
  SID_log("Fetching run parameters...",SID_LOG_OPEN);
  read_tree_run_parameters(filename_root_out,
                           &i_read_start,
                           &i_read_stop,
                           &i_read_step,
                           &n_search,
                           &flag_fix_bridges,
                           &flag_compute_fragmented,
                           &flag_compute_ghosts);
  SID_log("i_read_start=%d",SID_LOG_COMMENT,i_read_start);
  SID_log("i_read_stop =%d",SID_LOG_COMMENT,i_read_stop);
  SID_log("i_read_step =%d",SID_LOG_COMMENT,i_read_step);
  SID_log("Done.",SID_LOG_CLOSE);
  write_a_list(filename_snap_list_in,
               filename_root_out,
               i_read_start,
               i_read_stop,
               i_read_step);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
