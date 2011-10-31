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
  char        filename_root_matches[256];
  char        filename_root_out[256];
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
  int         n_files;
  int        *n_subgroups;
  int        *n_groups;
  int         flag_fix_bridges=TRUE;

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  strcpy(filename_halo_root_in,argv[1]);
  strcpy(filename_root_matches,argv[2]);
  i_read_start           =atoi(argv[3]);
  i_read_stop            =atoi(argv[4]);
  n_search               =atoi(argv[5]);
  i_read_step=1;

  SID_log("Constructing merger tree matches for snapshots #%d->#%d (n_search=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_start,i_read_stop,n_search);

  // Perform matching
  compute_trees_matches(filename_halo_root_in,
                        filename_root_matches,
                        i_read_start,
                        i_read_stop,
                        i_read_step,
                        &n_files,
                        &n_subgroups,
                        &n_groups,
                        n_search);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
