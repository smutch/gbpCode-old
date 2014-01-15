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
  int         i_read_1;
  int         i_read_2;
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
  i_read_1               =atoi(argv[3]);
  i_read_2               =atoi(argv[4]);

  SID_log("Constructing matches between snapshots #%d and #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_1,i_read_2);

  // Perform matching
  compute_single_matches(filename_halo_root_in,
                         filename_root_matches,
                         i_read_1,
                         i_read_2);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

