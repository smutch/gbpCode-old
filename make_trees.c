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
  int         i_read_start;
  int         i_read_stop;
  int         i_read_step;
  int         n_search;
  int         n_files_groups;
  int         n_files_subgroups;
  int         n_k_match=2;
  int         flag_clean=FALSE;

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  strcpy(filename_halo_root_in,argv[1]);
  strcpy(filename_cat_root_in, argv[2]);
  strcpy(filename_root_out,    argv[3]);
  i_read_start     =atoi(argv[4]);
  i_read_stop      =atoi(argv[5]);
  i_read_step      =atoi(argv[6]);
  n_search         =atoi(argv[7]);
  n_files_groups   =atoi(argv[8]);
  n_files_subgroups=atoi(argv[9]);

  SID_log("Constructing merger trees for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_start,i_read_stop,i_read_step);
  compute_trees_horizontal(filename_halo_root_in,
                           filename_cat_root_in,
                           filename_root_out,
                           i_read_start,
                           i_read_stop,
                           i_read_step,
                           n_search,
                           &flag_clean);

  //flag_clean=FALSE;
  compute_trees_vertical(filename_root_out,
                         i_read_start,
                         i_read_stop,
                         i_read_step,
                         n_search,
                         n_files_groups,
                         n_files_subgroups,
                         &flag_clean);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(0);
}
