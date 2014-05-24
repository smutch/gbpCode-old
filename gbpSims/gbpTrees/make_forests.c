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
  char        filename_root_out[256];

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  strcpy(filename_root_out,argv[1]);

  // Read the run parameters
  int i_read_start;
  int i_read_stop;
  int i_read_step;
  int n_search;
  int flag_fix_bridges;
  int flag_compute_fragmented;
  int flag_compute_ghosts;
  read_trees_run_parameters(filename_root_out,
                            &i_read_start,
                            &i_read_stop,
                            &i_read_step,
                            &n_search,
                            &flag_fix_bridges,
                            &flag_compute_fragmented,
                            &flag_compute_ghosts);

  // Default to scanning the whole range of snapshots when building forests
  int n_search_forests=i_read_stop;

  // Generate mapping
  compute_forests(filename_root_out,n_search_forests);

  SID_exit(ERROR_NONE);
}

