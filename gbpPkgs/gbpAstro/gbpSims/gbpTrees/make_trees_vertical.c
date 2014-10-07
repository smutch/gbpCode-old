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

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  char filename_SSimPL_dir[MAX_FILENAME_LENGTH];
  char filename_halo_version_root[MAX_FILENAME_LENGTH];
  char filename_trees_name[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_dir,       argv[1]);
  strcpy(filename_halo_version_root,argv[2]);
  strcpy(filename_trees_name,       argv[3]);
  double box_size     =(double)atof(argv[4]);
  int    n_dim_files  =        atoi(argv[5]);

  SID_log("Constructing vertical merger trees...",SID_LOG_OPEN|SID_LOG_TIMER);
  compute_trees_vertical(filename_SSimPL_dir,
                         filename_halo_version_root,
                         filename_trees_name,
                         box_size,
                         n_dim_files);

  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}

