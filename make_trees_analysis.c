#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  char filename_run_root[MAX_FILENAME_LENGTH];
  char filename_out_root[MAX_FILENAME_LENGTH];
  strcpy(filename_trees_root,argv[1]);
  strcpy(filename_halos_root,argv[2]);
  strcpy(filename_run_root,  argv[3]);
  strcpy(filename_out_root,  argv[4]);

  SID_log("Performing analysis of merger trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Perform analysis
  tree_info *trees;
  read_trees(filename_trees_root,
             filename_halos_root,
             filename_run_root,
             TREE_MODE_DEFAULT,
             &trees);

  // Compute merger rates ...
  generate_trees_analysis(trees,filename_out_root);

  // Clean-up 
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

