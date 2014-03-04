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
  char filename_SSimPL_dir[MAX_FILENAME_LENGTH];
  char filename_halo_version_root[MAX_FILENAME_LENGTH];
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_trees_name[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  char filename_out_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_dir,       argv[1]);
  strcpy(filename_halo_version_root,argv[2]);
  strcpy(filename_trees_name,       argv[3]);

  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_dir,filename_trees_name);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_dir,filename_halo_version_root);
  sprintf(filename_out_root,  "%s",filename_trees_name);

  SID_log("Performing analysis of merger trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Read the trees
  tree_info *trees;
  read_trees(filename_trees_root,
             filename_halos_root,
             TREE_MODE_DEFAULT,
             &trees);

  // Read ancillary data
  //read_trees_catalogs(trees,
  //                    filename_SSimPL_dir,
  //                    filename_halo_version_root,
  //                    READ_TREES_CATALOGS_GROUPS|READ_TREES_CATALOGS_SUBGROUPS);
  read_trees_match_scores(trees,
                          filename_SSimPL_dir,
                          READ_TREES_MATCH_SCORES_ALL);

  // Compute merger tree statistics
  compute_trees_analysis(trees,filename_out_root);

  // Clean-up
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

