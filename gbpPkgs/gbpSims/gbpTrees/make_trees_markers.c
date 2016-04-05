#define _MAIN
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  char filename_SSimPL_dir[MAX_FILENAME_LENGTH];
  char filename_halo_version_root[MAX_FILENAME_LENGTH];
  char filename_trees_name[MAX_FILENAME_LENGTH];
  char filename_output_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_dir,       argv[1]);
  strcpy(filename_halo_version_root,argv[2]);
  strcpy(filename_trees_name,       argv[3]);
  strcpy(filename_output_root,      argv[4]);

  // Set the halo and tree filename roots
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_dir,filename_trees_name);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_dir,filename_halo_version_root);

  SID_log("Generating treenode markers & analysis of merger trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Perform analysis
  tree_info *trees;
  read_trees(filename_SSimPL_dir,
             filename_halo_version_root,
             filename_trees_name,
             TREE_MODE_DEFAULT,
             &trees);

  // Read catalogs
  read_trees_catalogs(trees,READ_TREES_CATALOGS_BOTH);

  // Loop over the two halo types
  for(int i_type=0;i_type<2;i_type++){
     // Initialize markers (just one-at-a-time to save RAM)
     int mode;
     if(i_type==0)
        mode=PRECOMPUTE_TREENODE_MARKER_GROUPS;
     else
        mode=PRECOMPUTE_TREENODE_MARKER_SUBGROUPS;
     // Compute markers
     precompute_treenode_markers(trees,mode);
     // Write markers
     write_treenode_markers(trees,filename_output_root,mode);
     // Free markers
     free_precompute_treenode_markers(trees,mode);
  }

  // Clean-up
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

