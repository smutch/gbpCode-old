#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>
#include <gbpHighZ.h>

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  char filename_SSimPL_dir[MAX_FILENAME_LENGTH];
  char filename_halo_version_root[MAX_FILENAME_LENGTH];
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_trees_name[MAX_FILENAME_LENGTH];
  char filename_output_root[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_dir,       argv[1]);
  strcpy(filename_halo_version_root,argv[2]);
  strcpy(filename_trees_name,       argv[3]);
  strcpy(filename_output_root,      argv[4]);

  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_dir,filename_trees_name);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_dir,filename_halo_version_root);

  SID_log("Generating treenode markers & analysis of merger trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Perform analysis
  tree_info *trees;
  read_trees(filename_trees_root,
             filename_halos_root,
             TREE_MODE_DEFAULT,
             &trees);

  // Read catalogs
  read_trees_catalogs(trees,
                      filename_SSimPL_dir,
                      filename_halo_version_root,
                      READ_TREES_CATALOGS_BOTH);

  for(int i_type=0;i_type<2;i_type++){
     // Set some stuff which varies depending on whether we are
     //    analyzing groups or subgroups. 
     int mode;
     switch(i_type){
        case 0:
           mode=COMPUTE_ACCRETION_ANALYSIS_GROUPS;
           SID_log("Analyzing groups...",SID_LOG_OPEN|SID_LOG_TIMER);
           break;
        case 1:
           mode=COMPUTE_ACCRETION_ANALYSIS_SUBGROUPS;
           SID_log("Analyzing subgroups...",SID_LOG_OPEN|SID_LOG_TIMER);
           break;
     }
   
     // Generate markers
     tree_markers_info **markers=NULL;
     find_treenode_markers_all(trees,&markers,mode);
   
     // Write markers
     write_treenode_markers_all(trees,markers,filename_output_root,mode);
   
     // Compute merger tree statistics
     compute_marker_analysis(trees,markers,filename_trees_name,filename_output_root,mode);

     // Clean-up
     free_treenode_markers_all(trees,&markers);

     SID_log("Done.",SID_LOG_CLOSE);
  }

  // Clean-up
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}


