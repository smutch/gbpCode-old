#define _MAIN
#include <stdio.h>
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
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_trees_reference_root[MAX_FILENAME_LENGTH];
  char filename_trees_name[MAX_FILENAME_LENGTH];
  char filename_trees_reference_name[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  int  i_arg=1;
  strcpy(filename_SSimPL_dir,       argv[i_arg++]);
  strcpy(filename_halo_version_root,argv[i_arg++]);
  strcpy(filename_trees_name,       argv[i_arg++]);
  if(argc==8)
     strcpy(filename_trees_reference_name,argv[i_arg++]);
  else
     sprintf(filename_trees_reference_name,"");
  double logM_min=atof(argv[i_arg++]);
  double dlogM   =atof(argv[i_arg++]);
  int    n_logM  =atoi(argv[i_arg++]);

  // Set some filenames
  sprintf(filename_trees_root,          "%s/trees/%s",filename_SSimPL_dir,filename_trees_name);
  sprintf(filename_trees_reference_root,"%s/trees/%s",filename_SSimPL_dir,filename_trees_reference_name);
  sprintf(filename_halos_root,          "%s/halos/%s",filename_SSimPL_dir,filename_halo_version_root);

  SID_log("Performing analysis of merger trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Read trees
  tree_info *trees;
  read_trees(filename_SSimPL_dir,
             filename_halo_version_root,
             filename_trees_name,
             TREE_MODE_DEFAULT|TREE_READ_EXTENDED_POINTERS,
             &trees);

  // Read reference trees
  if(strcmp(filename_trees_reference_name,""))
     read_trees(filename_SSimPL_dir,
                filename_halo_version_root,
                filename_trees_reference_name,
                TREE_MODE_DEFAULT|TREE_READ_EXTENDED_POINTERS|TREE_MODE_REFERENCE,
                &(trees->trees_reference));

  // Read catalogs.  Use SHORT read to save RAM.
  read_trees_catalogs(trees,READ_TREES_CATALOGS_BOTH|READ_TREES_CATALOGS_SHORT);
  if(trees->trees_reference!=NULL)
     read_trees_catalogs(trees->trees_reference,READ_TREES_CATALOGS_BOTH|READ_TREES_CATALOGS_SHORT);

  //read_trees_match_scores(trees,
  //                        filename_SSimPL_dir,
  //                        READ_TREES_MATCH_SCORES_ALL);

  // Perform analysis
  compute_trees_analysis(trees,logM_min,dlogM,n_logM);

  // Clean-up
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

