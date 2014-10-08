#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void compute_trees_vertical(char   *filename_SSimPL_dir,
                            char   *filename_halo_version_root,
                            char   *filename_trees_name,
                            double  box_size,
                            int     n_dim_files){
  SID_log("Constructing vertical merger trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_dir,filename_trees_name);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_dir,filename_halo_version_root);

  // Read the horizontal trees into RAM
  tree_info *trees;
  read_trees(filename_SSimPL_dir,
             filename_halo_version_root,
             filename_trees_name,
             TREE_MODE_DEFAULT,
             &trees);

  // Read ancillary data
  read_trees_catalogs(trees,
                      filename_SSimPL_dir,
                      filename_halo_version_root,
                      READ_TREES_CATALOGS_BOTH|READ_TREES_CATALOGS_SAGE);

  // Build depth-first-index pointers
  finalize_trees_vertical(trees);

  // Write vertical trees
  write_trees_vertical(trees,
                       box_size,
                       n_dim_files,
                       filename_trees_root);

  // Clean-up
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
}

