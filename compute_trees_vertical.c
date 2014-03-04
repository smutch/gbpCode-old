#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void compute_trees_vertical(char       *filename_root_out,
                            int         n_dim_files){
  SID_log("Constructing vertical merger trees for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,
          trees->i_read_start,
          trees->i_read_stop,
          trees->i_read_step);

  // Read the horizontal trees into RAM
  tree_info *trees;
  read_trees(filename_trees_root,
             filename_halos_root,
             TREE_MODE_DEFAULT,
             &trees);

  // Read ancillary data
  read_trees_catalogs(trees,
                      filename_SSimPL_dir,
                      filename_halo_version_root,
                      READ_TREES_CATALOGS_GROUPS|READ_TREES_CATALOGS_SUBGROUPS);

  // Write vertical trees

  // Clean-up
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
}

