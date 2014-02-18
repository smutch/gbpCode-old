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

void init_trees(int         i_read_start,
                int         i_read_stop,
                int         i_read_step,
                int         n_forests,
                int         n_forests_local,
                tree_info **tree){

  SID_log("Initializing trees...",SID_LOG_OPEN);

  // Count the number of snaps being used
  int n_snaps=0;
  int i_read;
  int i_file;
  for(i_read=i_read_stop;i_read>=i_read_start;i_read-=i_read_step)
     n_snaps++;

  // Allocate the data structure
  (*tree)=(tree_info *)SID_malloc(sizeof(tree_info));

  // Set counts etc
  (*tree)->n_snaps        =n_snaps;
  (*tree)->i_read_start   =i_read_start;
  (*tree)->i_read_stop    =i_read_stop;
  (*tree)->i_read_step    =i_read_step;
  (*tree)->n_forests      =n_forests;
  (*tree)->n_forests_local=n_forests_local;

  // Create an array which maps the file numbers in the trees
  //   to the snapshot number (may differ from 1:1 if skipping snaps)
  (*tree)->snap_list=(int    *)SID_malloc(sizeof(int)   *n_snaps);
  (*tree)->z_list   =(double *)SID_malloc(sizeof(double)*n_snaps);
  (*tree)->t_list   =(double *)SID_malloc(sizeof(double)*n_snaps);
  for(i_read=i_read_stop,i_file=n_snaps-1;i_read>=i_read_start;i_read-=i_read_step,i_file--){
     (*tree)->snap_list[i_file]=i_read;
     (*tree)->z_list[i_file]   =0.;
     (*tree)->t_list[i_file]   =0.;
  }

  // Allocate counters
  (*tree)->n_groups_snap_local   =(int *)SID_calloc(sizeof(int)*n_snaps);
  (*tree)->n_subgroups_snap_local=(int *)SID_calloc(sizeof(int)*n_snaps);

  // Allocate pointers
  (*tree)->first_neighbour_groups   =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_snaps);
  (*tree)->first_neighbour_subgroups=(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_snaps);
  (*tree)->last_neighbour_groups    =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_snaps);
  (*tree)->last_neighbour_subgroups =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_snaps);
  (*tree)->first_in_forest_groups   =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_forests_local);
  (*tree)->first_in_forest_subgroups=(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_forests_local);
  (*tree)->last_in_forest_groups    =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_forests_local);
  (*tree)->last_in_forest_subgroups =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_forests_local);
  (*tree)->n_groups_forest_local    =(int             *)SID_malloc(sizeof(int)             *n_forests_local);
  (*tree)->n_subgroups_forest_local =(int             *)SID_malloc(sizeof(int)             *n_forests_local);

  // Initialize
  int i_snap;
  int i_forest;
  for(i_snap=0;i_snap<n_snaps;i_snap++){
    (*tree)->n_groups_snap_local[i_snap]      =0;
    (*tree)->n_subgroups_snap_local[i_snap]   =0;
    (*tree)->first_neighbour_groups[i_snap]   =NULL;
    (*tree)->first_neighbour_subgroups[i_snap]=NULL;
    (*tree)->last_neighbour_groups[i_snap]    =NULL;
    (*tree)->last_neighbour_subgroups[i_snap] =NULL;
  }
  for(i_forest=0;i_forest<n_forests_local;i_forest++){
    (*tree)->n_groups_forest_local[i_forest]    =0;
    (*tree)->n_subgroups_forest_local[i_forest] =0;
    (*tree)->first_in_forest_groups[i_forest]   =NULL;
    (*tree)->first_in_forest_subgroups[i_forest]=NULL;
    (*tree)->last_in_forest_groups[i_forest]    =NULL;
    (*tree)->last_in_forest_subgroups[i_forest] =NULL;
  }

  SID_log("Done.",SID_LOG_CLOSE);
}

