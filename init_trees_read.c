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

void init_trees_read(const char  *filename_tree_root,
                     tree_info  **tree){

  SID_log("Initializing trees...",SID_LOG_OPEN);

  // Allocate the data structure
  (*tree)=(tree_info *)SID_malloc(sizeof(tree_info));
  (*tree)->snap_list                   =NULL;
  (*tree)->z_list                      =NULL;
  (*tree)->t_list                      =NULL;
  (*tree)->n_groups_snap_local         =NULL;
  (*tree)->n_subgroups_snap_local      =NULL;
  (*tree)->n_groups_forest_local       =NULL;
  (*tree)->n_subgroups_forest_local    =NULL;
  (*tree)->first_neighbour_groups      =NULL;
  (*tree)->first_neighbour_subgroups   =NULL;
  (*tree)->last_neighbour_groups       =NULL;
  (*tree)->last_neighbour_subgroups    =NULL;
  (*tree)->first_in_forest_groups      =NULL;
  (*tree)->first_in_forest_subgroups   =NULL;
  (*tree)->last_in_forest_groups       =NULL;
  (*tree)->last_in_forest_subgroups    =NULL;
  (*tree)->tree2forest_mapping_group   =NULL;
  (*tree)->tree2forest_mapping_subgroup=NULL;

  // Initialize filename paths
  char filename_input_file_root[MAX_FILENAME_LENGTH];
  char filename_input_dir_horizontal[MAX_FILENAME_LENGTH];
  char filename_input_dir_horizontal_trees[MAX_FILENAME_LENGTH];
  strcpy(filename_input_file_root,filename_tree_root);
  strip_path(filename_input_file_root);
  sprintf(filename_input_dir_horizontal,      "%s/horizontal",filename_tree_root);
  sprintf(filename_input_dir_horizontal_trees,"%s/trees",     filename_input_dir_horizontal);

  // Read the tree snap-range/search/scan parameters
  read_trees_run_parameters(filename_tree_root,
                            &((*tree)->i_read_start),
                            &((*tree)->i_read_stop),
                            &((*tree)->i_read_step),
                            &((*tree)->n_search),
                            &((*tree)->flag_fix_bridges),
                            &((*tree)->flag_compute_fragmented),
                            &((*tree)->flag_compute_ghosts));

  // Read the final halo and tree totals from the header of the last tree file
  read_trees_final_totals(filename_input_dir_horizontal_trees,
                          (*tree)->i_read_start,
                          (*tree)->i_read_stop,
                          (*tree)->i_read_step,
                          &((*tree)->i_read_last),
                          &((*tree)->n_snaps),
                          &((*tree)->n_groups_max),
                          &((*tree)->n_subgroups_max),
                          &((*tree)->n_progenitors_max),
                          &((*tree)->n_trees_subgroup),
                          &((*tree)->n_trees_group));

  // Compute/fetch the mapping between horizontal tree IDs and forest IDs ...
  read_forests(filename_tree_root,
               &((*tree)->n_forests_group),
               &((*tree)->n_forests_subgroup),
               &((*tree)->n_forests_group_local),
               &((*tree)->n_forests_subgroup_local),
               &((*tree)->tree2forest_mapping_group),
               &((*tree)->tree2forest_mapping_subgroup),
               &((*tree)->n_trees_forest_groups_max),
               &((*tree)->n_trees_forest_subgroups_max),
               &((*tree)->forest_lo_group_local),
               &((*tree)->forest_hi_group_local),
               &((*tree)->forest_lo_subgroup_local),
               &((*tree)->forest_hi_subgroup_local),
               &((*tree)->n_groups_local),
               &((*tree)->n_subgroups_local),
               &((*tree)->n_groups_snap_alloc_local),
               &((*tree)->n_subgroups_snap_alloc_local));

  // Set counts etc.  We set the number of forests to be
  //    the number of subgroup forests, since this the
  //    most generic/natural constraint
  int n_forests_local     =(*tree)->n_forests_subgroup_local;
  (*tree)->n_forests      =(*tree)->n_forests_subgroup;
  (*tree)->n_forests_local=n_forests_local;
  (*tree)->n_wrap         =(*tree)->n_search+1;

  // Create an array which maps the file numbers in the trees
  //   to the snapshot number (may differ from 1:1 if skipping snaps)
  int i_read;
  int i_file;
  int n_snaps     =(*tree)->n_snaps;
  int i_read_start=(*tree)->i_read_start;
  int i_read_stop =(*tree)->i_read_stop;
  int i_read_step =(*tree)->i_read_step;
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

  // Initialize look-up information
  (*tree)->group_indices          =NULL;
  (*tree)->group_array            =NULL;
  (*tree)->subgroup_indices       =NULL;
  (*tree)->subgroup_array         =NULL;

  // Initialize arrays
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
  ADaPS_init(&((*tree)->data));
  (*tree)->group_properties   =NULL;
  (*tree)->subgroup_properties=NULL;

  // Read snapshot expansion factor list
  char        filename_alist_in[MAX_FILENAME_LENGTH];
  FILE       *fp_alist_in;
  int         n_alist_in;
  int         i_alist;
  double      a_in;
  cosmo_info *cosmo;
  char       *line=NULL;
  size_t      line_length=0;
  init_cosmo_std(&cosmo);
  sprintf(filename_alist_in,"%s/a_list.txt",filename_tree_root);
  fp_alist_in=fopen(filename_alist_in,"r");
  n_alist_in=count_lines_data(fp_alist_in);
  if(n_alist_in!=(*tree)->n_snaps)
    SID_trap_error("The number of entries in the a_list.txt file does not make sense (ie. %d!=%d)",ERROR_LOGIC,n_alist_in,(*tree)->n_snaps);
  for(i_alist=0;i_alist<(*tree)->n_snaps;i_alist++){
     grab_next_line_data(fp_alist_in,&line,&line_length);
     grab_double(line,1,&a_in);
     (*tree)->z_list[i_alist]=z_of_a(a_in);
     (*tree)->t_list[i_alist]=deltat_a(&cosmo,DELTAT_A_MIN_A,a_in);
  }
  fclose(fp_alist_in);
  SID_free(SID_FARG line);
  free_cosmo(&cosmo);

  SID_log("Done.",SID_LOG_CLOSE);
}

