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

void init_trees_read(const char  *filename_SSimPL_dir,
                     const char  *filename_trees_version,
                     int          mode,
                     tree_info  **tree){

  SID_log("Initializing trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Set the tree filename root
  char filename_trees_root[MAX_FILENAME_LENGTH];
  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_dir,filename_trees_version);

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
  char tree_name[MAX_FILENAME_LENGTH];
  strcpy(tree_name,filename_trees_root);
  strip_path(tree_name);
  sprintf((*tree)->filename_root,"%s",filename_trees_root);
  sprintf((*tree)->name,         "%s",tree_name);
  sprintf((*tree)->filename_root_analysis,        "%s/analysis",  (*tree)->filename_root);
  sprintf((*tree)->filename_root_horizontal,      "%s/horizontal",(*tree)->filename_root);
  sprintf((*tree)->filename_root_horizontal_trees,"%s/trees",     (*tree)->filename_root_horizontal);

  // Read the tree snap-range/search/scan parameters
  read_trees_run_parameters(filename_trees_root,
                            &((*tree)->i_read_start),
                            &((*tree)->i_read_stop),
                            &((*tree)->i_read_step),
                            &((*tree)->n_search),
                            &((*tree)->flag_fix_bridges),
                            &((*tree)->flag_compute_fragmented),
                            &((*tree)->flag_compute_ghosts));

  // Read the final halo and tree totals from the header of the last tree file
  read_trees_final_totals((*tree)->filename_root_horizontal_trees,
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
  (*tree)->n_wrap=(*tree)->n_search+1;

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

  // Initialize the cosmology
  char filename_cosmology[MAX_FILENAME_LENGTH];
  sprintf(filename_cosmology,"%s/run/cosmology.txt",filename_SSimPL_dir);
  read_gbpCosmo_file(&((*tree)->cosmo),filename_cosmology);
  cosmo_info *cosmo=(*tree)->cosmo;

  // Read snapshot expansion factor list
  char        filename_alist_in[MAX_FILENAME_LENGTH];
  FILE       *fp_alist_in;
  int         n_alist_in;
  int         i_alist;
  double      a_in;
  char       *line=NULL;
  size_t      line_length=0;
  sprintf(filename_alist_in,"%s/a_list.txt",filename_trees_root);
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

  // Initialize a bunch of stuff
  (*tree)->n_forests                =0;
  (*tree)->n_forests_local          =0;
  (*tree)->n_groups_trees           =0;
  (*tree)->n_groups_trees_local     =0;
  (*tree)->n_subgroups_trees_local  =0;
  (*tree)->n_subgroups_trees        =0;
  (*tree)->group_indices            =NULL;
  (*tree)->group_array              =NULL;
  (*tree)->subgroup_indices         =NULL;
  (*tree)->subgroup_array           =NULL;
  (*tree)->group_markers            =NULL;
  (*tree)->subgroup_markers         =NULL;
  (*tree)->group_properties         =NULL;
  (*tree)->subgroup_properties      =NULL;
  (*tree)->n_groups_catalog         =NULL;
  (*tree)->n_subgroups_catalog      =NULL;
  (*tree)->n_groups_snap_local      =NULL;
  (*tree)->n_subgroups_snap_local   =NULL;
  (*tree)->first_neighbour_groups   =NULL;
  (*tree)->first_neighbour_subgroups=NULL;
  (*tree)->last_neighbour_groups    =NULL;
  (*tree)->last_neighbour_subgroups =NULL;
  (*tree)->first_in_forest_groups   =NULL;
  (*tree)->first_in_forest_subgroups=NULL;
  (*tree)->last_in_forest_groups    =NULL;
  (*tree)->last_in_forest_subgroups =NULL;
  (*tree)->n_groups_forest_local    =NULL;
  (*tree)->n_subgroups_forest_local =NULL;
  ADaPS_init(&((*tree)->data));

  // If we only want the header information, then don't do the rest
  if(!check_mode_for_flag(mode,TREE_READ_HEADER_ONLY)){
     // Compute/fetch the mapping between horizontal tree IDs and forest IDs ...
     int n_forests_group_in;
     int n_forests_subgroup_in;
     int n_forests_group_local_in;
     int n_forests_subgroup_local_in;
     read_forests(filename_trees_root,
                  &n_forests_group_in,
                  &n_forests_subgroup_in,
                  &n_forests_group_local_in,
                  &n_forests_subgroup_local_in,
                  &((*tree)->tree2forest_mapping_group),
                  &((*tree)->tree2forest_mapping_subgroup),
                  &((*tree)->n_trees_forest_groups_max),
                  &((*tree)->n_trees_forest_subgroups_max),
                  &((*tree)->forest_lo_group_local),
                  &((*tree)->forest_hi_group_local),
                  &((*tree)->forest_lo_subgroup_local),
                  &((*tree)->forest_hi_subgroup_local),
                  &((*tree)->n_groups_raw_local),
                  &((*tree)->n_subgroups_raw_local),
                  &((*tree)->n_groups_snap_alloc_local),
                  &((*tree)->n_subgroups_snap_alloc_local));
     calc_sum_global(&((*tree)->n_groups_raw_local),   &((*tree)->n_groups_raw),   1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
     calc_sum_global(&((*tree)->n_subgroups_raw_local),&((*tree)->n_subgroups_raw),1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);

     // Set counts etc.  We set the number of forests to be
     //    the number of subgroup forests, since this is the
     //    most generic/natural constraint
     int n_forests_local             =n_forests_subgroup_local_in;
     (*tree)->n_forests              =n_forests_subgroup_in;
     (*tree)->n_forests_local        =n_forests_local;
     (*tree)->n_groups_trees         =0;
     (*tree)->n_groups_trees_local   =0;
     (*tree)->n_subgroups_trees_local=0;
     (*tree)->n_subgroups_trees      =0;

     // Allocate counters
     (*tree)->n_groups_catalog      =(int *)SID_calloc(sizeof(int)*n_snaps);
     (*tree)->n_subgroups_catalog   =(int *)SID_calloc(sizeof(int)*n_snaps);
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
  }

  SID_log("Done.",SID_LOG_CLOSE);
}

