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
#include <gbpCosmo.h>

void read_trees(char       *filename_tree_root,
                char       *filename_halo_root,
                char       *filename_run_root,
                int         mode_progenitor,
                tree_info **trees){
  SID_fp      fp_out;
  SID_fp      fp_out_MBP;
  char        filename_in[256];
  char        filename_input_file_root[256];
  char        filename_input_dir_horizontal[256];
  char        filename_input_dir_horizontal_trees[256];
  char        filename_input_dir_horizontal_groups[256];
  char        filename_input_dir_horizontal_subgroups[256];
  char        filename_input_halos_groups[256];
  char        filename_input_halos_subgroups[256];
  int         i_bin;
  int         j_bin;
  int         i_write;
  int         n_files;
  int         i_read;
  int         j_read;
  int         k_read;
  int         i_file;
  int         j_file;
  int         i_group;
  int         i_subgroup;
  int         j_subgroup;
  int         k_match;
  int         i_rank;
  int         i_search;
  int         i_forest;
  int         k_subgroup;
  int         group_descendant_id;
  int         subgroup_descendant_id;
  int         group_file_offset;
  int         group_file_index;
  int         n_groups;
  int         n_subgroups;
  int         n_trees_subgroup;
  int         n_trees_group;
  int         n_forests_subgroup;
  int         n_forests_group;
  int         n_tree_bins_groups;
  int         n_tree_bins_subgroups;
  int         n_tree_bins_groups_file;
  int         n_tree_bins_subgroups_file;
  int         n_tree_bins_groups_rank;
  int         n_tree_bins_subgroups_rank;
  int         n_subgroups_group;
  int         group_id;
  int         group_tree_case;
  int         group_tree_id;
  int         subgroup_id;
  int         subgroup_tree_case;
  int         subgroup_tree_id;
  int         min_sum;
  int         min_bin;
  int         n_write;
  int         flag_match_subgroups;
  int        *n_halos_tree_group;
  int        *n_halos_tree_subgroup;
  int         n_halos_groups;
  int         n_halos_subgroups;
  int         n_halos_used;
  int         n_halos_target;
  int        *n_halos_forest;
  int         tree_count_group_local;  
  int        *tree_count_group_rank=NULL;
  int         tree_count_subgroup_local;
  int        *tree_count_subgroup_rank=NULL;
  int         n_trees_subgroup_local;
  int         n_trees_group_local;
  int         subgroup_file_offset;
  int         subgroup_file_index;
  char        group_text_prefix[4];
  int         n_conjoined;
  int         n_conjoined_total;
  int         depth_first_index;
  int         flag_write_init;
  int         k_tree;
  int         flag_init;
  int         progenitor_score;
  int         flag;
  int         max_id=0;
  int         n_halos_written;
  int         descendant_snap;
  halo_MBP_info              halo_MBP;
  tree_node_info   *current=NULL;
  tree_node_info   *last   =NULL;
  tree_node_info   *next   =NULL;

  // Read the tree search/scan parameters
  int i_read_start;
  int i_read_stop;
  int i_read_step;
  int n_search;
  int flag_fix_bridges;
  int flag_compute_fragmented;
  int flag_compute_ghosts;
  read_tree_run_parameters(filename_tree_root,
                           &i_read_start,
                           &i_read_stop,
                           &i_read_step,
                           &n_search,
                           &flag_fix_bridges,
                           &flag_compute_fragmented,
                           &flag_compute_ghosts);

  // We need i_read_start,i_read_stop,i_read_step from above before we can write this status message
  SID_log("Reading merger trees for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_start,i_read_stop,i_read_step);

  // Initialize filename paths
  strcpy(filename_input_file_root,filename_tree_root);
  strip_path(filename_input_file_root);
  sprintf(filename_input_dir_horizontal,          "%s/horizontal",filename_tree_root);
  sprintf(filename_input_dir_horizontal_groups,   "%s/groups",    filename_input_dir_horizontal);
  sprintf(filename_input_dir_horizontal_subgroups,"%s/subgroups", filename_input_dir_horizontal);
  sprintf(filename_input_dir_horizontal_trees,    "%s/trees",     filename_input_dir_horizontal);

  // Read the final halo and tree totals from the header of the last tree file
  int i_read_last;
  int n_snaps;
  int n_groups_max;
  int n_subgroups_max;
  int n_progenitors_max;
  read_tree_final_totals(filename_input_dir_horizontal_trees,
                         i_read_start,
                         i_read_stop, 
                         i_read_step,
                         &i_read_last,
                         &n_snaps,
                         &n_groups_max,
                         &n_subgroups_max,
                         &n_progenitors_max,
                         &n_trees_subgroup,
                         &n_trees_group);

  // Compute/fetch the mapping between horizontal tree IDs and forest IDs ...
  int  *i_forest_group;
  int  *i_forest_subgroup;
  int  *n_halos_forest_group;
  int  *n_halos_forest_subgroup;
  int   n_trees_forest_groups_max;
  int   n_trees_forest_subgroups_max;
  int   n_groups_local;
  int   n_subgroups_local;
  int   n_groups_snap_alloc_local;
  int   n_subgroups_snap_alloc_local;
  int   n_forests_group_local;
  int   n_forests_subgroup_local;
  int   forest_lo_group_local;
  int   forest_hi_group_local;
  int   forest_lo_subgroup_local;
  int   forest_hi_subgroup_local;
  read_forests(filename_tree_root,
               n_trees_group,
               n_trees_subgroup,
               &n_forests_group,
               &n_forests_subgroup,
               &n_forests_group_local,
               &n_forests_subgroup_local,
               &i_forest_group,
               &i_forest_subgroup,
               &n_halos_forest_group,
               &n_halos_forest_subgroup,
               &n_trees_forest_groups_max,
               &n_trees_forest_subgroups_max,
               &forest_lo_group_local,
               &forest_hi_group_local,
               &forest_lo_subgroup_local,
               &forest_hi_subgroup_local,
               &n_groups_local,
               &n_subgroups_local,
               &n_groups_snap_alloc_local,
               &n_subgroups_snap_alloc_local);

  // Initialize tree data structure
  int n_forests      =n_forests_subgroup;
  int n_forests_local=n_forests_subgroup_local;
  init_trees(i_read_start,
             i_read_stop,
             i_read_step,
             n_forests,
             n_forests_local,
             trees);

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
  if(n_alist_in!=(*trees)->n_snaps)
    SID_trap_error("The number of entries in the a_list.txt file does not make sense (ie. %d!=%d)",ERROR_LOGIC,n_alist_in,(*trees)->n_snaps);
  for(i_alist=0;i_alist<(*trees)->n_snaps;i_alist++){
     grab_next_line_data(fp_alist_in,&line,&line_length);
     grab_double(line,1,&a_in);
     (*trees)->z_list[i_alist]=z_of_a(a_in);
     (*trees)->t_list[i_alist]=deltat_a(&cosmo,DELTAT_A_MIN_A,a_in);
  }
  fclose(fp_alist_in);
  SID_free(SID_FARG line);
  free_cosmo(&cosmo);

  // Create index->pointer look-up arrays.  These are temporary arrays
  //   needed only for reading and will be deallocated afterwards.
  int             **group_indices;
  tree_node_info ***group_array;
  int             **subgroup_indices;
  tree_node_info ***subgroup_array;
  int i_wrap;
  int n_wrap=n_search+1;
  group_indices   =(int             **)SID_malloc(sizeof(int             *)*n_wrap);
  group_array     =(tree_node_info ***)SID_malloc(sizeof(tree_node_info **)*n_wrap);
  subgroup_indices=(int             **)SID_malloc(sizeof(int             *)*n_wrap);
  subgroup_array  =(tree_node_info ***)SID_malloc(sizeof(tree_node_info **)*n_wrap);
  for(i_wrap=0;i_wrap<n_wrap;i_wrap++){
     group_indices[i_wrap]   =(int             *)SID_malloc(sizeof(int)             *n_groups_snap_alloc_local);
     group_array[i_wrap]     =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_groups_snap_alloc_local);
     subgroup_indices[i_wrap]=(int             *)SID_malloc(sizeof(int)             *n_subgroups_snap_alloc_local);
     subgroup_array[i_wrap]  =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_subgroups_snap_alloc_local);
  }
  
  // Loop over all the horizontal tree files in order of decreasing snapshot number, hanging halos on the trees as we go
  SID_log("Reading trees...",SID_LOG_OPEN|SID_LOG_TIMER);
  int n_subgroups_orphaned=0;
  for(i_read=i_read_stop,i_file=n_snaps-1,flag_init=TRUE;i_read>=i_read_start;i_read-=i_read_step,i_file--,flag_init=FALSE){
    SID_log("Processing snapshot %03d (%03d of %03d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read,i_file+1,n_snaps);

    // Open horizontal tree file
    SID_fp fp_trees_in;
    sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_input_dir_horizontal_trees,i_read);
    SID_fopen(filename_in,"r",&fp_trees_in);

    // Open the halo files and read their header (needed for halo particle counts)
    int    n_groups_cat;
    int    n_subgroups_cat;
    int    offset_size;
    SID_fp fp_groups_in;
    SID_fp fp_subgroups_in;
    sprintf(filename_input_halos_groups,   "%s_%03d.catalog_groups",   filename_halo_root,i_read);
    sprintf(filename_input_halos_subgroups,"%s_%03d.catalog_subgroups",filename_halo_root,i_read);
    SID_fopen(filename_input_halos_groups,"r",   &fp_groups_in);
    SID_fread_all(&n_groups_cat,sizeof(int),1,   &fp_groups_in);
    SID_fread_all(&offset_size, sizeof(int),1,   &fp_groups_in);
    SID_fopen(filename_input_halos_subgroups,"r",&fp_subgroups_in);
    SID_fread_all(&n_subgroups_cat,sizeof(int),1,&fp_subgroups_in);
    SID_fread_all(&offset_size,    sizeof(int),1,&fp_subgroups_in);

    // Read tree file header
    int n_step_in;
    int n_search_in;
    int n_groups_max_in;
    int n_subgroups_max_in;
    SID_fread_all(&n_step_in,         sizeof(int),1,&fp_trees_in);
    SID_fread_all(&n_search_in,       sizeof(int),1,&fp_trees_in);
    SID_fread_all(&n_groups,          sizeof(int),1,&fp_trees_in);
    SID_fread_all(&n_subgroups,       sizeof(int),1,&fp_trees_in);
    SID_fread_all(&n_groups_max_in,   sizeof(int),1,&fp_trees_in);
    SID_fread_all(&n_subgroups_max_in,sizeof(int),1,&fp_trees_in);
    SID_fread_all(&n_trees_subgroup,  sizeof(int),1,&fp_trees_in);
    SID_fread_all(&n_trees_group,     sizeof(int),1,&fp_trees_in);
    n_progenitors_max=MAX(n_groups_max_in,n_subgroups_max_in);
    if(n_step_in!=i_read_step) 
       SID_trap_error("Snapshot step sizes don't match (ie. %d!=%d)",ERROR_LOGIC,n_step_in,i_read_step);
    if(n_groups_cat!=n_groups)
       SID_trap_error("Group counts don't match between datasets (ie. %d!=%d)",ERROR_LOGIC,n_groups_cat,n_groups);
    if(n_subgroups_cat!=n_subgroups)
       SID_trap_error("Subgroup counts don't match between datasets (ie. %d!=%d)",ERROR_LOGIC,n_subgroups_cat,n_subgroups);
    
    // Read each group in turn
    int tree_read_buffer[7];
    int n_groups_unused        =0;
    int n_groups_added_multiply=0;
    for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){

      // Read horizontal trees for groups
      int n_particles_group;
      SID_fread_all(&n_particles_group, sizeof(int),1,&fp_groups_in);
      SID_fread_all(tree_read_buffer, 7*sizeof(int),1,&fp_trees_in);
      group_id           =tree_read_buffer[0];
      group_tree_case    =tree_read_buffer[1];
      group_descendant_id=tree_read_buffer[2];
      group_tree_id      =tree_read_buffer[3];
      group_file_offset  =tree_read_buffer[4];
      group_file_index   =tree_read_buffer[5];
      n_subgroups_group  =tree_read_buffer[6];

      // Do some interpreting of the group tree information
      int group_descendant_snap;
      int flag_group_added=FALSE;
      int group_forest_id;
      if(group_tree_id>=0)
        group_forest_id=i_forest_group[group_tree_id];
      else
        group_forest_id=-1;
      if(group_forest_id>=0){
         if(group_file_offset<=0)
           group_descendant_snap=-1; // Needed for halos in the root snapshot
         else
           group_descendant_snap=(i_file+group_file_offset);
      }

      // Read each subgroup in turn
      int i_subgroup_valid;
      tree_node_info *group_node;
      tree_node_info *subgroup_node;
      for(j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
        int n_particles_subgroup;
        SID_fread_all(&n_particles_subgroup,  sizeof(int),1,&fp_subgroups_in);
        SID_fread_all(tree_read_buffer,     6*sizeof(int),1,&fp_trees_in);
        subgroup_id           =tree_read_buffer[0];
        subgroup_tree_case    =tree_read_buffer[1];
        subgroup_descendant_id=tree_read_buffer[2];
        subgroup_tree_id      =tree_read_buffer[3];
        subgroup_file_offset  =tree_read_buffer[4];
        subgroup_file_index   =tree_read_buffer[5];

        // Ignore negative ids
        int subgroup_forest_id;
        if(subgroup_tree_id>=0)
          subgroup_forest_id=i_forest_subgroup[subgroup_forest_id];
        else
          subgroup_forest_id=-1;

        // Add nodes to trees
        if(subgroup_forest_id>=0){
          // If this subgroup belongs to a local forest ...
          i_forest=subgroup_forest_id-forest_lo_subgroup_local;
          if(i_forest>=0 && i_forest<n_forests_subgroup_local){ 
            // ... add the group ...
            if(!flag_group_added){
               add_node_to_trees((*trees),               // The tree datastructure
                                 i_forest,               // Local forest index
                                 group_tree_case,        // Halo's TREE_CASE BWS
                                 n_particles_group,      // Number of particles
                                 i_file,                 // Halo's tree snapshot number
                                 i_group,                // Halo's file index
                                 group_descendant_snap,  // Descendant's snap
                                 group_file_index,       // Descendant's index
                                 group_indices,          // Rolling look-up array for index->pointer mapping
                                 group_array,            // Rolling look-up array for index->pointer mapping
                                 n_wrap,                 // Rolling look-up array size
                                 NULL,                   // Pointer to the new node's group
                                 &group_node);           // Pointer to the new node
               flag_group_added=TRUE;
            }
            // ... create a new branch and add it to its tree ...
            int subgroup_descendant_snap;
            if(subgroup_file_offset<=0)
              subgroup_descendant_snap=-1; // Needed for halos in the root snapshot
            else
              subgroup_descendant_snap=(i_file+subgroup_file_offset);
            // ... add this halo to the trees ...
            add_node_to_trees((*trees),                 // The tree datastructure
                              i_forest,                 // Local forest index
                              subgroup_tree_case,       // Halo's TREE_CASE BWS
                              n_particles_subgroup,     // Number of particles
                              i_file,                   // Halo's tree snapshot number
                              i_subgroup,               // Halo's file index
                              subgroup_descendant_snap, // Descendant's snap
                              subgroup_file_index,      // Descendant's index
                              subgroup_indices,         // Rolling look-up array for index->pointer mapping
                              subgroup_array,           // Rolling look-up array for index->pointer mapping
                              n_wrap,                   // Rolling look-up array size
                              group_node,               // Pointer to the new node's group
                              &subgroup_node);          // Pointer to the new node
          }
        }
      }

      // Check how many ranks have used this group.  Should be
      //   just one in all but the most pathological cases.
      SID_Allreduce(SID_IN_PLACE,&flag_group_added,1,SID_INT,SID_SUM,SID.COMM_WORLD);
      if(flag_group_added==0)
         n_groups_unused++;
      else if(flag_group_added>1)
         n_groups_added_multiply++;

    } // i_group
    SID_fclose(&fp_trees_in);
    SID_fclose(&fp_groups_in);
    SID_fclose(&fp_subgroups_in);

    // Update the temporary look-up arrays
    int             i_halo;
    tree_node_info *current;
    i_halo =0;
    current=(*trees)->first_neighbour_groups[i_file];
    while(current!=NULL){
       group_indices[i_file%n_wrap][i_halo]=current->file_index;
       group_array[i_file%n_wrap][i_halo]  =current;
       i_halo++;
       current=current->next_neighbour;
    }
    i_halo =0;
    current=(*trees)->first_neighbour_subgroups[i_file];
    while(current!=NULL){
       subgroup_indices[i_file%n_wrap][i_halo]=current->file_index;
       subgroup_array[i_file%n_wrap][i_halo]  =current;
       i_halo++;
       current=current->next_neighbour;
    }
 
    // Report some group statistics
    if(n_groups_unused>0)
       SID_log("%d groups unused...",SID_LOG_CONTINUE,n_groups_unused);
    if(n_groups_added_multiply>0)
       SID_log("%d groups used by multiple cores...",SID_LOG_CONTINUE,n_groups_added_multiply);

    SID_log("Done.",SID_LOG_CLOSE);
  } // i_read
  SID_log("Done.",SID_LOG_CLOSE);

  // Finalize trees
  finalize_trees((*trees),mode_progenitor);

  // Clean-up
  for(i_wrap=0;i_wrap<n_wrap;i_wrap++){
     SID_free(SID_FARG group_indices[i_wrap]);
     SID_free(SID_FARG group_array[i_wrap]);
     SID_free(SID_FARG subgroup_indices[i_wrap]);
     SID_free(SID_FARG subgroup_array[i_wrap]);
  }
  SID_free(SID_FARG group_indices);
  SID_free(SID_FARG group_array);
  SID_free(SID_FARG subgroup_indices);
  SID_free(SID_FARG subgroup_array);
  SID_free(SID_FARG i_forest_group);
  SID_free(SID_FARG i_forest_subgroup);
  SID_free(SID_FARG n_halos_forest_group);
  SID_free(SID_FARG n_halos_forest_subgroup);

  SID_log("Done.",SID_LOG_CLOSE);
}

