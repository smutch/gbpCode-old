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
#include <gbpCosmo.h>

void read_trees(const char *filename_SSimPL_root,
                const char *filename_halos_version,
                const char *filename_trees_version,
                int         read_mode,
                tree_info **trees){

  // Set the halo and tree filename roots
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_root,filename_trees_version);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_root,filename_halos_version);

  // Read the tree search/scan parameters so we can write the opening log statement
  int i_read_start_temp;
  int i_read_stop_temp;
  int i_read_step_temp;
  int n_search_temp;
  int flag_fix_bridges_temp;
  int flag_compute_fragmented_temp;
  int flag_compute_ghosts_temp;
  read_trees_run_parameters(filename_trees_root,
                            &i_read_start_temp,
                            &i_read_stop_temp,
                            &i_read_step_temp,
                            &n_search_temp,
                            &flag_fix_bridges_temp,
                            &flag_compute_fragmented_temp,
                            &flag_compute_ghosts_temp);

  // We need i_read_start,i_read_stop,i_read_step from above before we can write this status message
  SID_log("Reading merger trees for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_start_temp,i_read_stop_temp,i_read_step_temp);

  // Initialize tree data structure and populate it
  //   with various pieces of header information
  init_trees_read(filename_SSimPL_root,filename_halos_version,filename_trees_version,TREE_READ_DEFAULT,trees);

  // Set the mode
  (*trees)->mode=read_mode;

  // To make the code look a little cleaner, create some aliases
  int i_read_start=(*trees)->i_read_start;
  int i_read_stop =(*trees)->i_read_stop;
  int i_read_step =(*trees)->i_read_step;
  int i_read_last =(*trees)->i_read_last;
  int n_snaps     =(*trees)->n_snaps;
  int n_search    =(*trees)->n_search;

  // Initialize filename paths
  char filename_input_file_root[MAX_FILENAME_LENGTH];
  char filename_input_dir_horizontal[MAX_FILENAME_LENGTH];
  char filename_input_dir_horizontal_groups[MAX_FILENAME_LENGTH];
  char filename_input_dir_horizontal_subgroups[MAX_FILENAME_LENGTH];
  char filename_input_dir_horizontal_trees[MAX_FILENAME_LENGTH];
  strcpy(filename_input_file_root,filename_trees_root);
  strip_path(filename_input_file_root);
  sprintf(filename_input_dir_horizontal,          "%s/horizontal",filename_trees_root);
  sprintf(filename_input_dir_horizontal_groups,   "%s/groups",    filename_input_dir_horizontal);
  sprintf(filename_input_dir_horizontal_subgroups,"%s/subgroups", filename_input_dir_horizontal);
  sprintf(filename_input_dir_horizontal_trees,    "%s/trees",     filename_input_dir_horizontal);

  // Create index->pointer look-up arrays.  These are temporary arrays
  //   needed only for reading and will be deallocated afterwards, unless
  //   these trees are being read with the TREE_MODE_REFERENCE flag,
  //   in which case they will be kept for doing reference halo lookups
  init_trees_lookup((*trees));

  // Initialize the extended tree pointers if we are reading them
  int               flag_read_extended_pointers        =FALSE;
  tree_node_info ***forematch_pointers_groups_local    =NULL;
  tree_node_info ***backmatch_pointers_groups_local    =NULL;
  tree_node_info ***forematch_pointers_subgroups_local =NULL;
  tree_node_info ***backmatch_pointers_subgroups_local =NULL;
  float           **descendant_pointers_groups_score   =NULL;
  float           **progenitor_pointers_groups_score   =NULL;
  float           **forematch_pointers_groups_score    =NULL;
  float           **backmatch_pointers_groups_score    =NULL;
  float           **descendant_pointers_subgroups_score=NULL;
  float           **progenitor_pointers_subgroups_score=NULL;
  float           **forematch_pointers_subgroups_score =NULL;
  float           **backmatch_pointers_subgroups_score =NULL;
  if(check_mode_for_flag(read_mode,TREE_READ_EXTENDED_POINTERS)){
     // Since the halo counts for each snap haven't been set yet, each of these
     //    just allocate an array of pointers, with one element per snap set to NULL.
     //    The arrays for each snap will need to be allocated later.
     init_trees_data((*trees),(void ***)(&forematch_pointers_groups_local),    0,INIT_TREE_DATA_GROUPS,   "forematch_pointers_groups_local");
     init_trees_data((*trees),(void ***)(&backmatch_pointers_groups_local),    0,INIT_TREE_DATA_GROUPS,   "backmatch_pointers_groups_local");
     init_trees_data((*trees),(void ***)(&forematch_pointers_subgroups_local), 0,INIT_TREE_DATA_SUBGROUPS,"forematch_pointers_subgroups_local");
     init_trees_data((*trees),(void ***)(&backmatch_pointers_subgroups_local), 0,INIT_TREE_DATA_SUBGROUPS,"backmatch_pointers_subgroups_local");
     init_trees_data((*trees),(void ***)(&descendant_pointers_groups_score),   0,INIT_TREE_DATA_GROUPS,   "descendant_pointers_groups_score");
     init_trees_data((*trees),(void ***)(&progenitor_pointers_groups_score),   0,INIT_TREE_DATA_GROUPS,   "progenitor_pointers_groups_score");
     init_trees_data((*trees),(void ***)(&forematch_pointers_groups_score),    0,INIT_TREE_DATA_GROUPS,   "forematch_pointers_groups_score");
     init_trees_data((*trees),(void ***)(&backmatch_pointers_groups_score),    0,INIT_TREE_DATA_GROUPS,   "backmatch_pointers_groups_score");
     init_trees_data((*trees),(void ***)(&descendant_pointers_subgroups_score),0,INIT_TREE_DATA_SUBGROUPS,"descendant_pointers_subgroups_score");
     init_trees_data((*trees),(void ***)(&progenitor_pointers_subgroups_score),0,INIT_TREE_DATA_SUBGROUPS,"progenitor_pointers_subgroups_score");
     init_trees_data((*trees),(void ***)(&forematch_pointers_subgroups_score), 0,INIT_TREE_DATA_SUBGROUPS,"forematch_pointers_subgroups_score");
     init_trees_data((*trees),(void ***)(&backmatch_pointers_subgroups_score), 0,INIT_TREE_DATA_SUBGROUPS,"backmatch_pointers_subgroups_score");
     (*trees)->group_backmatch_pointers   =backmatch_pointers_groups_local;
     (*trees)->subgroup_backmatch_pointers=backmatch_pointers_subgroups_local;
     (*trees)->group_forematch_pointers   =forematch_pointers_groups_local;
     (*trees)->subgroup_forematch_pointers=forematch_pointers_subgroups_local;
     (*trees)->group_descendant_score     =descendant_pointers_groups_score;
     (*trees)->subgroup_descendant_score  =descendant_pointers_subgroups_score;
     (*trees)->group_progenitor_score     =progenitor_pointers_groups_score;
     (*trees)->subgroup_progenitor_score  =progenitor_pointers_subgroups_score;
     (*trees)->group_backmatch_score      =backmatch_pointers_groups_score;
     (*trees)->subgroup_backmatch_score   =backmatch_pointers_subgroups_score;
     (*trees)->group_forematch_score      =forematch_pointers_groups_score;
     (*trees)->subgroup_forematch_score   =forematch_pointers_subgroups_score;
     flag_read_extended_pointers=TRUE;
  }

  // Check if the halo files have substructure information.  If they do, we will read it
  //    as we read the trees to set the substructure hierarchy pointers.
  char  filename_in[MAX_FILENAME_LENGTH];
  int   flag_read_sub_pointers=FALSE;
  if(SID.I_am_Master){
     FILE *fp_test=NULL;
     sprintf(filename_in,"%s/halos/%s_%03d.catalog_subgroups",filename_SSimPL_root,filename_halos_version,i_read_stop);
     fp_test=fopen(filename_in,"r");
     int n_subgroups_in;fread_verify(&n_subgroups_in,sizeof(int),1,fp_test);
     int offset_size;   fread_verify(&offset_size,   sizeof(int),1,fp_test);
     fseeko(fp_test,n_subgroups_in*(sizeof(int)+offset_size),SEEK_CUR);
     int test;          fread(&offset_size,sizeof(int),1,fp_test); // make sure not to verify here!
     flag_read_sub_pointers=!feof(fp_test);
     fclose(fp_test);
  }
  SID_Bcast(&flag_read_sub_pointers,sizeof(int),MASTER_RANK,SID.COMM_WORLD);
  if(flag_read_sub_pointers){
     SID_log("Substructure hierarchy pointers present and will be used.",SID_LOG_COMMENT);
     (*trees)->mode|=TREE_MODE_SUBSTRUCTURE_HIERARCHY_ON;
  }
  else
     SID_log("Substructure hierarchy pointers not present and will not be used.",SID_LOG_COMMENT);

  // Maintain a list of substructure pointers for
  //    setting substructure hierarchy pointers
  tree_node_info **group_substructure_pointers=NULL;
  int n_hierarchy_alloc=0;
  if(flag_read_sub_pointers){
     n_hierarchy_alloc          =32;
     group_substructure_pointers=(tree_node_info **)SID_malloc(n_hierarchy_alloc*sizeof(tree_node_info *));
  }

  // Loop over all the horizontal tree files in order of decreasing snapshot number, hanging halos on the trees as we go
  //    For the back match pointer reading, the halos being pointed to must be in place before the read, so
  //    we have to wait until (n_search+1) snapshots have been processed before we start reading those files.
  SID_log("Reading trees...",SID_LOG_OPEN|SID_LOG_TIMER);
  int n_subgroups_orphaned=0;
  int i_read           =i_read_stop;
  int i_read_bridge    =i_read_stop;
  int i_read_backmatch =i_read_stop+n_search*i_read_step;
  int i_read_progenitor=i_read_stop+n_search*i_read_step;
  int i_file           =n_snaps-1;
  int i_file_bridge    =n_snaps-1;
  int i_file_backmatch =n_snaps-1+n_search;
  int i_file_progenitor=n_snaps-1+n_search;
  for(;
      i_read>=i_read_start;
      i_read           -=i_read_step,
      i_read_bridge    -=i_read_step,
      i_read_backmatch -=i_read_step,
      i_read_progenitor-=i_read_step,
      i_file--,
      i_file_bridge--,
      i_file_backmatch--,
      i_file_progenitor--){
    SID_log("Processing snapshot %03d (%03d of %03d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read,i_file+1,n_snaps);

    // Open horizontal tree file
    SID_fp fp_trees_in;
    sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_input_dir_horizontal_trees,i_read);
    SID_fopen(filename_in,"r",&fp_trees_in);

    // Open the halo files and read their header (needed for halo particle counts)
    char   filename_input_halos_groups[MAX_FILENAME_LENGTH];
    char   filename_input_halos_subgroups[MAX_FILENAME_LENGTH];
    int    n_groups_cat;
    int    n_subgroups_cat;
    int    offset_size;
    SID_fp fp_groups_in;
    SID_fp fp_subgroups_in;
    sprintf(filename_input_halos_groups,   "%s_%03d.catalog_groups",   filename_halos_root,i_read);
    sprintf(filename_input_halos_subgroups,"%s_%03d.catalog_subgroups",filename_halos_root,i_read);
    SID_fopen(filename_input_halos_groups,   "r",&fp_groups_in);
    SID_fopen(filename_input_halos_subgroups,"r",&fp_subgroups_in);
    SID_fread_all(&n_groups_cat,   sizeof(int),1,&fp_groups_in);
    SID_fread_all(&offset_size,    sizeof(int),1,&fp_groups_in);
    SID_fread_all(&n_subgroups_cat,sizeof(int),1,&fp_subgroups_in);
    SID_fread_all(&offset_size,    sizeof(int),1,&fp_subgroups_in);
    (*trees)->n_groups_catalog[i_file]   =n_groups_cat;
    (*trees)->n_subgroups_catalog[i_file]=n_subgroups_cat;

    // Set a pointer to the substructure hierarchy pointers if they are being read
    SID_fp fp_hierarchy_in;
    if(flag_read_sub_pointers){
       int n_subgroups_in;
       int offset_size_in;
       SID_fopen(filename_input_halos_subgroups,"r",&fp_hierarchy_in);
       SID_fread_all(&n_subgroups_in,sizeof(int),1,&fp_hierarchy_in);
       SID_fread_all(&offset_size_in,sizeof(int),1,&fp_hierarchy_in);
       SID_fseek(&fp_hierarchy_in,sizeof(int),   n_subgroups_in,SID_SEEK_CUR);
       SID_fseek(&fp_hierarchy_in,offset_size_in,n_subgroups_in,SID_SEEK_CUR);
    }

    // Read tree file header
    int n_step_in;
    int n_search_in;
    int n_groups;
    int n_subgroups;
    int n_groups_max_in;
    int n_subgroups_max_in;
    int n_trees_subgroup;
    int n_trees_group;
    int n_progenitors_max;
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

    // Initialize read buffers
    SID_fp_buffer *fp_hierarchy_in_buffer=NULL;
    SID_fp_buffer *fp_subgroups_in_buffer=NULL;
    SID_fp_buffer *fp_groups_in_buffer   =NULL;
    SID_fp_buffer *fp_trees_in_buffer    =NULL;
    size_t n_bytes_trees                 =sizeof(int)*((8*(size_t)n_groups)+(7*(size_t)n_subgroups));
    init_SID_fp_buffer(&fp_subgroups_in,(size_t)n_subgroups*sizeof(int),SIZE_OF_MEGABYTE,&fp_subgroups_in_buffer);
    init_SID_fp_buffer(&fp_groups_in,   (size_t)n_groups   *sizeof(int),SIZE_OF_MEGABYTE,&fp_groups_in_buffer);
    init_SID_fp_buffer(&fp_trees_in,    n_bytes_trees,                  SIZE_OF_MEGABYTE,&fp_trees_in_buffer);
    if(flag_read_sub_pointers)
       init_SID_fp_buffer(&fp_hierarchy_in,(size_t)n_subgroups*sizeof(int),SIZE_OF_MEGABYTE,&fp_hierarchy_in_buffer);

    // Read each group in turn
    int    n_groups_unused        =0;
    int    n_groups_added_multiply=0;
    int    i_group                =0;
    int    i_subgroup             =0;
    for(;i_group<n_groups;i_group++){

      // Read horizontal trees for groups
      int n_particles_group;     SID_fread_all_buffer(&n_particles_group,     sizeof(int),1,fp_groups_in_buffer);
      int group_id;              SID_fread_all_buffer(&group_id,              sizeof(int),1,fp_trees_in_buffer);
      int group_tree_case;       SID_fread_all_buffer(&group_tree_case,       sizeof(int),1,fp_trees_in_buffer);
      int group_descendant_id;   SID_fread_all_buffer(&group_descendant_id,   sizeof(int),1,fp_trees_in_buffer);
      int group_tree_id;         SID_fread_all_buffer(&group_tree_id,         sizeof(int),1,fp_trees_in_buffer);
      int group_file_offset;     SID_fread_all_buffer(&group_file_offset,     sizeof(int),1,fp_trees_in_buffer);
      int group_file_index;      SID_fread_all_buffer(&group_file_index,      sizeof(int),1,fp_trees_in_buffer);
      int group_n_particles_peak;SID_fread_all_buffer(&group_n_particles_peak,sizeof(int),1,fp_trees_in_buffer);
      int n_subgroups_group;     SID_fread_all_buffer(&n_subgroups_group,     sizeof(int),1,fp_trees_in_buffer);

      // Add offset to snap
      int group_descendant_snap;
      int flag_group_added=FALSE;
      if(group_file_offset<=0)
        group_descendant_snap=-1; // Needed for halos in the root snapshot
      else
        group_descendant_snap=(i_file+group_file_offset);

      // If we are reading the substructure hieararchy, maintain a list
      //    of pointers for this group which we can use to set these pointers
      if(flag_read_sub_pointers && n_subgroups_group>n_hierarchy_alloc){
         SID_free(SID_FARG group_substructure_pointers);
         n_hierarchy_alloc=MAX(2*n_hierarchy_alloc,n_subgroups_group);
         group_substructure_pointers=(tree_node_info **)SID_malloc(n_hierarchy_alloc*sizeof(tree_node_info *));
      }

      // Read each subgroup in turn
      tree_node_info *group_node   =NULL;
      tree_node_info *subgroup_node=NULL;
      for(int j_subgroup=0;j_subgroup<n_subgroups_group;j_subgroup++){
         int n_particles_subgroup;     SID_fread_all_buffer(&n_particles_subgroup,     sizeof(int),1,fp_subgroups_in_buffer);
         int subgroup_id;              SID_fread_all_buffer(&subgroup_id,              sizeof(int),1,fp_trees_in_buffer);
         int subgroup_tree_case;       SID_fread_all_buffer(&subgroup_tree_case,       sizeof(int),1,fp_trees_in_buffer);
         int subgroup_descendant_id;   SID_fread_all_buffer(&subgroup_descendant_id,   sizeof(int),1,fp_trees_in_buffer);
         int subgroup_tree_id;         SID_fread_all_buffer(&subgroup_tree_id,         sizeof(int),1,fp_trees_in_buffer);
         int subgroup_file_offset;     SID_fread_all_buffer(&subgroup_file_offset,     sizeof(int),1,fp_trees_in_buffer);
         int subgroup_file_index;      SID_fread_all_buffer(&subgroup_file_index,      sizeof(int),1,fp_trees_in_buffer);
         int subgroup_n_particles_peak;SID_fread_all_buffer(&subgroup_n_particles_peak,sizeof(int),1,fp_trees_in_buffer);

         // Add offset to snap
         int subgroup_descendant_snap;
         if(subgroup_file_offset<=0)
            subgroup_descendant_snap=-1; // Needed for halos in the root snapshot
         else
            subgroup_descendant_snap=(i_file+subgroup_file_offset);

         // Ignore halos with undefined tree_ids
         int subgroup_forest_id;
         int i_forest;
         if(subgroup_tree_id>=0){
            subgroup_forest_id=(*trees)->tree2forest_mapping_subgroup[subgroup_tree_id];
            i_forest          =subgroup_forest_id-(*trees)->forest_lo_subgroup_local;
         }
         else{
            subgroup_forest_id=-1;
            i_forest          =-1;
         }

         // Add node to trees if this subgroup belongs to a local forest ...
         subgroup_node=NULL;
         if(i_forest>=0 && i_forest<(*trees)->n_forests_local){ 
            // ... add the group ...
            if(!flag_group_added){
               add_node_to_trees((*trees),               // The tree datastructure
                                 i_forest,               // Local forest index
                                 group_tree_case,        // Halo's TREE_CASE BWS
                                 n_particles_group,      // Number of particles
                                 group_id,               // Halo's tree ID
                                 i_file,                 // Halo's tree snapshot number
                                 i_group,                // Halo's file index
                                 group_descendant_snap,  // Descendant's snap
                                 group_file_index,       // Descendant's index
                                 NULL,                   // Pointer to this halo's group.  NULL if halo is a group.
                                 &group_node);           // Pointer to the new node
               group_node->n_particles_peak=group_n_particles_peak;
               flag_group_added=TRUE;
            }
            // ... add the subgroup ...
            add_node_to_trees((*trees),                 // The tree datastructure
                              i_forest,                 // Local forest index
                              subgroup_tree_case,       // Halo's TREE_CASE BWS
                              n_particles_subgroup,     // Number of particles
                              subgroup_id,              // Halo's tree ID
                              i_file,                   // Halo's tree snapshot number
                              i_subgroup+j_subgroup,    // Halo's file index
                              subgroup_descendant_snap, // Descendant's snap
                              subgroup_file_index,      // Descendant's index
                              group_node,               // Pointer to this halo's group.  NULL if halo is a group.
                              &subgroup_node);          // Pointer to the new node
            subgroup_node->n_particles_peak=subgroup_n_particles_peak;
         }
         // Add new pointer to substructure list (if needed)
         if(flag_read_sub_pointers)
            group_substructure_pointers[j_subgroup]=subgroup_node;
         else
            add_node_to_substructure_hierarchy((*trees),      // The tree datastructure
                                               subgroup_node, // Pointer to the substructure being added
                                               group_node);   // Pointer to its parent
      }

      // Check how many ranks have used this group.  Should be just one.
      SID_Allreduce(SID_IN_PLACE,&flag_group_added,1,SID_INT,SID_SUM,SID.COMM_WORLD);
      if(flag_group_added==0)
         n_groups_unused++;
      else if(flag_group_added>1)
         n_groups_added_multiply++;

      // Build substructure hierarchy
      if(flag_read_sub_pointers){
         // Process each of this group's substructures
         int flag_central_set=FALSE;
         for(int j_subgroup=0;j_subgroup<n_subgroups_group;j_subgroup++){
            // Read and set substructure hierarchy pointers
            int parent_pointer=-1;
            SID_fread_all_buffer(&parent_pointer,sizeof(int),1,fp_hierarchy_in_buffer);
            tree_node_info *parent_node=NULL;
            // If a substructure is labeled as belonging to '0' or itself,
            //    then it belongs to the top substructure level.
            //    Everything else is a descendant of these.
            if(parent_pointer<=0 || parent_pointer==j_subgroup){
               parent_node     =group_node;
               flag_central_set=TRUE;
            }
            else{
               parent_node=group_substructure_pointers[parent_pointer];
               if(parent_node==NULL)
                  SID_trap_error("Failed to set a substructure pointer (rank=%d,i_subgroup=%d,parent=%d)",ERROR_LOGIC,SID.My_rank,i_subgroup+j_subgroup,parent_pointer);
            }
            // Add halo to the hierarchy
            tree_node_info *current_subgroup=group_substructure_pointers[j_subgroup];
            if(current_subgroup!=NULL){
               add_node_to_substructure_hierarchy((*trees),         // The tree datastructure
                                                  current_subgroup, // Pointer to the substructure being added
                                                  parent_node);     // Pointer to its parent
            }
         }
         // Sanity Check
         if(!flag_central_set && n_subgroups_group>0)
            SID_trap_error("Failed to identify a central substructure for a group with substructure (rank=%d,i_group=%d,i_subgroup=%d,n_subgroups_group=%d)",
                           ERROR_LOGIC,SID.My_rank,i_group,i_subgroup,n_subgroups_group);
         // Build inclusive particle counts
         for(int j_subgroup=0;j_subgroup<n_subgroups_group;j_subgroup++)
            compute_inclusive_particle_count_recursive((*trees),group_substructure_pointers[j_subgroup],NULL);
      }

      // Increment global subgroup counter
      i_subgroup+=n_subgroups_group;
    } // i_group

    // Free the buffers and perform sanity checks
    SID_fclose(&fp_trees_in);
    SID_fclose(&fp_groups_in);
    SID_fclose(&fp_subgroups_in);
    free_SID_fp_buffer(&fp_subgroups_in_buffer);
    free_SID_fp_buffer(&fp_groups_in_buffer);
    free_SID_fp_buffer(&fp_trees_in_buffer);
    if(flag_read_sub_pointers){
       SID_fclose(&fp_hierarchy_in);
       free_SID_fp_buffer(&fp_hierarchy_in_buffer);
    }

    // Update the temporary look-up arrays
    update_trees_lookup((*trees),i_file);

    // Read extended pointer set (optional)
    SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
    if(flag_read_extended_pointers){
       // Read forematch pointers
       if(i_file_bridge>=0 && i_file_bridge<n_snaps)
          read_trees_pointers((*trees),
                              filename_input_dir_horizontal_trees,
                              i_file_bridge,
                              i_read_bridge,
                              READ_TREES_POINTERS_BRIDGE_FOREMATCH);
       // Read backmatch pointers
       if(i_file_backmatch>=0 && i_file_backmatch<n_snaps)
          read_trees_pointers((*trees),
                              filename_input_dir_horizontal_trees,
                              i_file_backmatch,
                              i_read_backmatch,
                              READ_TREES_POINTERS_BRIDGE_BACKMATCH);
    }
    SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

    // Report some group statistics
    //if(n_groups_unused>0)
    //   SID_log("%d groups unused...",SID_LOG_CONTINUE,n_groups_unused);
    if(n_groups_added_multiply>0)
       SID_log("%d groups used by multiple cores...",SID_LOG_CONTINUE,n_groups_added_multiply);

    SID_log("Done.",SID_LOG_CLOSE);
  } // i_read

  // Finish reading pointers 
  if(flag_read_extended_pointers){
     SID_log("Reading extended pointers...",SID_LOG_OPEN);
     SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
     for(;i_read_backmatch>=i_read_start;i_read_backmatch-=i_read_step,i_file_backmatch--){
        read_trees_pointers((*trees),
                            filename_input_dir_horizontal_trees,
                            i_file_backmatch,
                            i_read_backmatch,
                            READ_TREES_POINTERS_BRIDGE_BACKMATCH);
     }
     SID_log("Done.",SID_LOG_CLOSE);
  }
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  SID_log("Done.",SID_LOG_CLOSE);

  // Create halo sums
  calc_sum_global(&((*trees)->n_groups_trees_local),   &((*trees)->n_groups_trees),   1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
  calc_sum_global(&((*trees)->n_subgroups_trees_local),&((*trees)->n_subgroups_trees),1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);

  // Finalize trees
  finalize_trees((*trees),read_mode);

  // Compute some useful array size maxima
  calc_max((*trees)->n_groups_snap_local,     &((*trees)->max_n_groups_snap_local),     (*trees)->n_snaps,        SID_INT,CALC_MODE_DEFAULT);
  calc_max((*trees)->n_subgroups_snap_local,  &((*trees)->max_n_subgroups_snap_local),  (*trees)->n_snaps,        SID_INT,CALC_MODE_DEFAULT);
  calc_max((*trees)->n_groups_forest_local,   &((*trees)->max_n_groups_forest_local),   (*trees)->n_forests_local,SID_INT,CALC_MODE_DEFAULT);
  calc_max((*trees)->n_subgroups_forest_local,&((*trees)->max_n_subgroups_forest_local),(*trees)->n_forests_local,SID_INT,CALC_MODE_DEFAULT);

  // Clean-up
  SID_free(SID_FARG group_substructure_pointers);
  if(!check_mode_for_flag((*trees)->mode,TREE_MODE_REFERENCE))
     free_trees_lookup((*trees));

  SID_log("Done.",SID_LOG_CLOSE);
}

