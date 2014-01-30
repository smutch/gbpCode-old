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

void read_trees(char *filename_tree_root,
                char *filename_cat_root,
                char *filename_snap_list_in,
                int   n_files_groups,
                int   n_files_subgroups,
                int   k_match){
  SID_fp      fp_in;
  SID_fp      fp_out;
  SID_fp      fp_out_MBP;
  char        filename_in[256];
  char        filename_out[256];
  char        filename_out_MBP[256];
  char        filename_output_file_root[256];
  char        filename_output_dir_horizontal[256];
  char        filename_output_dir_horizontal_trees[256];
  char        filename_output_dir_vertical[256];
  char        filename_output_dir_horizontal_groups_properties[256];
  char        filename_output_dir_horizontal_groups[256];
  char        filename_output_dir_horizontal_subgroups_properties[256];
  char        filename_output_dir_horizontal_subgroups[256];
  char        filename_output_vertical_root[256];
  char       *filename_output_dir_horizontal_properties;
  char       *line=NULL;
  size_t      line_length=0;
  int         i_bin;
  int         j_bin;
  int         n_snap;
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
  int         n_forests_local;
  int         i_rank;
  int         i_search;
  int         i_forest;
  int         i_tree;
  int         k_subgroup;
  int         group_descendant_id;
  int         subgroup_descendant_id;
  int         group_file_offset;
  int         group_file_index;
  int         forest_lo_rank;
  int         n_groups;
  int         n_subgroups;
  int         n_progenitors_max;
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
  int         group_type;
  int         group_tree_id;
  int         subgroup_id;
  int         subgroup_type;
  int         subgroup_tree_id;
  int         min_sum;
  int         min_bin;
  int         n_write;
  int        *forest_lo_file=NULL;  
  int        *forest_hi_file=NULL;
  int        *tree_count_file=NULL;
  int         n_trees_file;
  int         flag_match_subgroups;
  int        *n_halos_tree_group;
  int        *n_halos_tree_subgroup;
  int         n_halos_groups;
  int         n_halos_subgroups;
  int         n_halos_used;
  int         n_halos_target;
  int        *n_halos_forest;
  int        *forest_lo_group_file=NULL;  
  int        *forest_hi_group_file=NULL;  
  int         forest_lo_group_local;  
  int         forest_hi_group_local;  
  int        *forest_lo_group_rank=NULL;
  int        *forest_hi_group_rank=NULL;
  int        *tree_count_group_file=NULL; 
  int         tree_count_group_local;  
  int        *tree_count_group_rank=NULL;
  int        *forest_lo_subgroup_file=NULL;
  int        *forest_hi_subgroup_file=NULL;
  int         forest_lo_subgroup_local;  
  int         forest_hi_subgroup_local;  
  int        *forest_lo_subgroup_rank=NULL;
  int        *forest_hi_subgroup_rank=NULL;
  int        *tree_count_subgroup_file=NULL;
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
  int         halo_snap,descendant_snap;
  halo_MBP_info     halo_MBP;
  tree_info       **trees;
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
  write_a_list(filename_snap_list_in,
               filename_tree_root,
               i_read_start,
               i_read_stop,
               i_read_step);

  // We need i_read_start,i_read_stop,i_read_step from above before we can write this status message
  SID_log("Constructing vertical merger trees for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_start,i_read_stop,i_read_step);

  // Determine the directory that the horizontal tree files are in
  strcpy(filename_output_file_root,filename_tree_root);
  strip_path(filename_output_file_root);
  sprintf(filename_output_dir_horizontal,                     "%s/horizontal",filename_tree_root);
  sprintf(filename_output_dir_horizontal_groups,              "%s/groups",    filename_output_dir_horizontal);
  sprintf(filename_output_dir_horizontal_subgroups,           "%s/subgroups", filename_output_dir_horizontal);
  sprintf(filename_output_dir_horizontal_trees,               "%s/trees",     filename_output_dir_horizontal);
  sprintf(filename_output_dir_vertical,                       "%s/vertical",  filename_tree_root);
  sprintf(filename_output_dir_horizontal_groups_properties,   "%s/properties",filename_output_dir_horizontal_groups);
  sprintf(filename_output_dir_horizontal_subgroups_properties,"%s/properties",filename_output_dir_horizontal_subgroups);
  mkdir(filename_output_dir_vertical,02755);

  // Compute a histogram of tree occupation by tree id
  for(i_read=i_read_stop,n_snap=0;i_read>=i_read_start;i_read-=i_read_step) n_snap++;

  // Figure-out what the last-used snapshot is so we can read its header
  int i_read_last;
  for(i_read=i_read_stop;i_read>=i_read_start;i_read-=i_read_step) i_read_last=i_read;

  // Determining iso-tree id mappings (onto tree structures, ranks and files).  Read header info
  //    from the i_read_start'th file since that is where the final tree tally gets written.
  int n_step_in;
  int n_search_in;
  int n_groups_max_in;
  int n_subgroups_max_in;
  sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_output_dir_horizontal_trees,i_read_last);
  SID_fopen(filename_in,"r",&fp_in);
  SID_fread_all(&n_step_in,         sizeof(int),1,&fp_in);
  SID_fread_all(&n_search_in,       sizeof(int),1,&fp_in);
  SID_fread_all(&n_groups,          sizeof(int),1,&fp_in);
  SID_fread_all(&n_subgroups,       sizeof(int),1,&fp_in);
  SID_fread_all(&n_groups_max_in,   sizeof(int),1,&fp_in);
  SID_fread_all(&n_subgroups_max_in,sizeof(int),1,&fp_in);
  SID_fread_all(&n_trees_subgroup,  sizeof(int),1,&fp_in);
  SID_fread_all(&n_trees_group,     sizeof(int),1,&fp_in);
  SID_fclose(&fp_in);
  n_progenitors_max=MAX(n_groups_max_in,n_subgroups_max_in);
  if(n_step_in!=i_read_step) SID_trap_error("Snapshot step sizes don't match (ie. %d!=%d)",ERROR_LOGIC,n_step_in,i_read_step);

  // Force the forest construction to use all snapshots
  int n_search_forests=i_read_stop;

  // Read forest information
  int  *i_forest_group;
  int  *i_forest_subgroup;
  int  *n_halos_forest_group;
  int  *n_halos_forest_subgroup;
  int   n_trees_forest_groups_max;
  int   n_trees_forest_subgroups_max;
  read_forests(filename_tree_root,
               n_trees_group,
               n_trees_subgroup,
               &n_forests_group,
               &n_forests_subgroup,
               &i_forest_group,
               &i_forest_subgroup,
               &n_halos_forest_group,
               &n_halos_forest_subgroup,
               &n_trees_forest_groups_max,
               &n_trees_forest_subgroups_max);

  // DOMAIN DECOMPOSITION STARTS HERE

  // Assign trees to files
  SID_log("Assigning trees to files and ranks...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Initialize some arrays needed for domain decomposition
  forest_lo_group_rank    =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
  forest_hi_group_rank    =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
  tree_count_group_rank   =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
  forest_lo_group_file    =(int *)SID_malloc(sizeof(int)*n_files_groups);
  forest_hi_group_file    =(int *)SID_malloc(sizeof(int)*n_files_groups);
  tree_count_group_file   =(int *)SID_malloc(sizeof(int)*n_files_groups);
  forest_lo_subgroup_rank =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
  forest_hi_subgroup_rank =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
  tree_count_subgroup_rank=(int *)SID_malloc(sizeof(int)*SID.n_proc);  
  forest_lo_subgroup_file =(int *)SID_malloc(sizeof(int)*n_files_subgroups);  
  forest_hi_subgroup_file =(int *)SID_malloc(sizeof(int)*n_files_subgroups);  
  tree_count_subgroup_file=(int *)SID_malloc(sizeof(int)*n_files_subgroups);

  // Determine how subgroup trees will be distributed between ranks 
  for(i_rank=0,n_halos_used=0,i_tree=0;i_rank<SID.n_proc;i_rank++,i_tree++){
    forest_lo_subgroup_rank[i_rank] =i_tree;
    forest_hi_subgroup_rank[i_rank] =i_tree;
    tree_count_subgroup_rank[i_rank]=0;
    if(i_tree<n_forests_subgroup){
      n_halos_target                =(n_halos_subgroups-n_halos_used)/(SID.n_proc-i_rank);
      tree_count_subgroup_rank[i_rank]+=n_halos_forest_subgroup[i_tree];
      while(tree_count_subgroup_rank[i_rank]<n_halos_target && i_tree<(n_forests_subgroup-1)){
        i_tree++;
        forest_hi_subgroup_rank[i_rank]  =i_tree;
        tree_count_subgroup_rank[i_rank]+=n_halos_forest_subgroup[i_tree];
      }
      n_halos_used+=tree_count_subgroup_rank[i_rank];
    }
  }
  while(i_tree<(n_forests_subgroup-1)){
    i_tree++;
    forest_hi_subgroup_rank[SID.n_proc-1]  =i_tree;
    tree_count_subgroup_rank[SID.n_proc-1]+=n_halos_forest_subgroup[i_tree];
  }

  // Determine how group trees will be distributed between ranks 
  for(i_rank=0,n_halos_used=0,i_tree=0;i_rank<SID.n_proc;i_rank++,i_tree++){
    forest_lo_group_rank[i_rank]   =i_tree;
    forest_hi_group_rank[i_rank]   =i_tree;
    tree_count_group_rank[i_rank]=0;
    if(i_tree<n_forests_group){
      n_halos_target                =(n_halos_groups-n_halos_used)/(SID.n_proc-i_rank);
      tree_count_group_rank[i_rank]+=n_halos_forest_group[i_tree];
      while(tree_count_group_rank[i_rank]<n_halos_target && i_tree<(n_forests_group-1)){
        i_tree++;
        forest_hi_group_rank[i_rank]    =i_tree;
        tree_count_group_rank[i_rank]+=n_halos_forest_group[i_tree];
      }
      n_halos_used+=tree_count_group_rank[i_rank];
    }
  }
  while(i_tree<(n_forests_group-1)){
    i_tree++;
    forest_hi_group_rank[SID.n_proc-1]    =i_tree;
    tree_count_group_rank[SID.n_proc-1]+=n_halos_forest_group[i_tree];
  }

  // Determine how subgroup trees will be distributed between files
  for(i_file=0,n_halos_used=0,i_tree=0;i_file<n_files_subgroups;i_file++,i_tree++){
    forest_lo_subgroup_file[i_file]   =i_tree;
    forest_hi_subgroup_file[i_file]   =i_tree;
    tree_count_subgroup_file[i_file]=0;
    if(i_tree<n_forests_subgroup){
      n_halos_target                   =(n_halos_subgroups-n_halos_used)/(n_files_subgroups-i_file);
      tree_count_subgroup_file[i_file]+=n_halos_forest_subgroup[i_tree];
      while(tree_count_subgroup_file[i_file]<n_halos_target && i_tree<(n_forests_subgroup-1)){
        i_tree++;
        forest_hi_subgroup_file[i_file]    =i_tree;
        tree_count_subgroup_file[i_file]+=n_halos_forest_subgroup[i_tree];
      }
      n_halos_used+=tree_count_subgroup_file[i_file];
    }
  }
  while(i_tree<(n_forests_subgroup-1)){
    i_tree++;
    forest_hi_subgroup_file[n_files_subgroups-1]    =i_tree;
    tree_count_subgroup_file[n_files_subgroups-1]+=n_halos_forest_subgroup[i_tree];
  }

  // Determine how group trees will be distributed between files
  for(i_file=0,n_halos_used=0,i_tree=0;i_file<n_files_groups;i_file++,i_tree++){
    forest_lo_group_file[i_file]   =i_tree;
    forest_hi_group_file[i_file]   =i_tree;
    tree_count_group_file[i_file]=0;
    if(i_tree<n_forests_group){
      n_halos_target                =(n_halos_groups-n_halos_used)/(n_files_groups-i_file);
      tree_count_group_file[i_file]+=n_halos_forest_group[i_tree];
      while(tree_count_group_file[i_file]<n_halos_target && i_tree<(n_forests_group-1)){
        i_tree++;
        forest_hi_group_file[i_file]    =i_tree;
        tree_count_group_file[i_file]+=n_halos_forest_group[i_tree];
      }
      n_halos_used+=tree_count_group_file[i_file];
    }
  }
  while(i_tree<(n_forests_group-1)){
    i_tree++;
    forest_hi_group_file[n_files_groups-1]    =i_tree;
    tree_count_group_file[n_files_groups-1]+=n_halos_forest_group[i_tree];
  }

  // Set aliases for local values
  forest_lo_subgroup_local =forest_lo_subgroup_rank[SID.My_rank];
  forest_hi_subgroup_local =forest_hi_subgroup_rank[SID.My_rank];
  n_trees_subgroup_local   =forest_hi_subgroup_local-forest_lo_subgroup_local+1;
  tree_count_subgroup_local=tree_count_subgroup_rank[SID.My_rank];
  forest_lo_group_local    =forest_lo_group_rank[SID.My_rank];
  forest_hi_group_local    =forest_hi_group_rank[SID.My_rank];
  n_trees_group_local      =forest_hi_group_local-forest_lo_group_local+1;
  tree_count_group_local   =tree_count_group_rank[SID.My_rank];

  // Report decomposition results
  if(n_files_subgroups>1){
    for(i_bin=0;i_bin<n_files_subgroups;i_bin++)
      SID_log("Subgroup file #%4d will store %5d subgroup trees (%8d halos in total)",
               SID_LOG_COMMENT,i_bin+1,forest_hi_subgroup_file[i_bin]-forest_lo_subgroup_file[i_bin]+1,tree_count_subgroup_file[i_bin]);    
  }
  if(n_files_groups>1){
    for(i_bin=0;i_bin<n_files_groups;i_bin++)
      SID_log("Group    file #%4d will store %5d group    trees (%8d halos in total)",
              SID_LOG_COMMENT,i_bin+1,forest_hi_group_file[i_bin]-forest_lo_group_file[i_bin]+1,tree_count_group_file[i_bin]);    
  }
  if(SID.n_proc>1){
    for(i_bin=0;i_bin<SID.n_proc;i_bin++)
      SID_log("Rank     #%4d will process %5d subgroup trees (%8d halos in total)",
              SID_LOG_COMMENT,i_bin+1,forest_hi_subgroup_rank[i_bin]-forest_lo_subgroup_rank[i_bin]+1,tree_count_subgroup_rank[i_bin]);
    for(i_bin=0;i_bin<SID.n_proc;i_bin++)
      SID_log("Rank     #%4d will process %5d group    trees (%8d halos in total)",
              SID_LOG_COMMENT,i_bin+1,forest_hi_group_rank[i_bin]-forest_lo_group_rank[i_bin]+1,tree_count_group_rank[i_bin]);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Process subgroup trees (k_match==K_MATCH_SUBGROUPS) and then group trees (k_match==K_MATCH_GROUPS)
  //   If you want to change the order in which these are done, then reverse the preprocessor settings
  fp_catalog_info  fp_group_properties;
  fp_catalog_info  fp_subgroup_properties;
  halo_info       *group_properties;
  halo_info       *subgroup_properties;
  group_properties   =(halo_info *)SID_calloc(sizeof(halo_info));
  subgroup_properties=(halo_info *)SID_calloc(sizeof(halo_info));

  // Process subgroups and then groups
  switch(k_match){
    case K_MATCH_SUBGROUPS:
      filename_output_dir_horizontal_properties=filename_output_dir_horizontal_subgroups_properties;
      sprintf(filename_output_vertical_root,"%s",filename_output_dir_vertical);
      n_halos_forest   =n_halos_forest_subgroup;
      n_forests_local  =n_trees_subgroup_local;
      forest_lo_rank   =forest_lo_subgroup_rank[SID.My_rank];
      tree_count_file=tree_count_subgroup_file;
      forest_lo_file   =forest_lo_subgroup_file;
      forest_hi_file   =forest_hi_subgroup_file;
      n_write        =n_files_subgroups;
      sprintf(group_text_prefix,"sub");
      break;
    case K_MATCH_GROUPS:
      filename_output_dir_horizontal_properties=filename_output_dir_horizontal_groups_properties;
      sprintf(filename_output_vertical_root,"%s",filename_output_dir_vertical);
      n_halos_forest   =n_halos_forest_group;
      n_forests_local  =n_trees_group_local;
      forest_lo_rank   =forest_lo_group_rank[SID.My_rank];
      tree_count_file=tree_count_group_file;
      forest_lo_file   =forest_lo_group_file;
      forest_hi_file   =forest_hi_group_file;
      n_write        =n_files_groups;
      sprintf(group_text_prefix,"");
      break;      
  }
    
  // Initialize tree arrays
  trees=(tree_info **)SID_malloc(sizeof(tree_info *)*n_forests_local);
  for(i_tree=0;i_tree<n_forests_local;i_tree++)
    init_trees_vertical(n_snap,&trees[i_tree]);
  
  // Print a status message
  SID_log("Reading %sgroup trees...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix);

  // Loop over all the horizontal tree files in order of decreasing snapshot number, hanging halos on the trees as we go
  int n_subgroups_orphaned=0;
  for(i_read=i_read_stop,i_file=n_snap-1,flag_init=TRUE;i_read>=i_read_start;i_read-=i_read_step,i_file--,flag_init=FALSE){
    SID_log("Processing snapshot %03d (%03d of %03d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read,i_file+1,n_snap);

    // This counter makes sure that the halo snap index in the trees
    //   is continuous, even if we are skipping snapshots
    halo_snap=i_file;

    // Initialize this snapshot's halo list here
    for(i_tree=0;i_tree<n_forests_local;i_tree++){
      trees[i_tree]->neighbour_halos[i_file]=NULL;
      trees[i_tree]->n_neighbours[i_file]   =0;
    }
    sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_output_dir_horizontal_trees,i_read);

    // Open properties catalog
    fopen_catalog(filename_cat_root,
                  i_read,
                  READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES,
                  &fp_group_properties);
    if(k_match==K_MATCH_SUBGROUPS)
       fopen_catalog(filename_cat_root,
                     i_read,
                     READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES,
                     &fp_subgroup_properties);

    // Open horizontal tree file
    SID_fopen(filename_in,"r",&fp_in);

    // Read header
    SID_fread_all(&n_step_in,         sizeof(int),1,&fp_in);
    SID_fread_all(&n_search_in,       sizeof(int),1,&fp_in);
    SID_fread_all(&n_groups,          sizeof(int),1,&fp_in);
    SID_fread_all(&n_subgroups,       sizeof(int),1,&fp_in);
    SID_fread_all(&n_groups_max_in,   sizeof(int),1,&fp_in);
    SID_fread_all(&n_subgroups_max_in,sizeof(int),1,&fp_in);
    SID_fread_all(&n_trees_subgroup,  sizeof(int),1,&fp_in);
    SID_fread_all(&n_trees_group,     sizeof(int),1,&fp_in);
    n_progenitors_max=MAX(n_groups_max_in,n_subgroups_max_in);
    if(n_step_in!=i_read_step) 
       SID_trap_error("Snapshot step sizes don't match (ie. %d!=%d)",ERROR_LOGIC,n_step_in,i_read_step);
    
    // Read each group in turn
    int tree_read_buffer[7];
    for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){

      // Read horizontal trees for groups
      SID_fread_all(tree_read_buffer,7*sizeof(int),1,&fp_in);
      group_id           =tree_read_buffer[0];
      group_type         =tree_read_buffer[1];
      group_descendant_id=tree_read_buffer[2];
      group_tree_id      =tree_read_buffer[3];
      group_file_offset  =tree_read_buffer[4];
      group_file_index   =tree_read_buffer[5];
      n_subgroups_group  =tree_read_buffer[6];

      if(group_tree_id>=0)
        group_tree_id=i_forest_group[group_tree_id];
      else
        group_tree_id=-1;

      // Read halo information from catalog files (needed even for subgroup trees; we need FoF masses for the most massive progenitors)
      fread_catalog_file(&fp_group_properties,group_properties,NULL,NULL,i_group);

      // If we are processing subgroup trees ...
      if(k_match==K_MATCH_SUBGROUPS){
        // Read each subgroup in turn
        int i_subgroup_valid;
        for(j_subgroup=0,i_subgroup_valid=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
          //SID_fread_all(&(subgroup_id),           sizeof(int),1,&fp_in);
          //SID_fread_all(&(subgroup_type),         sizeof(int),1,&fp_in);
          //SID_fread_all(&(subgroup_descendant_id),sizeof(int),1,&fp_in);
          //SID_fread_all(&(subgroup_tree_id),      sizeof(int),1,&fp_in);
          //SID_fread_all(&(subgroup_file_offset),  sizeof(int),1,&fp_in);
          //SID_fread_all(&(subgroup_file_index),   sizeof(int),1,&fp_in);
          SID_fread_all(tree_read_buffer,6*sizeof(int),1,&fp_in);
          subgroup_id           =tree_read_buffer[0];
          subgroup_type         =tree_read_buffer[1];
          subgroup_descendant_id=tree_read_buffer[2];
          subgroup_tree_id      =tree_read_buffer[3];
          subgroup_file_offset  =tree_read_buffer[4];
          subgroup_file_index   =tree_read_buffer[5];

          // Ignore negative ids
          if(subgroup_tree_id>=0)
            subgroup_tree_id=i_forest_subgroup[subgroup_tree_id];
          else
            subgroup_tree_id=-1;
          if(subgroup_id>=0 && subgroup_tree_id>=0){
            // If this subgroup belongs to a local tree ...
            i_tree=subgroup_tree_id-forest_lo_rank;
            if(i_tree>=0 && i_tree<n_forests_local){ 
              // ... create a new branch and add it to its tree ...
              if(subgroup_file_offset<=0)
                descendant_snap=-1; // Needed for halos in the root snapshot
              else
                descendant_snap=(i_file+subgroup_file_offset);
              // ... read halo information from catalog files ...
              fread_catalog_file(&fp_subgroup_properties,subgroup_properties,NULL,NULL,i_subgroup);
              // ... set the most massive progenitor's mass to the FoF mass ...
              if(j_subgroup==0)
                 subgroup_properties->M_vir=group_properties->M_vir;
              // ... adjust snap counter to make things work with skipped snaps ...
              subgroup_properties->snap_num=halo_snap;
              // ... add this halo to the trees ...
              int flag_found_group=
              add_node_to_tree(trees[i_tree],
                               subgroup_type,
                               subgroup_id,
                               group_id,
                               subgroup_descendant_id,
                               halo_snap,
                               descendant_snap,
                               subgroup_properties);
              if(i_subgroup_valid>0 && !flag_found_group)
                 n_subgroups_orphaned++;
              i_subgroup_valid++;
            }
          }
        }
      }
      // ... else, process group trees 
      else{
        SID_fseek(&fp_in,sizeof(int),n_subgroups_group*6,SID_SEEK_CUR);
        // Ignore negative IDs
        if(group_id>=0 && group_tree_id>=0){
          i_tree=group_tree_id-forest_lo_rank;
          // If this group belongs to a local tree ...
          if(i_tree>=0 && i_tree<n_forests_local){ 
            // ... if so, create a new branch and add it to the tree ...
            if(group_file_offset<=0)
              descendant_snap=-1; // Needed for halos in the root snapshot
            else
              descendant_snap=(i_file+group_file_offset);
            // ... adjust snap counter to make things work with skipped snaps ...
            group_properties->snap_num=halo_snap;
            // ... add this halo to the trees ...
            add_node_to_tree(trees[i_tree],
                             group_type,
                             group_id,
                             group_id,
                             group_descendant_id,
                             halo_snap,
                             descendant_snap,                      
                             group_properties);
          }
        }
      } // if k_match
    } // i_group
    SID_fclose(&fp_in);
    fclose_catalog(&fp_group_properties);
    if(k_match==K_MATCH_SUBGROUPS)
       fclose_catalog(&fp_subgroup_properties);

    SID_log("Done.",SID_LOG_CLOSE);
  } // i_read
  SID_log("Done.",SID_LOG_CLOSE);

  // Report any subgroups that become orphaned from their parent halo
  if(k_match==K_MATCH_SUBGROUPS){
     SID_Allreduce(SID_IN_PLACE,&n_subgroups_orphaned,1,SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_log("No. of subgroups that failed group assignment = %d",SID_LOG_COMMENT,n_subgroups_orphaned);
  }

  // Finalize trees
  finalize_trees_vertical(trees,&(n_halos_forest[forest_lo_rank]),n_forests_local,n_snap,TREE_PROGENITOR_ORDER_DELUCIA);

  // Free trees
  for(i_tree=0;i_tree<n_forests_local;i_tree++)
    free_trees_vertical(&trees[i_tree]);
  SID_free((void **)&trees);

  SID_free(SID_FARG group_properties);
  SID_free(SID_FARG subgroup_properties);
  
  // Clean-up
  SID_free(SID_FARG line);
  SID_free(SID_FARG i_forest_group);
  SID_free(SID_FARG i_forest_subgroup);
  SID_free(SID_FARG n_halos_forest_group);
  SID_free(SID_FARG n_halos_forest_subgroup);
  SID_free(SID_FARG forest_lo_group_rank);
  SID_free(SID_FARG forest_hi_group_rank);
  SID_free(SID_FARG tree_count_group_rank);
  SID_free(SID_FARG forest_lo_subgroup_rank);
  SID_free(SID_FARG forest_hi_subgroup_rank);
  SID_free(SID_FARG tree_count_subgroup_rank);
  SID_free(SID_FARG forest_lo_group_file);
  SID_free(SID_FARG forest_hi_group_file);
  SID_free(SID_FARG tree_count_group_file);
  SID_free(SID_FARG forest_lo_subgroup_file);
  SID_free(SID_FARG forest_hi_subgroup_file);
  SID_free(SID_FARG tree_count_subgroup_file);

  SID_log("Done.",SID_LOG_CLOSE);
}

