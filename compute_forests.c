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

void compute_forests(char *filename_root_out,int n_search_forests){
  SID_fp      fp_in;
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
  int         k_match;
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
  int  flag_match_subgroups;
  int  group_forest_array_min;
  int  subgroup_forest_array_min;
  int *n_halos_tree_group;
  int *n_halos_tree_subgroup;
  int  n_halos_groups;
  int  n_halos_subgroups;
  int  n_halos_used;
  int  n_halos_target;
  int *n_halos_forest;
  int *forest_lo_group_file=NULL;  
  int *forest_hi_group_file=NULL;  
  int  forest_lo_group_local;  
  int  forest_hi_group_local;  
  int *forest_lo_group_rank=NULL;
  int *forest_hi_group_rank=NULL;
  int *tree_count_group_file=NULL; 
  int  tree_count_group_local;  
  int *tree_count_group_rank=NULL;
  int *forest_lo_subgroup_file=NULL;
  int *forest_hi_subgroup_file=NULL;
  int  forest_lo_subgroup_local;  
  int  forest_hi_subgroup_local;  
  int *forest_lo_subgroup_rank=NULL;
  int *forest_hi_subgroup_rank=NULL;
  int *tree_count_subgroup_file=NULL;
  int  tree_count_subgroup_local;
  int *tree_count_subgroup_rank=NULL;
  int  n_trees_subgroup_local;
  int  n_trees_group_local;
  int  subgroup_file_offset;
  int  subgroup_file_index;
  char group_text_prefix[4];
  int  n_conjoined;
  int  n_conjoined_total;
  halo_MBP_info     halo_MBP;
  tree_info       **trees;
  tree_node_info   *current=NULL;
  tree_node_info   *last   =NULL;
  tree_node_info   *next   =NULL;
  int               depth_first_index;
  int flag_write_init;
  int k_tree;
  int flag_init;
  int progenitor_score;
  int flag;
  int max_id=0;
  int n_halos_written;
  int halo_snap,descendant_snap;

  if(SID.I_am_Master){
     // Read the tree search/scan parameters
     int i_read_start;
     int i_read_stop;
     int i_read_step;
     int n_search;
     int flag_fix_bridges;
     int flag_compute_fragmented;
     int flag_compute_ghosts;
     read_tree_run_parameters(filename_root_out,
                              &i_read_start,
                              &i_read_stop,
                              &i_read_step,
                              &n_search,
                              &flag_fix_bridges,
                              &flag_compute_fragmented,
                              &flag_compute_ghosts);

     // Compute a histogram of tree occupation by tree id
     for(i_read=i_read_stop,n_snap=0;i_read>=i_read_start;i_read-=i_read_step) n_snap++;

     // Determine the directory that the horizontal tree files are in
     strcpy(filename_output_file_root,filename_root_out);
     strip_path(filename_output_file_root);
     sprintf(filename_output_dir_horizontal,                     "%s/horizontal",filename_root_out);
     sprintf(filename_output_dir_horizontal_groups,              "%s/groups",    filename_output_dir_horizontal);
     sprintf(filename_output_dir_horizontal_subgroups,           "%s/subgroups", filename_output_dir_horizontal);
     sprintf(filename_output_dir_horizontal_trees,               "%s/trees",     filename_output_dir_horizontal);
     sprintf(filename_output_dir_vertical,                       "%s/vertical",  filename_root_out);
     sprintf(filename_output_dir_horizontal_groups_properties,   "%s/properties",filename_output_dir_horizontal_groups);
     sprintf(filename_output_dir_horizontal_subgroups_properties,"%s/properties",filename_output_dir_horizontal_subgroups);
     mkdir(filename_output_dir_vertical,02755);

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
     if(n_step_in!=i_read_step) SID_trap_error("Snapshot step sizes don't match (ie. %d!=%d)",ERROR_LOGIC,n_step_in,i_read_step);

     // Determine the mapping between the horizontal tree IDs and the
     //   forest IDs that indicate how the isotrees are packaged into
     //   vertical trees.
     // Read-in all the tree halos and join trees into forests
     //    when substructure tree IDs start mixing between their initial values.
     //    Perform scan only over n_search to prevent over-linking.

     // Initialize forest arrays
     int *group_forest_array;
     int *subgroup_forest_array;
     group_forest_array   =(int  *)SID_malloc(sizeof(int)*n_trees_group);
     subgroup_forest_array=(int  *)SID_malloc(sizeof(int)*n_trees_subgroup);
     n_halos_tree_group   =(int  *)SID_malloc(sizeof(int)*n_trees_group);
     n_halos_tree_subgroup=(int  *)SID_malloc(sizeof(int)*n_trees_subgroup);
     for(i_tree=0;i_tree<n_trees_group;i_tree++){
       group_forest_array[i_tree] = i_tree; // Initially, every group tree is given it's own forest
       n_halos_tree_group[i_tree] = 0;
     }
     for(i_tree=0;i_tree<n_trees_subgroup;i_tree++){
       subgroup_forest_array[i_tree]=-1;    // This will be set later to point to a group tree id
       n_halos_tree_subgroup[i_tree]= 0;
     }

     // Scan snapshots, forming forests as we go
     int n_trees_group_i;
     int n_trees_subgroup_i;
     int n_halos_groups_unused=0;
     int n_halos_subgroups_unused=0;
     int tree_read_buffer[7];
     SID_log("Generating mapping of trees to forests...",SID_LOG_OPEN|SID_LOG_TIMER);
     for(i_read=i_read_stop,n_halos_groups=0,n_halos_subgroups=0,j_read=0;
         i_read>=i_read_start;
         i_read-=i_read_step,j_read+=i_read_step){
       SID_log("Processing snapshot #%03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read);
       if(!(j_read<=n_search_forests))
          SID_log("(forest joining is off)...",SID_LOG_CONTINUE);
       else
          SID_log("(forest joining is on)...",SID_LOG_CONTINUE);
       sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_output_dir_horizontal_trees,i_read);
       SID_fopen(filename_in,"r",&fp_in);
       SID_fread_all(&n_step_in,         sizeof(int),1,&fp_in);
       SID_fread_all(&n_search_in,       sizeof(int),1,&fp_in);
       SID_fread_all(&n_groups,          sizeof(int),1,&fp_in);
       SID_fread_all(&n_subgroups,       sizeof(int),1,&fp_in);
       SID_fread_all(&n_groups_max_in,   sizeof(int),1,&fp_in);
       SID_fread_all(&n_subgroups_max_in,sizeof(int),1,&fp_in);
       SID_fread_all(&n_trees_subgroup_i,sizeof(int),1,&fp_in);
       SID_fread_all(&n_trees_group_i,   sizeof(int),1,&fp_in);
       if(n_step_in!=i_read_step) SID_trap_error("Snapshot step sizes don't match (ie. %d!=%d)",ERROR_LOGIC,n_step_in,i_read_step);
       // Loop over groups
       for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
          // Read group
          SID_fread_all(tree_read_buffer,7*sizeof(int),1,&fp_in);
          group_id           =tree_read_buffer[0];
          group_type         =tree_read_buffer[1];
          group_descendant_id=tree_read_buffer[2];
          group_tree_id      =tree_read_buffer[3];
          group_file_offset  =tree_read_buffer[4];
          group_file_index   =tree_read_buffer[5];
          n_subgroups_group  =tree_read_buffer[6];

          // Count groups and decide if this group is valid for forest building
          int flag_valid_group=TRUE;
          if(group_tree_id>=0){
            if(group_id>=0){
              n_halos_tree_group[group_tree_id]++;
              n_halos_groups++;
            }
            else
              n_halos_groups_unused++;
          }
          else{
            flag_valid_group=FALSE;
            n_halos_groups_unused++;
          }

          // Loop over subgroups
          for(j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
            // Read subgroup
            SID_fread_all(tree_read_buffer,6*sizeof(int),1,&fp_in);
            subgroup_id           =tree_read_buffer[0];
            subgroup_type         =tree_read_buffer[1];
            subgroup_descendant_id=tree_read_buffer[2];
            subgroup_tree_id      =tree_read_buffer[3];
            subgroup_file_offset  =tree_read_buffer[4];
            subgroup_file_index   =tree_read_buffer[5];
            if(subgroup_id>=0 && subgroup_tree_id>=0){
              n_halos_tree_subgroup[subgroup_tree_id]++;
              n_halos_subgroups++;
              if(flag_valid_group){
                 int subgroup_forest_array_i=subgroup_forest_array[subgroup_tree_id];
                 if(subgroup_forest_array_i<0){
                    subgroup_forest_array[subgroup_tree_id]=group_tree_id;
                    subgroup_forest_array_i                =group_tree_id;
                 }
                 // If we have tied this subgroup to a legit group_tree_id (do this only for the first 
                 //    n_search_forests snapshots; we still have to sacn over all snapshots though, 
                 //    to make sure that all subgroups are initialized) ...
                 if(subgroup_forest_array_i>=0 && j_read<=n_search_forests){
                    // ... then join group trees into forests if this subgroup demands it.

                    // To make things robust against multiple group tree linking events, we have to make sure we're
                    //    always working with the lowest-level link in the chain linking things together into forests.
                    //    Only in this way will previously linked trees be linked consistantly during subsequent links.
                    while(group_forest_array[group_tree_id]          !=group_tree_id)           group_tree_id          =group_forest_array[group_tree_id];
                    while(group_forest_array[subgroup_forest_array_i]!=subgroup_forest_array_i) subgroup_forest_array_i=group_forest_array[subgroup_forest_array_i];

                    // Perform any necessary linking.
                    if(subgroup_forest_array_i<group_tree_id)
                       group_forest_array[group_tree_id]=subgroup_forest_array_i;
                    else if(subgroup_forest_array_i>group_tree_id)
                       group_forest_array[subgroup_forest_array_i]=group_tree_id;
                 }
              }
            }
            else
               n_halos_subgroups_unused++;
          }
       }
       SID_fclose(&fp_in);
       SID_log("Done.",SID_LOG_CLOSE);
     }

     // Due to linking, there may be gaps in the forrest arrays at this point. Collapse them to fix this.
     SID_log("Collapsing results and performing halo counts...",SID_LOG_OPEN);

     // If an isotree has been joined to another to form a forest, one of them will
     //    have a reduced group_forest_array value (subgroup_forest_arrays point to this array at
     //    at this point, except for undefined subgroup trees which have 
     //    values > n_trees_group).  Follow the trail to the lowest value so that
     //    all isotrees meant for a single forrest have the same value.  This can be done
     //    by starting at the lowest tree ID and propagating results upwards.
     for(i_tree=0;i_tree<n_trees_group;i_tree++){
        if(group_forest_array[i_tree]>0)
           group_forest_array[i_tree]=group_forest_array[group_forest_array[i_tree]];
     }

     //  Eliminate gaps in the list of group forest IDs and perform halo/tree counts
     size_t *group_forest_array_index;
     int    *n_trees_forest_groups=(int *)SID_calloc(sizeof(int)*n_trees_group);
     int    *n_halos_forest_groups=(int *)SID_calloc(sizeof(int)*n_trees_group);
     int     n_trees_forest_groups_max=0;
     int     n_halos_forest_groups_max=0;
     merge_sort(group_forest_array,(size_t)n_trees_group,&group_forest_array_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
     i_tree=0;
     while(group_forest_array[group_forest_array_index[i_tree]]<0 && i_tree<(n_trees_group-1)){
        group_forest_array[group_forest_array_index[i_tree]]=-1;
        i_tree++;
     }
     if(group_forest_array[group_forest_array_index[i_tree]]<0){
        group_forest_array[group_forest_array_index[i_tree]]=-1;
        i_tree++;
     }
     for(i_forest=0;i_tree<n_trees_group;i_forest++){
       int j_forest=group_forest_array[group_forest_array_index[i_tree]];
       do{
         group_forest_array[group_forest_array_index[i_tree]]=i_forest;
         n_trees_forest_groups[i_forest]++;
         n_halos_forest_groups[i_forest]+=n_halos_tree_group[group_forest_array_index[i_tree]];
         i_tree++;
         if(i_tree>=n_trees_group) break;
       }
       while(group_forest_array[group_forest_array_index[i_tree]]==j_forest);
       if(n_trees_forest_groups[i_forest]<=0)
          SID_trap_error("A group forest (%d) has been assigned no trees.",ERROR_LOGIC,i_forest);
     }
     n_forests_group=i_forest+1;
     calc_max(n_trees_forest_groups,&n_trees_forest_groups_max,n_forests_group,SID_INT,CALC_MODE_DEFAULT);
     calc_max(n_halos_forest_groups,&n_halos_forest_groups_max,n_forests_group,SID_INT,CALC_MODE_DEFAULT);
     SID_free(SID_FARG group_forest_array_index);

     // At this point, the subgroups just point to group IDs.  Give each
     //   subgroup to that group ID's forest. 
     for(i_tree=0;i_tree<n_trees_subgroup;i_tree++){
        group_tree_id=subgroup_forest_array[i_tree];
        if(group_tree_id>=0 && group_tree_id<n_trees_group)
           subgroup_forest_array[i_tree]=group_forest_array[group_tree_id];
        else
           subgroup_forest_array[i_tree]=-1;
     }

     //  Eliminate gaps in the list of forest IDs and perform halo/tree counts
     size_t *subgroup_forest_array_index;
     int    *n_trees_forest_subgroups=(int *)SID_calloc(sizeof(int)*n_trees_subgroup);
     int    *n_halos_forest_subgroups=(int *)SID_calloc(sizeof(int)*n_trees_subgroup);
     int     n_trees_forest_subgroups_max=0;
     int     n_halos_forest_subgroups_max=0;
     merge_sort(subgroup_forest_array,(size_t)n_trees_subgroup,&subgroup_forest_array_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
     i_tree=0;
     while(subgroup_forest_array[subgroup_forest_array_index[i_tree]]<0 && i_tree<(n_trees_subgroup-1)){
        subgroup_forest_array[subgroup_forest_array_index[i_tree]]=-1;
        i_tree++;
     }
     if(subgroup_forest_array[subgroup_forest_array_index[i_tree]]<0){
        subgroup_forest_array[subgroup_forest_array_index[i_tree]]=-1;
        i_tree++;
     }
     for(i_forest=0;i_tree<n_trees_subgroup;i_forest++){
       int j_forest=subgroup_forest_array[subgroup_forest_array_index[i_tree]];
       do{
         subgroup_forest_array[subgroup_forest_array_index[i_tree]]=i_forest;
         n_trees_forest_subgroups[i_forest]++;
         n_halos_forest_subgroups[i_forest]+=n_halos_tree_subgroup[subgroup_forest_array_index[i_tree]];
         i_tree++;
         if(i_tree>=n_trees_subgroup) break;
       }
       while(subgroup_forest_array[subgroup_forest_array_index[i_tree]]==j_forest);
       if(n_trees_forest_subgroups[i_forest]<=0)
          SID_trap_error("A subgroup forest (%d) has been assigned no trees.",ERROR_LOGIC,i_forest);
     }
     n_forests_subgroup=i_forest+1;
     calc_max(n_trees_forest_subgroups,&n_trees_forest_subgroups_max,n_forests_subgroup,SID_INT,CALC_MODE_DEFAULT);
     calc_max(n_halos_forest_subgroups,&n_halos_forest_subgroups_max,n_forests_subgroup,SID_INT,CALC_MODE_DEFAULT);
     SID_free(SID_FARG subgroup_forest_array_index);
     
     SID_log("Done.",SID_LOG_CLOSE);

     // Write results for groups
     SID_log("Writing results...",SID_LOG_OPEN);
     FILE *fp_out=NULL;
     int   i_column;
     sprintf(filename_out,"%s/tree2forest_mapping_groups.txt",filename_root_out);
     fp_out=fopen(filename_out,"w");
     i_column=1;
     fprintf(fp_out,"# Tree->forest mappings for groups.\n");
     fprintf(fp_out,"# No. of groups                         = %d\n",n_halos_groups);
     fprintf(fp_out,"# No. of groups unused                  = %d\n",n_halos_groups_unused);
     fprintf(fp_out,"# No. of group isotrees                 = %d\n",n_trees_group);
     fprintf(fp_out,"# No. of group forests                  = %d\n",n_forests_group);
     fprintf(fp_out,"# Max. No. of trees in one group forest = %d\n",n_trees_forest_groups_max);
     fprintf(fp_out,"# Max. No. of halos in one group forest = %d\n",n_halos_forest_groups_max);
     fprintf(fp_out,"# No. of group forests                  = %d\n",n_forests_group);
     fprintf(fp_out,"# Column (%02d): Horizontal tree ID\n",  i_column++);
     fprintf(fp_out,"#        (%02d): Forest ID\n",           i_column++);
     fprintf(fp_out,"#        (%02d): No. of halos in tree\n",i_column++);
     for(i_tree=0;i_tree<n_trees_group;i_tree++)
        fprintf(fp_out,"%d %d %d\n",i_tree,group_forest_array[i_tree],n_halos_tree_group[i_tree]);
     fclose(fp_out);
     sprintf(filename_out,"%s/forest_info_groups.txt",filename_root_out);
     fp_out=fopen(filename_out,"w");
     i_column=1;
     fprintf(fp_out,"# Info about group forests.\n");
     fprintf(fp_out,"# No. of groups                         = %d\n",n_halos_groups);
     fprintf(fp_out,"# No. of groups unused                  = %d\n",n_halos_groups_unused);
     fprintf(fp_out,"# No. of group isotrees                 = %d\n",n_trees_group);
     fprintf(fp_out,"# No. of group forests                  = %d\n",n_forests_group);
     fprintf(fp_out,"# Max. No. of trees in one group forest = %d\n",n_trees_forest_groups_max);
     fprintf(fp_out,"# Max. No. of halos in one group forest = %d\n",n_halos_forest_groups_max);
     fprintf(fp_out,"# No. of group forests                  = %d\n",n_forests_group);
     fprintf(fp_out,"# Column (%02d): Forest ID\n",             i_column++);
     fprintf(fp_out,"#        (%02d): No. of trees in forest\n",i_column++);
     fprintf(fp_out,"#        (%02d): No. of halos in forest\n",i_column++);
     for(i_forest=0;i_forest<n_forests_group;i_forest++)
        fprintf(fp_out,"%d %d %d\n",i_forest,n_trees_forest_groups[i_forest],n_halos_forest_groups[i_forest]);
     fclose(fp_out);

     // Write results for subgroups
     sprintf(filename_out,"%s/tree2forest_mapping_subgroups.txt",filename_root_out);
     fp_out=fopen(filename_out,"w");
     i_column=1;
     fprintf(fp_out,"# Tree->forest mappings for subgroups.\n");
     fprintf(fp_out,"# No. of subgroups                         = %d\n",n_halos_subgroups);
     fprintf(fp_out,"# No. of subgroups unused                  = %d\n",n_halos_subgroups_unused);
     fprintf(fp_out,"# No. of subgroup isotrees                 = %d\n",n_trees_subgroup);
     fprintf(fp_out,"# No. of subgroup forests                  = %d\n",n_forests_subgroup);
     fprintf(fp_out,"# Max. No. of trees in one subgroup forest = %d\n",n_trees_forest_subgroups_max);
     fprintf(fp_out,"# Max. No. of halos in one subgroup forest = %d\n",n_halos_forest_subgroups_max);
     fprintf(fp_out,"# No. of subgroup forests                  = %d\n",n_forests_subgroup);
     fprintf(fp_out,"# Column (%02d): Horizontal tree ID\n",  i_column++);
     fprintf(fp_out,"#        (%02d): Forest ID\n",           i_column++);
     fprintf(fp_out,"#        (%02d): No. of halos in tree\n",i_column++);
     for(i_tree=0;i_tree<n_trees_subgroup;i_tree++)
        fprintf(fp_out,"%d %d %d\n",i_tree,subgroup_forest_array[i_tree],n_halos_tree_subgroup[i_tree]);
     fclose(fp_out);
     sprintf(filename_out,"%s/forest_info_subgroups.txt",filename_root_out);
     fp_out=fopen(filename_out,"w");
     i_column=1;
     fprintf(fp_out,"# Info about subgroup forests.\n");
     fprintf(fp_out,"# No. of subgroups                         = %d\n",n_halos_subgroups);
     fprintf(fp_out,"# No. of subgroups unused                  = %d\n",n_halos_subgroups_unused);
     fprintf(fp_out,"# No. of subgroup isotrees                 = %d\n",n_trees_subgroup);
     fprintf(fp_out,"# No. of subgroup forests                  = %d\n",n_forests_subgroup);
     fprintf(fp_out,"# Max. No. of trees in one subgroup forest = %d\n",n_trees_forest_subgroups_max);
     fprintf(fp_out,"# Max. No. of halos in one subgroup forest = %d\n",n_halos_forest_subgroups_max);
     fprintf(fp_out,"# No. of subgroup forests                  = %d\n",n_forests_subgroup);
     fprintf(fp_out,"# Column (%02d): Forest ID\n",             i_column++);
     fprintf(fp_out,"#        (%02d): No. of trees in forest\n",i_column++);
     fprintf(fp_out,"#        (%02d): No. of halos in forest\n",i_column++);
     for(i_forest=0;i_forest<n_forests_subgroup;i_forest++)
        fprintf(fp_out,"%d %d %d\n",i_forest,n_trees_forest_subgroups[i_forest],n_halos_forest_subgroups[i_forest]);
     fclose(fp_out);

     SID_log("Done.",SID_LOG_CLOSE);

     // Report some statistics
     SID_log("No. of groups            = %d",SID_LOG_COMMENT,n_halos_groups);
     SID_log("No. of groups unused     = %d",SID_LOG_COMMENT,n_halos_groups_unused);
     SID_log("No. of group isotrees    = %d",SID_LOG_COMMENT,n_trees_group);
     SID_log("No. of group forests     = %d",SID_LOG_COMMENT,n_forests_group);
     SID_log("No. of subgroups         = %d",SID_LOG_COMMENT,n_halos_subgroups);
     SID_log("No. of subgroups unused  = %d",SID_LOG_COMMENT,n_halos_subgroups_unused);
     SID_log("No. of subgroup isotrees = %d",SID_LOG_COMMENT,n_trees_subgroup);
     SID_log("No. of subgroup forests  = %d",SID_LOG_COMMENT,n_forests_subgroup);

     // Clean-up
     SID_free(SID_FARG group_forest_array);
     SID_free(SID_FARG subgroup_forest_array);
     SID_free(SID_FARG n_trees_forest_groups);
     SID_free(SID_FARG n_trees_forest_subgroups);
     SID_free(SID_FARG n_halos_tree_group);
     SID_free(SID_FARG n_halos_tree_subgroup);
     SID_free(SID_FARG n_halos_forest_groups);
     SID_free(SID_FARG n_halos_forest_subgroups);

     SID_log("Done.",SID_LOG_CLOSE);
  } // master rank only
}

