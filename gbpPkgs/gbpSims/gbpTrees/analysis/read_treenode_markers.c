#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

void read_treenode_markers(tree_info *trees,const char *filename_input_root,int mode){
   SID_log("Reading markers...",SID_LOG_OPEN|SID_LOG_TIMER);

SID_trap_error("This function is not working yet.  It needs to be debugged.",ERROR_LOGIC);

   // Generate the markers starting recursively from each tree root
   char filename_input_group_text[16];
   int  flag_process_groups=FALSE;
   tree_markers_info ***markers;
   int                 *n_halos_total;
   int                 *n_halos_local;
   tree_node_info     **first_neighbour;
   if(check_mode_for_flag(mode,PRECOMPUTE_TREENODE_MARKER_GROUPS)){
      sprintf(filename_input_group_text,"groups");
      flag_process_groups=TRUE;
      markers        =&(trees->group_markers);
      n_halos_total  =trees->n_groups_catalog;
      n_halos_local  =trees->n_groups_snap_local;
      first_neighbour=trees->first_neighbour_groups;
   }
   else if(check_mode_for_flag(mode,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS)){
      sprintf(filename_input_group_text,"subgroups");
      markers      =&(trees->subgroup_markers);
      n_halos_total=trees->n_subgroups_catalog;
      n_halos_local=trees->n_subgroups_snap_local;
      first_neighbour=trees->first_neighbour_subgroups;
   }
   else
      SID_trap_error("group/subgroup mode has not been properly specified in read_treenode_markers().",ERROR_LOGIC);

   // Allocate memory for the markers
   SID_log("Creating look-up arrays...",SID_LOG_OPEN|SID_LOG_TIMER);
   (*markers)                              =(tree_markers_info **)SID_malloc(sizeof(tree_markers_info *)*trees->n_snaps);
   tree_node_info ***nodes_local           =(tree_node_info   ***)SID_malloc(sizeof(tree_node_info   **)*trees->n_snaps);
   int             **file_index_local      =(int               **)SID_malloc(sizeof(int               *)*trees->n_snaps);
   size_t          **file_index_local_index=(size_t            **)SID_malloc(sizeof(size_t            *)*trees->n_snaps);
   for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
      // Allocate RAM for this snapshot
      (*markers)[i_snap]      =(tree_markers_info *)SID_malloc(sizeof(tree_markers_info)*n_halos_local[i_snap]);
      nodes_local[i_snap]     =(tree_node_info   **)SID_malloc(sizeof(tree_node_info *) *n_halos_local[i_snap]);
      file_index_local[i_snap]=(int               *)SID_malloc(sizeof(int)              *n_halos_local[i_snap]);
      // Create file index and treenode look-up arrays
      tree_node_info *current_node=first_neighbour[i_snap];
      int             i_halo      =0;
      while(current_node!=NULL){
         nodes_local[i_snap][i_halo]     =current_node;
         file_index_local[i_snap][i_halo]=current_node->file_index;
         current_node                    =current_node->next_neighbour;
         i_halo++;
      }
      // Sort file indices
      merge_sort(file_index_local[i_snap],(size_t)(n_halos_local[i_snap]),&(file_index_local_index[i_snap]),SID_INT,SORT_COMPUTE_INDEX,FALSE);
   }
   SID_log("Done.",SID_LOG_CLOSE);

   // Perform the read
   char filename_in_dir[MAX_FILENAME_LENGTH];
   sprintf(filename_in_dir,"%s_markers",filename_input_root);
   for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
      SID_log("Processing snapshot #%03d...",SID_LOG_OPEN|SID_LOG_TIMER,trees->snap_list[i_snap]);
      SID_log("(%d halos)...",SID_LOG_CONTINUE,n_halos_total[i_snap]);

      // Open file and read header
      int  i_snap_file;
      int  n_halos_total_file;
      char filename_in[MAX_FILENAME_LENGTH];
      sprintf(filename_in,"%s/%s_%03d.dat",filename_in_dir,filename_input_group_text,trees->snap_list[i_snap]);
      SID_fp fp_in;
      SID_fopen(filename_in,"r",&fp_in);
      SID_fread_all(&i_snap_file,       sizeof(int),1,&fp_in);
      SID_fread_all(&n_halos_total_file,sizeof(int),1,&fp_in);
      if(i_snap_file!=i_snap)
         SID_trap_error("Snapshot numbers don't match what's in the file (ie. %d!=%d).",ERROR_LOGIC,i_snap_file,i_snap);
      if(n_halos_total_file!=n_halos_total[i_snap])
         SID_trap_error("Halo counts don't match what's in the file (ie. %d!=%d).",ERROR_LOGIC,n_halos_total_file,n_halos_total);

      // Perform the (buffered) read of this file
      int            i_found         =0;
      int            n_markers       =5;
      int            bytes_per_marker=4;
      int            bytes_per_halo  =n_markers*bytes_per_marker;
      SID_fp_buffer *fp_in_buffer    =NULL;
      init_SID_fp_buffer(&fp_in,(size_t)(bytes_per_halo*n_halos_total[i_snap])*sizeof(int),SIZE_OF_MEGABYTE,&fp_in_buffer);
      for(int i_halo=0;i_halo<n_halos_total[i_snap];i_halo++){
         int                flag_keep    =FALSE;
         size_t             i_found_index=file_index_local_index[i_snap][i_found];
         tree_markers_info *markers_i    =&((*markers)[i_snap][i_found_index]);
         if(file_index_local[i_snap][i_found_index]==i_halo)
            flag_keep=TRUE;
         for(int i_marker=0;i_marker<n_markers;i_marker++){
            int halo_snap; SID_fread_all_buffer(&halo_snap, sizeof(int),1,fp_in_buffer);
            int halo_index;SID_fread_all_buffer(&halo_index,sizeof(int),1,fp_in_buffer);
            if(flag_keep){
               tree_node_info **marker;
               switch(i_marker){
                  case 0:
                     marker=&(markers_i->branch_leaf);
                     break;
                  case 1:
                     marker=&(markers_i->branch_root);
                     break;
                  case 2:
                     marker=&(markers_i->first_became_satellite);
                     break;
                  case 3:
                     marker=&(markers_i->joined_current_parent_top);
                     break;
                  case 4:
                     marker=&(markers_i->peak_mass);
                     break;
                  case 5:
                     marker=&(markers_i->half_peak_mass);
                     break;
                  case 6:
                     marker=&(markers_i->merger_33pc_remnant);
                     break;
                  case 7:
                     marker=&(markers_i->merger_33pc_host);
                     break;
                  case 8:
                     marker=&(markers_i->merger_33pc_merger);
                     break;
                  case 9:
                     marker=&(markers_i->merger_10pc_remnant);
                     break;
                  case 10:
                     marker=&(markers_i->merger_10pc_host);
                     break;
                  case 11:
                     marker=&(markers_i->merger_10pc_merger);
                     break;
               }
               if(halo_index>=0){
                  // Set the marker that we have just read
                  int marker_index;
                  int test[5];
                  marker_index=find_index_int(test,
                                              halo_index,
                                              (size_t)(n_halos_local[halo_snap]),
                                              file_index_local_index[halo_snap]);

                  //marker_index=find_index_int(file_index_local[halo_snap],
                  //                            halo_index,
                  //                            (size_t)(n_halos_local[halo_snap]),
                  //                            file_index_local_index[halo_snap]);
                  (*marker)=nodes_local[halo_snap][marker_index];
                  // Sanity check
                  if((*marker)->file_index!=halo_index)
                     SID_trap_error("The halo identified as a marker does not have the correct file index (ie %d!=%d)",ERROR_LOGIC,
                                    (*marker)->file_index,halo_index);
                  if((*marker)->snap_tree!=halo_snap)
                     SID_trap_error("The halo identified as a marker does not have the correct snapshot (ie %d!=%d)",ERROR_LOGIC,
                                    (*marker)->snap_tree,halo_snap);
               }
            }
         }
         if(flag_keep)
            i_found++;
      }
      
      // Close the file and clean-up 
      free_SID_fp_buffer(&fp_in_buffer);
      SID_fclose(&fp_in);

      // Sanity check
      if(i_found!=n_halos_local[i_snap])
         SID_trap_error("Failed to load all local halos (ie. %d!=%d).",ERROR_LOGIC,i_found,n_halos_local[i_snap]);

      SID_log("Done.",SID_LOG_CLOSE);
   }

   // Clean-up
   for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
      SID_free(SID_FARG nodes_local[i_snap]);
      SID_free(SID_FARG file_index_local[i_snap]);
      SID_free(SID_FARG file_index_local_index[i_snap]);
   }
   SID_free(SID_FARG nodes_local);
   SID_free(SID_FARG file_index_local);
   SID_free(SID_FARG file_index_local_index);

   SID_log("Done.",SID_LOG_CLOSE);
}

