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

void write_treenode_markers(tree_info *trees,const char *filename_output_root,int mode){

   // Set some stuff specific to whether we're processing groups or subgroups
   char filename_output_group_text[16];
   int  flag_process_groups=FALSE;
   tree_markers_info **markers;
   int                *n_halos_array;
   if(check_mode_for_flag(mode,PRECOMPUTE_TREENODE_MARKER_GROUPS)){
      sprintf(filename_output_group_text,"groups");
      flag_process_groups=TRUE;
      markers      =trees->group_markers;
      n_halos_array=trees->n_groups_snap_local;
   }
   else if(check_mode_for_flag(mode,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS)){
      sprintf(filename_output_group_text,"subgroups");
      markers      =trees->subgroup_markers;
      n_halos_array=trees->n_subgroups_snap_local;
   }
   else
      SID_trap_error("group/subgroup mode has not been properly specified in write_treenode_markers().",ERROR_LOGIC);

   // Perform the write
   SID_log("Writing markers...",SID_LOG_OPEN|SID_LOG_TIMER);
   char filename_out_dir[MAX_FILENAME_LENGTH];
   sprintf(filename_out_dir,"%s_markers",filename_output_root);
   mkdir(filename_out_dir,02755);
   for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
      int n_halos_local;
      int n_halos_total;
      if(flag_process_groups){
         n_halos_total=trees->n_groups_catalog[i_snap];
         n_halos_local=trees->n_groups_snap_local[i_snap];
      }
      else{
         n_halos_total=trees->n_subgroups_catalog[i_snap];
         n_halos_local=trees->n_subgroups_snap_local[i_snap];
      }
      SID_log("Processing snapshot #%03d...",SID_LOG_OPEN|SID_LOG_TIMER,trees->snap_list[i_snap]);
      SID_log("(%d halos)...",SID_LOG_CONTINUE,n_halos_total);
      // Open file for writing and write header
      char filename_out[MAX_FILENAME_LENGTH];
      sprintf(filename_out,"%s/%s_%03d.dat",filename_out_dir,filename_output_group_text,trees->snap_list[i_snap]);
      FILE *fp_out;
      if(SID.I_am_Master){
         fp_out=fopen(filename_out,"w");
         fwrite(&i_snap,       sizeof(int),1,fp_out);
         fwrite(&n_halos_total,sizeof(int),1,fp_out);
      }
      // Write in catalog order
      int i_halo;
      int j_halo;
      int i_buffer;
      int buffered_count_local;
      int invalid_entry=-100;
      int n_marker_write   =12;
      int n_buffer_per_halo=1+2*n_marker_write;
      int n_buffer_halos;
      int n_buffer_halos_max=1024*1024;
      int n_buffer;
      int *buffer=(int *)SID_malloc(sizeof(int)*n_buffer_halos_max*n_buffer_per_halo);
      for(i_halo=0,buffered_count_local=0;i_halo<n_halos_total;i_halo+=n_buffer_halos){
         // Decide this buffer iteration's size
         n_buffer_halos=MIN(n_buffer_halos_max,n_halos_total-i_halo);
         n_buffer      =n_buffer_halos*n_buffer_per_halo;
         // Set the buffer to a default value smaller than the smallest possible data size
         for(i_buffer=0;i_buffer<n_buffer;i_buffer++)
            buffer[i_buffer]=invalid_entry; 
         // Determine if any of the local data is being used for this buffer
         tree_node_info *current_halo;
         if(flag_process_groups)
            current_halo=trees->first_neighbour_groups[i_snap];
         else
            current_halo=trees->first_neighbour_subgroups[i_snap];
         while(current_halo!=NULL){
            j_halo=current_halo->neighbour_index;
            // ... if so, set the appropriate buffer value
            int index_test=current_halo->file_index-i_halo;
            if(index_test>=0 && index_test<n_buffer_halos){
               int i_block=0;
               tree_markers_info *markers_i=&(markers[i_snap][j_halo]);
               buffer[index_test*n_buffer_per_halo+(i_block++)]=markers_i->flag_halo_is_main_progenitor;
               for(int i_marker=0;i_marker<n_marker_write;i_marker++){
                  tree_node_info *marker;
                  switch(i_marker){
                     case 0:
                        marker=markers_i->branch_leaf;
                        break;
                     case 1:
                        marker=markers_i->branch_root;
                        break;
                     case 2:
                        marker=markers_i->first_became_satellite;
                        break;
                     case 3:
                        marker=markers_i->joined_current_parent;
                        break;
                     case 4:
                        marker=markers_i->peak_mass;
                        break;
                     case 5:
                        marker=markers_i->half_peak_mass;
                        break;
                     case 6:
                        marker=markers_i->merger_33pc_remnant;
                        break;
                     case 7:
                        marker=markers_i->merger_33pc_host;
                        break;
                     case 8:
                        marker=markers_i->merger_33pc_merger;
                        break;
                     case 9:
                        marker=markers_i->merger_10pc_remnant;
                        break;
                     case 10:
                        marker=markers_i->merger_10pc_host;
                        break;
                     case 11:
                        marker=markers_i->merger_10pc_merger;
                        break;
                     default:
                        SID_trap_error("Marker count (%d) has exceeded the specified write count (%d).",ERROR_LOGIC,i_marker,n_marker_write);
                  }
                  if(marker!=NULL){
                     buffer[index_test*n_buffer_per_halo+(i_block++)]=marker->snap_tree;
                     buffer[index_test*n_buffer_per_halo+(i_block++)]=marker->file_index;
                  }
                  else{
                     buffer[index_test*n_buffer_per_halo+(i_block++)]=-1;
                     buffer[index_test*n_buffer_per_halo+(i_block++)]=-1;
                  }
               } // loop over markers
               if(i_block!=n_buffer_per_halo) 
                  SID_trap_error("block count (%d) does not agree with the allocation size (%d).",ERROR_LOGIC,i_block,n_buffer_per_halo);
               buffered_count_local++;
            } // if halo is in the current buffer
            current_halo=current_halo->next_neighbour;
         }
         // Doing a global max on the buffer yields the needed buffer on all ranks
         SID_Allreduce(SID_IN_PLACE,buffer,n_buffer,SID_INT,SID_MAX,SID.COMM_WORLD);
         // Perform the write
         if(SID.I_am_Master){
            // Sanity check
            for(i_buffer=0;i_buffer<n_buffer;i_buffer++){
               if(buffer[i_buffer]==invalid_entry){
                  buffer[i_buffer]=-1; //because we skip some groups, this can actually happen
                  //SID_trap_error("Illegal buffer value (%d).",ERROR_LOGIC,buffer[i_buffer]);
               }
            }
            // Write the buffer
            fwrite(buffer,sizeof(int),(size_t)n_buffer,fp_out);
         }
      } // Loop over buffer-sized chuncks
      if(SID.I_am_Master)
         fclose(fp_out);
      SID_free(SID_FARG buffer);
      SID_log("Done.",SID_LOG_CLOSE);
   }
   SID_log("Done.",SID_LOG_CLOSE);

}

