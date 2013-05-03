#include <gbpLib.h>
#include <gbpTrees.h>

void construct_progenitors(tree_horizontal_info **halos,
                           tree_horizontal_info  *halos_i,
                           int   **n_subgroups_group,
                           int     n_halos_i,
                           int    *match_id,
                           float  *match_score,
                           size_t *match_index,
                           int    *n_particles,
                           int     i_file,
                           int     i_read,
                           int     i_read_start,
                           int     i_read_stop,
                           int     i_read_step,
                           int     n_search,
                           int     n_wrap,
                           int     n_files,
                           int     flag_fix_bridges,
                           int    *max_id,
                           int    *n_halos_1_matches,
                           int    *n_halos_2_matches,
                           char   *filename_root_matches,
                           char   *group_text_prefix,
                           int     flag_match_subgroups){
   SID_log("Constructing progenitors from forward-matching...",SID_LOG_OPEN|SID_LOG_TIMER);
   int j_file_1;
   int j_file_2;
   int j_read_1;
   int j_read_2;
   int i_search;
   int i_halo;
   for(j_file_1  =i_file,
         j_file_2=i_file+1,
         j_read_1=i_read,
         j_read_2=i_read+i_read_step,
         i_search=0;
       j_read_2<=i_read_stop && i_search<n_search;
       j_file_2++,
         j_read_2+=i_read_step,
         i_search++){

      // Read forward-matching
      read_matches(filename_root_matches,
                   j_read_1,j_read_2,
                   flag_match_subgroups,
                   n_halos_1_matches,
                   n_halos_2_matches,
                   n_particles,
                   NULL,
                   n_subgroups_group[j_file_1%n_wrap],
                   n_subgroups_group[j_file_2%n_wrap],
                   match_id,
                   match_score,
                   match_index);

      // Store halo sizes
      if(!flag_fix_bridges){
         if(i_search==0){
            for(i_halo=0;i_halo<(*n_halos_1_matches);i_halo++)
               halos[i_file%n_wrap][i_halo].n_particles=n_particles[i_halo];
         }
      }

      // Perform matching for all the halos in i_file_1.  This loop should deal completely with
      //   all simple matches and dropped halos.  It also identifies matches to bridges, which
      //   require special treatment in the loop that follows (to look for matches to emergent halos)
      //   and at the end of the loop over j_read/i_search to finalize those not matched to emerged halos.
      for(i_halo=0;i_halo<(*n_halos_1_matches);i_halo++){
         // If this halo hasn't been processed during earlier searches ...
         if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_UNPROCESSED)){
            int my_descendant_index;
            my_descendant_index=match_id[i_halo];
            // If this halo has been matched to something in i_file_2 ...
            if(my_descendant_index>=0){
               if(my_descendant_index<(*n_halos_2_matches))
                  set_halo_and_descendant(halos,
                                          i_file,
                                          i_halo,
                                          j_file_2,
                                          my_descendant_index,
                                          match_score[i_halo],
                                          max_id,
                                          n_wrap);
               else
                  SID_log_warning("descendant ID out of bounds (ie. %d>%d) in snapshot %03d -> snapshot %03d %sgroup matching for i_halo=%d.",
                                  SID_WARNING_DEFAULT,my_descendant_index,(*n_halos_2_matches)-1,j_read_1,j_read_2,group_text_prefix,i_halo);
            }
         }
      }

      // Try to match halos-matched-to-bridges to the candidate emergent halos identified in the bridge-lists
      for(i_halo=0;i_halo<(*n_halos_1_matches);i_halo++){
         if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED)){
            // Loop over all the emergent halos identified with the bridge that i_halo has been matched to
            int          n_list;
            bridge_info *bridges;
            if(halos_i[i_halo].bridge_forematch.halo==NULL)
              SID_trap_error("Bridge match not defined during emerged halo search.",ERROR_LOGIC);
            bridges=halos[(halos_i[i_halo].bridge_forematch.halo->file)%n_wrap][halos_i[i_halo].bridge_forematch.halo->index].bridges;
            if(bridges==NULL)
              SID_trap_error("Bridges not defined during emerged halo search.",ERROR_LOGIC);
            n_list=halos[(halos_i[i_halo].bridge_forematch.halo->file)%n_wrap][halos_i[i_halo].bridge_forematch.halo->index].n_bridges;
            int k_halo;
            for(k_halo=0;k_halo<n_list;k_halo++){
               if(bridges[k_halo].halo->file==j_file_2 && match_id[i_halo]==bridges[k_halo].halo->index){
                  set_halo_and_descendant(halos,
                                          i_file,
                                          i_halo,
                                          bridges[k_halo].halo->file,
                                          bridges[k_halo].halo->index,
                                          match_score[i_halo],
                                          max_id,
                                          n_wrap);
               }
            }
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);

}

