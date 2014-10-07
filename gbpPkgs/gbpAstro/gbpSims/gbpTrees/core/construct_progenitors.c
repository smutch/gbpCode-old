#include <gbpLib.h>
#include <gbpTrees_build.h>

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
                           int     n_halos_max,
                           int     n_files,
                           int     flag_fix_bridges,
                           int    *max_id,
                           int    *n_halos_1_matches,
                           int    *n_halos_2_matches,
                           char   *filename_root_matches,
                           char   *group_text_prefix,
                           int     flag_match_subgroups){
   SID_log("Constructing progenitors...",SID_LOG_OPEN|SID_LOG_TIMER);
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
                   j_read_1,j_read_2,n_halos_max,
                   flag_match_subgroups,
                   n_halos_1_matches,
                   n_halos_2_matches,
                   n_particles,
                   NULL,
                   n_subgroups_group[j_file_1%n_wrap],
                   n_subgroups_group[j_file_2%n_wrap],
                   match_id,
                   match_score,
                   match_index,
                   F_GOODNESS_OF_MATCH);

      // Store halo sizes
      if(!flag_fix_bridges && i_search==0){
         for(i_halo=0;i_halo<(*n_halos_1_matches);i_halo++)
            halos[i_file%n_wrap][i_halo].n_particles=n_particles[i_halo];
      }

      // Perform matching for all the halos in i_file_1.  This loop first tries to identify a
      //   first/default match for each halo in this snapshot, ultimately dealing with all simple 
      //   matches and dropped halos.  Matches to emerged halos (either coming from the default match, or to 
      //   bridged halos in it's descendant line), get dealt with subsequently.
      tree_horizontal_info *halos_i=halos[j_file_1%n_wrap];
      tree_horizontal_info *halos_j=halos[j_file_2%n_wrap];
      for(i_halo=0;i_halo<(*n_halos_1_matches);i_halo++){
         // Check for and set first matches
         tree_horizontal_info *forematch_i=halos_i[i_halo].forematch_default.halo;
         if(forematch_i==NULL){
            // If this halo has been matched to something in i_file_2 ...
            int my_descendant_index;
            my_descendant_index=match_id[i_halo];
            if(my_descendant_index>=0){
               if(my_descendant_index<(*n_halos_2_matches)){
                  halos_i[i_halo].forematch_first.halo   =&(halos_j[my_descendant_index]);
                  halos_i[i_halo].forematch_first.score  =match_score[i_halo];
                  halos_i[i_halo].forematch_default.halo =&(halos_j[my_descendant_index]);
                  halos_i[i_halo].forematch_default.score=match_score[i_halo];
                  // If there is a back match from this descendant, then we are done with this halo
                  if(check_if_descendant_is_back_matched(&(halos_i[i_halo]),&(halos_j[my_descendant_index])))
                     halos_i[i_halo].type&=(~TREE_CASE_UNPROCESSED); 
               }
               else
                  SID_log_warning("descendant ID out of bounds (ie. %d>%d) in snapshot %03d -> snapshot %03d %sgroup matching for i_halo=%d.",
                                  SID_WARNING_DEFAULT,my_descendant_index,(*n_halos_2_matches)-1,j_read_1,j_read_2,group_text_prefix,i_halo);
            }
         }
         // Scan the candidate emerged halos in the descendant line
         //    of the current default match. We don't have to address
         //    bridged halos matched above in this iteration of the loop
         //    since all emerged halos will be from the match files we will
         //    subsequently read (hence the 'else' is ok).  The check on 
         //    TREE_CASE_UNPROCESSED will fail if the default match is to
         //    a halo back matched to this one (this check is done above)
         //    or if/once we have determined that there are no more emerged 
         //    candidates to be checked for this halo.
         else if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_UNPROCESSED)){
            int n_unchecked=0;
            int flag_unchanged=TRUE;
            while(forematch_i!=NULL && flag_unchanged){
               if(check_mode_for_flag(forematch_i->type,TREE_CASE_BRIDGED)){
                  // Loop over all the candidate emerged halos identified with this
                  //   halo in the descendnat line of the current default match
                  back_match_info *back_matches  =halos[(forematch_i->file)%n_wrap][forematch_i->index].back_matches;
                  int              n_back_matches=halos[(forematch_i->file)%n_wrap][forematch_i->index].n_back_matches;
                  for(int k_halo=0;k_halo<n_back_matches;k_halo++){
                     tree_horizontal_info *current_back_match=back_matches[k_halo].halo;
                     if(current_back_match->file==j_file_2 && current_back_match->index==match_id[i_halo]){
                        if(check_validity_of_emerged_match(&(halos_i[i_halo]),&(back_matches[k_halo]),n_search)){
                           halos_i[i_halo].forematch_default.halo =&(halos[j_file_2%n_wrap][current_back_match->index]);
                           halos_i[i_halo].forematch_default.score=match_score[i_halo];
                           flag_unchanged=FALSE;
                        }
                     }
                     else if((back_matches[k_halo].halo->file)>j_file_2)
                        n_unchecked++;
                  }
               }
               forematch_i=forematch_i->descendant.halo;
               if(forematch_i!=NULL){
                  if((forematch_i->file-i_file)>=n_search || forematch_i->snap==i_read_stop)
                     forematch_i=NULL;
               }
            }
            if(n_unchecked==0 && flag_unchanged)
               halos_i[i_halo].type&=(~TREE_CASE_UNPROCESSED); 
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

