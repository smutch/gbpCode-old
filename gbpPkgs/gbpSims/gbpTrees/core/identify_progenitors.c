#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpTrees_build.h>

void identify_progenitors(tree_horizontal_info **halos,
                          tree_horizontal_info  *halos_i,
                          int   **n_subgroups_group,
                          int     n_halos_i,
                          int    *match_id,
                          float  *match_score,
                          size_t *match_index,
                          char   *match_flag_two_way,
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
   SID_log("Identifying progenitors...",SID_LOG_OPEN|SID_LOG_TIMER);

   // Loop over the scan range
   int j_file_1=i_file;
   int j_file_2=i_file+1;
   int j_read_1=i_read;
   int j_read_2=i_read+i_read_step;
   int i_search=0;
   for(;
       j_read_2<=i_read_stop && i_search<n_search;
       j_file_2++,j_read_2+=i_read_step,i_search++){

      // Read forward-matching for searching for dropped and emerged matches
      //    (or for setting descendants if bridge fixing is turned off)
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
                   match_flag_two_way,
                   F_GOODNESS_OF_MATCH);

      // Determine the back match which would work best as a descendant
      if(i_search==0){
         if(flag_fix_bridges){
            for(int i_halo=0;i_halo<(*n_halos_1_matches);i_halo++){
               tree_horizontal_info *halo_i=&(halos_i[i_halo]);
               int n_back_matches=halo_i->n_back_matches;
               if(n_back_matches>0){
                  match_info *back_matches=halo_i->back_matches;
                  // Loop over back matches
                  for(int k_halo=0;k_halo<n_back_matches;k_halo++){
                     // Choose the most immediate back match with the best score
                     match_info *current_back_match=&(back_matches[k_halo]);
                     if(k_halo==0)
                        memcpy(&(halo_i->forematch_first),current_back_match,sizeof(match_info));
                     else if(check_if_match_is_better(halo_i,&(halo_i->forematch_first),current_back_match))
                        memcpy(&(halo_i->forematch_first),current_back_match,sizeof(match_info));
                  }
                  // Initialize the default pointer to the best back match pointer
                  memcpy(&(halo_i->forematch_default),&(halo_i->forematch_first),sizeof(match_info));
               }
            }
         }
         // Store halo sizes.  If flag_fix_bridges is on, then this has already been done.
         else{
            for(int i_halo=0;i_halo<(*n_halos_1_matches);i_halo++)
               halos[i_file%n_wrap][i_halo].n_particles=n_particles[i_halo];
         }
      }

      // Examine the matching for all the halos in i_file_1.  This loop scans the progenitor line of the
      //   current 'default' match of each halo (initialized to be the first match identified above), checking 
      //   the back matches to the descendant line for back matches to halos
      //   marked as emerged halo candidates.  Such matches (either coming from emerged candidates of the
      //   first match, or to emerged candidates of bridged halos in it's descendant line) get dealt with 
      //   subsequently and will be called the 'default' match.  Decisions about whether to use the 'first'
      //   match or the 'default' match get made later.
      if(flag_fix_bridges){
         tree_horizontal_info *halos_j=halos[j_file_2%n_wrap];
         for(int i_halo=0;i_halo<(*n_halos_1_matches);i_halo++){
            // The check on TREE_CASE_UNPROCESSED will fail if/once we have scanned
            //    all the possible emerged candidates (and emerged candidates of good 
            //    emerged matches, etc) for this halo.
            tree_horizontal_info *halo_i     =&(halos_i[i_halo]);
            tree_horizontal_info *forematch_i=halo_i->forematch_default.halo;
            // Only check unprocessed halos that have been given a descendant
            if(forematch_i!=NULL && check_mode_for_flag(halo_i->type,TREE_CASE_UNPROCESSED)){
               // Only scan for emerged halos if we are working with
               //    matches past the current default match
               if(j_file_2>forematch_i->file){
                  // Walk the descendant line of the default match until a change happens (or we reach the extent of the search interval)
                  int flag_unchecked=FALSE;
                  while(forematch_i!=NULL){
                     int n_back_matches=forematch_i->n_back_matches;
                     if(n_back_matches>0){
                        match_info *back_matches=forematch_i->back_matches;
                        // Loop over all the candidate emerged halos identified with this
                        //   halo in the descendant line of the current default match
                        for(int k_halo=0;k_halo<n_back_matches && !flag_unchecked;k_halo++){
                           tree_horizontal_info *current_back_match=back_matches[k_halo].halo;
                           // If this back match is in the snapshot we are currently checking
                           if(current_back_match->file==j_file_2){
                              if(current_back_match->index==match_id[i_halo]){
                                 // If this is a good match to an emerged halo, then set a new default.  Do not remove the UNPROCESSED flag
                                 //    though, because we want to continue checking emerged candidates in the descendant line of the new match
                                 if(check_validity_of_emerged_match(halo_i,&(back_matches[k_halo]),n_search)){
                                    halo_i->forematch_default.halo           =current_back_match;
                                    halo_i->forematch_default.score          =match_score[i_halo];
                                    halo_i->forematch_default.flag_two_way   =match_flag_two_way[i_halo];
                                    halo_i->forematch_default.flag_back_match=FALSE;
                                    flag_unchecked=TRUE; // we want to make sure we check the back matches of the new match
                                 }
                              }
                           }
                           // Determine if there are back matches still to be checked in later snapshots
                           //    This check works because the back matches have been sorted by their file number
                           else if(current_back_match->file>j_file_2)
                              flag_unchecked=TRUE;
                        }
                     }
                     forematch_i=forematch_i->descendant.halo;
                     // Check if we've reached the extent of the search interval
                     if(forematch_i!=NULL){
                        if((forematch_i->file-i_file)>=n_search || forematch_i->snap==i_read_stop)
                           forematch_i=NULL;
                     }
                  }
                  // Stop processing this halo if there are no more emerged candaidates to check
                  if(!flag_unchecked)
                     halo_i->type&=(~TREE_CASE_UNPROCESSED); 
               }
            }
         }
      }
      // ... else, set progenitors if bridge fixing is switched off ...
      else{
         tree_horizontal_info *halos_j=halos[j_file_2%n_wrap];
         for(int i_halo=0;i_halo<(*n_halos_1_matches);i_halo++){
             int j_halo=match_id[i_halo];
             if(j_halo>=0){
                tree_horizontal_info *halo_i=&(halos[j_file_1%n_wrap][i_halo]);
                if(halo_i->forematch_first.halo==NULL){
                   tree_horizontal_info *halo_j           =&(halos[j_file_2%n_wrap][j_halo]);
                   halo_i->forematch_first.halo           =halo_j;
                   halo_i->forematch_first.score          =match_score[i_halo];
                   halo_i->forematch_first.flag_two_way   =match_flag_two_way[i_halo];
                   halo_i->forematch_first.flag_back_match=FALSE;
                }
             }
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}
