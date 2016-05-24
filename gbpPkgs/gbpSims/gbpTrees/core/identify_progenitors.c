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

   // Now that all first and default forematch pointers have been set,
   //    determine the best match to any halos matched to from this snapshot
   //    that don't have a progenitor yet.
   // We need to wait for both the first and default pointers to be properly 
   //    set to do this, so we need to do this now as a seperate and subsequent 
   //    loop.  When matches are finalized, these best pointers will be the first 
   //    choice for constructing the trees.  The default pointers get used otherwise.
   // All this is done to deal with the following situation:
   //    It may be that the core of the halo we choose to be the progenitor here 
   //    is ejected from the merger remnant.  In this case, the optimal choice
   //    would be to label the remnant as a large fragmented halo, and catch
   //    the ejected core as an emerged halo.  Unfortunately, *any* other merger
   //    to the remnant - no matter how small - will get called a main progenitor,
   //    turn the otherwise fragmented remnant into a normal halo (probably with
   //    a BIG mass jump) and ultimately lead to one merger being counted as two.
   //    The choice we make below is a necessary compromise to maximize the chance
   //    of getting the count (and merger ratio) of major mergers correct.
   if(flag_fix_bridges){
      // Loop once to choose the best progenitors ...
      for(int i_halo=0;i_halo<(*n_halos_1_matches);i_halo++){
         tree_horizontal_info *halo_i=&(halos_i[i_halo]);
         tree_horizontal_info *halo_j=halo_i->forematch_first.halo;
         // ... if this halo has been matched to something ...
         if(halo_j!=NULL){
            match_info forematch_new;
            forematch_new.halo        =halo_i;
            forematch_new.score       =halo_i->forematch_first.score;
            forematch_new.flag_two_way=halo_i->forematch_first.flag_two_way;
            // If this is the first match to the halo it's matched to, 
            //    choose it as the best progenitor by default ...
            if(halo_j->forematch_best.halo==NULL)
               memcpy(&(halo_j->forematch_best),&forematch_new,sizeof(match_info));
            // ... else, if this is a subsequent match, decide if it is better.
            //     Make sure that the most massive halo is a progenitor.  An
            //     attempt was made to use a halo not matched to an emerged
            //     halo for the progenitor here.  Unfortunately, this does not
            //     work in cases where the cores of both halos are ejected and
            //     another small halo is present to become a progenitor.  This
            //     leads to one major merger being turned into two.  Simply
            //     choosing the largest halo here conserves merger count at least.
            //     If it's core becomes fragmented, it will be labled ejected below.
            else if(check_if_match_is_bigger(halo_j,&(halo_j->forematch_best),&forematch_new))
               memcpy(&(halo_j->forematch_best),&forematch_new,sizeof(match_info));
         }
      }
      // Loop a second time to set TREE_CASE_FRAGMENTED_EJECTED flags ...
      for(int i_halo=0;i_halo<(*n_halos_1_matches);i_halo++){
         tree_horizontal_info *halo_i=&(halos_i[i_halo]);
         tree_horizontal_info *halo_j=halo_i->forematch_first.halo;
         // ... if this halo has been matched to something ...
         if(halo_j!=NULL){
            tree_horizontal_info *halo_k=halo_j->forematch_best.halo;
            int flag_is_best_progenitor=((halo_j->forematch_best.halo)==(halo_i));
            // ... and it has been selected as a best match above ...
            if(flag_is_best_progenitor){
               // ... but was identified as a possible progenitor of an en emerged halo,
               //     then mark that ejected halo as an EJECTED fragment.  This may get
               //     overridden later if a progenitor is found in subsequent snapshots.
               int flag_matched_to_emerged=(!((halo_i->forematch_first.halo)==(halo_i->forematch_default.halo)));
               if(flag_matched_to_emerged)
                  halo_i->forematch_default.halo->type|=TREE_CASE_FRAGMENTED_EJECTED;
            }
         }
      }
   }

   SID_log("Done.",SID_LOG_CLOSE);
}
