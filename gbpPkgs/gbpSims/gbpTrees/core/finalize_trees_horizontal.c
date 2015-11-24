#include <gbpLib.h>
#include <gbpTrees_build.h>

void finalize_trees_horizontal(int                    n_halos_1_matches,
                               int                    n_halos_i,
                               tree_horizontal_info **halos,
                               tree_horizontal_info  *halos_i,
                               int                    i_file,
                               int                    n_search,
                               int                    n_wrap,
                               int                   *max_id,
                               int                   *max_tree_id){
   SID_log("Finalizing...",SID_LOG_OPEN|SID_LOG_TIMER);

   // ... finalize the sucessful matches ...
   for(int i_halo=0;i_halo<n_halos_1_matches;i_halo++){
      if((halos_i[i_halo].forematch_default.halo)!=NULL){
         // First, if this halo's initial match is the best from this
         //   snapshot to a halo with no progenitors yet, use it to
         //   make sure that every halo with a match gets a progenitor 
         //   (rather then let all matches get converted to emerged 
         //   halos later).  If a halo has a default match, it must
         //   have a first match (by construction) so no need to check
         //   for it.
         if((halos_i[i_halo].forematch_first.halo->forematch_best.halo)==(&(halos_i[i_halo]))){
            add_progenitor_to_halo(halos,
                                   i_file,
                                   i_halo,
                                   halos_i[i_halo].forematch_first.halo->file,
                                   halos_i[i_halo].forematch_first.halo->index,
                                   halos_i[i_halo].forematch_first.score,
                                   max_id,
                                   n_wrap);
            // Add 2WAY match flags
            if(halos_i[i_halo].forematch_first.flag_two_way)
               halos_i[i_halo].type|=TREE_CASE_2WAY_MATCH;
         }
         // ... else, use the default pointer.
         else{
            add_progenitor_to_halo(halos,
                                   i_file,
                                   i_halo,
                                   halos_i[i_halo].forematch_default.halo->file,
                                   halos_i[i_halo].forematch_default.halo->index,
                                   halos_i[i_halo].forematch_default.score,
                                   max_id,
                                   n_wrap);
            // Add 2WAY match flags
            if(halos_i[i_halo].forematch_default.flag_two_way)
               halos_i[i_halo].type|=TREE_CASE_2WAY_MATCH;

            // Turn off TREE_CASE_SET_BY_BACKMATCH if the default match used != first match 
            if(halos_i[i_halo].forematch_default.halo!=halos_i[i_halo].forematch_first.halo)
               halos_i[i_halo].type&=(~TREE_CASE_SET_BY_BACKMATCH);
         }
      }
   }

   // ... then assign flags for halos not successfully processed.  Call them strays.
   //     These will include halos which have stopped existing without merging with
   //     anything over the search range.  Do some other final cleaning of the flags as well.
   for(int i_halo=0;i_halo<n_halos_i;i_halo++){
      if(halos_i[i_halo].descendant.halo==NULL){
         halos_i[i_halo].type   |=TREE_CASE_STRAYED;
         halos_i[i_halo].type   &=(~TREE_CASE_UNPROCESSED);
         halos_i[i_halo].id      =(*max_id)++;
         halos_i[i_halo].tree_id =(*max_tree_id)++;
      }
      // Turn off TREE_CASE_EMERGED_CANDIDATE if TREE_CASE_EMERGED or TREE_CASE_FRAGMENTED_NEW are on.
      if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_EMERGED) || check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_FRAGMENTED_NEW))
         halos_i[i_halo].type&=(~TREE_CASE_EMERGED_CANDIDATE);
   }

   SID_log("Done.",SID_LOG_CLOSE);
}

