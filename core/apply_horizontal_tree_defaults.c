#include <gbpLib.h>
#include <gbpTrees_build.h>

void apply_horizontal_tree_defaults(int                    n_halos_1_matches,
                                    int                    n_halos_i,
                                    tree_horizontal_info **halos,
                                    tree_horizontal_info  *halos_i,
                                    int                    i_file,
                                    int                    n_wrap,
                                    int                   *max_id,
                                    int                   *max_tree_id){
   SID_log("Applying defaults to unprocessed halos...",SID_LOG_OPEN|SID_LOG_TIMER);
   // ... first deal with unprocessed matches to bridged halos (apply default behavior)
   int i_halo;
   for(i_halo=0;i_halo<n_halos_1_matches;i_halo++){
      if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED)){
         // The descendant info is already set.  Just increment it's progenitor counter and set this halo's info.
         //   Leave the target flag untouched so we can later identify which BRIDGE_PROGENITORS were found
         halos_i[i_halo].type|=(TREE_CASE_BRIDGE_FINALIZE|TREE_CASE_BRIDGE_DEFAULT);
         set_halo_and_descendant(halos,
                                 i_file,
                                 i_halo,
                                 halos_i[i_halo].bridge_forematch.halo->file,
                                 halos_i[i_halo].bridge_forematch.halo->index,
                                 halos_i[i_halo].bridge_forematch.score,
                                 max_id,
                                 n_wrap);
      }
   }
   // ... then assign flags for halos not successfully processed.  Call them strays.
   //     These will include halos which have stopped existing without merging with
   //     anything over the search range.
   for(i_halo=0;i_halo<n_halos_i;i_halo++){
      if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_UNPROCESSED)){
         halos_i[i_halo].type   |=TREE_CASE_STRAYED;
         halos_i[i_halo].type   &=(~TREE_CASE_UNPROCESSED);
         halos_i[i_halo].id      =(*max_id)++;
         halos_i[i_halo].tree_id =(*max_tree_id)++;
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

