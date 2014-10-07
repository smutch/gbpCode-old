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

   // ... first, finalize the sucessful matches ...
   int i_halo;
   for(i_halo=0;i_halo<n_halos_1_matches;i_halo++){
      if((halos_i[i_halo].forematch_default.halo)!=NULL)
         add_progenitor_to_halo(halos,
                                i_file,
                                i_halo,
                                halos_i[i_halo].forematch_default.halo->file,
                                halos_i[i_halo].forematch_default.halo->index,
                                halos_i[i_halo].forematch_default.score,
                                max_id,
                                n_wrap);
   }

   // ... then assign flags for halos not successfully processed.  Call them strays.
   //     These will include halos which have stopped existing without merging with
   //     anything over the search range.
   for(i_halo=0;i_halo<n_halos_i;i_halo++){
      if(halos_i[i_halo].descendant.halo==NULL){
         halos_i[i_halo].type   |=TREE_CASE_STRAYED;
         halos_i[i_halo].type   &=(~TREE_CASE_UNPROCESSED);
         halos_i[i_halo].id      =(*max_id)++;
         halos_i[i_halo].tree_id =(*max_tree_id)++;
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

