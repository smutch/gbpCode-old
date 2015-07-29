#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void identify_emerged_halo_candidates(tree_horizontal_info *halos_i,
                                      int                   n_halos_i,
                                      int                   n_search){
   SID_log("Identify candidate emerged halos...",SID_LOG_OPEN|SID_LOG_TIMER);
   // Check all the halos in this snapshot ...
   for(int i_halo=0;i_halo<n_halos_i;i_halo++){
      if((halos_i[i_halo].n_back_matches)>1){
         // ... check all the halos back-matched to it ...
         back_match_info *back_matches=halos_i[i_halo].back_matches;
         for(int i_back_match=0;i_back_match<(halos_i[i_halo].n_back_matches);i_back_match++){
            // ... all but the one identified as it's main progenitor are candidate emerged halos.
            if(back_matches[i_back_match].halo->id!=halos_i[i_halo].id)
               back_matches[i_back_match].halo->type|=TREE_CASE_EMERGED_CANDIDATE;
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

