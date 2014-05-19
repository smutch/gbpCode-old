#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void finalize_bridged_halo_list(tree_horizontal_info *halos_i,
                                int                   n_halos_i,
                                int                   i_file,
                                int                   n_search,
                                int                   n_files){
   SID_log("Removing main progenitors from candidate emerged halo lists...",SID_LOG_OPEN|SID_LOG_TIMER);
   int i_halo;
   int k_file;
   int l_file;
   for(i_halo=0;i_halo<n_halos_i;i_halo++){
      // Check all bridged halos ...
      if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGED)){
         // Since we may have removed item
         if(halos_i[i_halo].n_back_matches<=1){
            halos_i[i_halo].type&=(~TREE_CASE_BRIDGED);
            SID_free(SID_FARG halos_i[i_halo].back_matches);
            halos_i[i_halo].n_back_matches=0;
         }
         // Label the halos in the back-match list as candidate emerged halos
         else{
            back_match_info *back_match;
            int          j_halo;
            for(j_halo=0;j_halo<halos_i[i_halo].n_back_matches;j_halo++){
               back_match = &(halos_i[i_halo].back_matches[j_halo]);
               back_match->halo->type|=TREE_CASE_EMERGED_CANDIDATE;
            }
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

