#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

// This must not be called until the main progenitors for this snapshot have been finalized.
void identify_emerged_halo_candidates(tree_horizontal_info *halos_i,
                                      int                   n_halos_i,
                                      int                   n_search){
   SID_log("Identify candidate emerged halos...",SID_LOG_OPEN|SID_LOG_TIMER);
   // Check each halo in this snapshot ...
   for(int i_halo=0;i_halo<n_halos_i;i_halo++){
      int n_back_matches=halos_i[i_halo].n_back_matches;
      if(n_back_matches>1){
         // ... check all the halos back-matched to it ...
         back_match_info *back_matches=halos_i[i_halo].back_matches;
         for(int i_back_match=0;i_back_match<n_back_matches;i_back_match++){
            tree_horizontal_info *back_match_i_halo=back_matches[i_back_match].halo;
            // ... only consider back matches which have not been placed 
            //     in the bridged halo's main progenitor line ...
            if(back_match_i_halo->id!=halos_i[i_halo].id){
               // Check that this backmatch does not have a progenitor at least as immediate as this snapshot.
               //    This can happen when that halo is not a 2way match with it's progenitor.
               //int flag_pass=TRUE;
               //if(back_match_i_halo->first_progenitor.halo!=NULL){
               //   int offset=(back_match_i_halo->file)-(back_match_i_halo->first_progenitor.halo->file);
               //   if(offset==1)
               //      flag_pass=FALSE;
               //}
               //if(flag_pass)
                  back_matches[i_back_match].halo->type|=TREE_CASE_EMERGED_CANDIDATE;
            }
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

