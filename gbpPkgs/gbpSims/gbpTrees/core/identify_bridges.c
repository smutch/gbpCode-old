#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

// This must not be called until the main progenitors for this snapshot have been finalized.
void identify_bridges(tree_horizontal_info *halos_i,
                      int                   n_halos_i,
                      int                   n_search){
   SID_log("Identify bridges and candidate emerged halos...",SID_LOG_OPEN|SID_LOG_TIMER);
   // Check each halo in this snapshot ...
   for(int i_halo=0;i_halo<n_halos_i;i_halo++){
      int n_back_matches=halos_i[i_halo].n_back_matches;
      if(n_back_matches>0){
         // ... check all the halos back-matched to it ...
         match_info *back_matches       =halos_i[i_halo].back_matches;
         int         n_emerged_candidate=0;
         for(int i_back_match=0;i_back_match<n_back_matches;i_back_match++){
            tree_horizontal_info *back_match_i_halo=back_matches[i_back_match].halo;
            // ... only consider back matches which have not been placed 
            //     in the bridged halo's main progenitor line, or which
            //     have not been given progenitors (possible in asymetric
            //     matches, where a halo is a back match to one halo,
            //     but has a fore match from another, for example).
            //     Both of these cases can be checked for by just seeing
            //     if the halo does not have a progenitor.
            if(back_match_i_halo->first_progenitor.halo==NULL){
               back_match_i_halo->type|=TREE_CASE_EMERGED_CANDIDATE;
               n_emerged_candidate++;
            }
         }
         // Set bridge halo flag
         if(n_emerged_candidate>0)
            halos_i[i_halo].type|=TREE_CASE_BRIDGED;
         else
            halos_i[i_halo].type&=(~TREE_CASE_BRIDGED);
      }
      else
         halos_i[i_halo].type&=(~TREE_CASE_BRIDGED);
   }
   SID_log("Done.",SID_LOG_CLOSE);
}
