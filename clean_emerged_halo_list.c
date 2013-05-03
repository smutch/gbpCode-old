#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees.h>

void clean_emerged_halo_list(tree_horizontal_info *halos_i,
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
         // ... and check all of their descendants ...
         tree_horizontal_info *current;
         current=halos_i[i_halo].descendant.halo;
         if(current!=NULL)
            k_file=current->file;
         l_file=k_file;
         while(current!=NULL && k_file>=l_file && k_file<MIN(n_files,i_file+(n_search+1))){
            bridge_info *bridge;
            int          j_halo;
            for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;){
               bridge=&(halos_i[i_halo].bridges[j_halo]);
               if(bridge->halo==current){
                  // Remove halos by decrementing the counter ...
                  halos_i[i_halo].n_bridges--;
                  // ... and sliding all the halos down ...
                  int k_halo;
                  for(k_halo=j_halo;k_halo<halos_i[i_halo].n_bridges;k_halo++)
                     memcpy(&(halos_i[i_halo].bridges[k_halo]),&(halos_i[i_halo].bridges[k_halo+1]),sizeof(bridge_info));
               }
               else
                 j_halo++; // We only need to increment the counter if we don't find a match
            }
            current=current->descendant.halo;
            l_file=k_file;
            if(current!=NULL)
               k_file =current->file;
         }

         // Since we may have removed items, we might not have a bridged halo any more.
         if(halos_i[i_halo].n_bridges<1){
            halos_i[i_halo].type&=(~TREE_CASE_BRIDGED);
            SID_free(SID_FARG halos_i[i_halo].bridges);
            halos_i[i_halo].n_bridges=0;
         }
         // If this halo is still bridged, label the halos in the remaining back-match list as candidate emerged halos
         else{
            bridge_info *bridge;
            int          j_halo;
            for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;j_halo++){
               bridge = &(halos_i[i_halo].bridges[j_halo]);
               bridge->halo->type|=TREE_CASE_EMERGED_CANDIDATE;
            }
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

