#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

int find_treenode_Mpeak(tree_info       *trees,
                        tree_node_info  *halo,
                        tree_node_info **halo_peak){
   SID_trap_error("find_treenode_Mpeak() not working.",ERROR_LOGIC);
   /*
   (*halo_peak)=halo;
   if(halo!=NULL){
      tree_node_info *current_halo=halo;
      // Scan to the halo's leaf
      while((current_halo->progenitor_first)!=NULL){
         current_halo=current_halo->progenitor_first;
      }
      // Follow the dominant halo line back.  Wait, this wont bring us 
      //    back to halo in general ... need to think about this ...
      return(TRUE);
   }
   */
   return(FALSE);
}

