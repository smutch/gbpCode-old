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
   (*halo_peak)=halo;
   if(halo!=NULL){
      tree_node_info *current_halo=halo;
      // Scan to the halo's leaf
      while((current_halo->progenitor_dominant)!=NULL){
         current_halo=current_halo->progenitor_first;
      }
      // Follow the main progenitor line back.
      return(TRUE);
   }
   return(FALSE);
}

