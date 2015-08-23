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

int find_treenode_formation(tree_info       *trees,
                            tree_node_info  *halo,
                            double           fraction,
                            tree_node_info **halo_formation){
   (*halo_formation)=halo;
   if(halo!=NULL){
      double M_peak_halo=fetch_treenode_Mpeak(trees,halo);
      double M_target   =fraction*M_peak_halo;
      tree_node_info *current_halo=halo;
      while(current_halo!=NULL){
         double M_peak_current=fetch_treenode_Mpeak(trees,current_halo);
         if(M_peak_halo>M_target)
            (*halo_formation)=current_halo;
         else break;
         current_halo=current_halo->progenitor_first;
      }
      return(TRUE);
   }
   return(FALSE);
}

