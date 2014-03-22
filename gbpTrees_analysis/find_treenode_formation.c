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

int find_treenode_formation(tree_info *trees,tree_node_info *halo,double f,tree_node_info **formation_progenitor){
   if(halo!=NULL){
      halo_properties_info *properties=fetch_treenode_properties(trees,halo);
      tree_node_info *current_halo=halo;
      (*formation_progenitor)     =halo;
      double M_target;
      double M_vir;
      double M_peak=0.;
      while(current_halo!=NULL){
         M_vir   =properties->M_vir;
         M_peak  =MAX(M_peak,M_vir);
         M_target=f*M_peak;
         if(M_vir>M_target)
            (*formation_progenitor)=current_halo;
         current_halo=current_halo->progenitor_first;
      }
      return(TRUE);
   }
   (*formation_progenitor)=NULL;
   return(FALSE);
}

