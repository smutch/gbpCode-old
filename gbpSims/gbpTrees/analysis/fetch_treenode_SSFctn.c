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

double fetch_treenode_SSFctn(tree_info *trees,tree_node_info *halo){
   double SSFctn=0.;
   if(halo!=NULL){
      // Sum the halo's substructure mass
      double M_vir_substructure=0.;
      tree_node_info *current_substructure=halo->substructure_first;
      while(current_substructure!=NULL){
         M_vir_substructure+=fetch_treenode_Mvir(trees,current_substructure);
         current_substructure=current_substructure->substructure_next;
      }
      return(M_vir_substructure/fetch_treenode_Mvir(trees,halo));
   }
   return(SSFctn);
}

