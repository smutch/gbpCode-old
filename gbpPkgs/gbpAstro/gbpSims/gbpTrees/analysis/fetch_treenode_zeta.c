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

double fetch_treenode_zeta(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL){
      tree_node_info *descendant=halo->descendant;
      if(descendant!=NULL){
         double M_i=fetch_treenode_Mpeak(trees,halo);
         double M_1=fetch_treenode_Mpeak(trees,descendant->progenitor_first);
         return(M_i/M_1);
      }
   }
   return(-1.);
}

