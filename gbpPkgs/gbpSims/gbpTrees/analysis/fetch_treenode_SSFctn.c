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
      int np_parent      =fetch_treenode_n_particles(trees,halo);
      int np_substructure=0;
      int np_most_massive=0;
      int np_i;
      int i_sub=0;
      tree_node_info *current_substructure=halo->substructure_first;
SID_trap_error("fetch_treenode_SSFctn() needs to be checked!",ERROR_NONE);
      while(current_substructure!=NULL){
         np_i            =fetch_treenode_n_particles(trees,current_substructure);
         np_substructure+=np_i;
         if(np_i>np_most_massive)
            np_most_massive=np_i;
         i_sub++;
         current_substructure=current_substructure->substructure_next;
      }
      np_substructure-=np_most_massive;
      SSFctn          =(double)np_substructure/(double)np_parent;
   }
   return(SSFctn);
}

