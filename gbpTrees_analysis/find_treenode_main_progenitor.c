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

int find_treenode_main_progenitor(tree_info *trees,tree_node_info *halo,tree_node_info **main_progenitor){
   tree_node_info *last_descendant;
   (*main_progenitor)=NULL;
   if(find_treenode_last_descendant(trees,halo,&last_descendant)){
      tree_node_info *track_halo =last_descendant->descendant;
      int             target_snap=halo->snap_tree;
      while(track_halo!=NULL){
         if(track_halo->snap_tree==target_snap) break;
         track_halo=track_halo->progenitor_first;
      }
      if(track_halo!=NULL){
         if(track_halo->snap_tree==target_snap){
            // return NULL if the halo is it's own main progenitor
            if(track_halo!=halo)
               (*main_progenitor)=track_halo;
            return(TRUE);
         }
         else
            return(FALSE);
      }
      else
         return(FALSE);
   }
   else
      return(FALSE);
}

