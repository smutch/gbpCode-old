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

tree_node_info *fetch_treenode_reference(tree_info *trees,tree_node_info *halo){
   tree_node_info *reference_halo=halo;
   if(trees->trees_reference!=NULL){
      int snapshot_halo      =trees->snap_list[halo->snap_tree];
      int file_index_halo    =fetch_treenode_file_index(trees,halo);
      int snap_tree_reference=find_treesnap_snap(trees->trees_reference,snapshot_halo);
      if(!find_tree_node(trees->trees_reference,snap_tree_reference,file_index_halo,halo->parent==NULL,&reference_halo))
         SID_trap_error("Could not find halo in reference trees (snapshot=%d file_index=%d).",ERROR_LOGIC,snapshot_halo,file_index_halo);
   }
   return(reference_halo);
}

