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

int find_treenode_markers(tree_info *trees,tree_node_info *halo,tree_markers_info *markers){
   find_treenode_leaf(           trees,halo,&(markers->leaf));
   find_treenode_last_descendant(trees,halo,&(markers->descendant_last));
   find_treenode_main_progenitor(trees,halo,&(markers->progenitor_main));
   find_treenode_accretion(      trees,halo,&(markers->progenitor_accretion_first),
                                            &(markers->progenitor_accretion_last));
   find_treenode_formation(      trees,halo,0.5,
                                            &(markers->progenitor_formation));
   find_treenode_main_progenitor(trees,halo,&(markers->progenitor_main));
   markers->flag_halo_is_main_progenitor=check_treenode_if_main_progenitor(halo);
}

