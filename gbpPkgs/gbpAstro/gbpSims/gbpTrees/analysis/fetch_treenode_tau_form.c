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

double fetch_treenode_tau_form(tree_info *trees,tree_node_info *halo){
   double tau=-1.; 
   if(halo!=NULL){
      tree_markers_info *markers=fetch_treenode_precomputed_markers(trees,halo);
      tree_node_info    *marker =markers->half_peak_mass;
      if(marker!=NULL){
         int    halo_snap  =fetch_treenode_snap_tree(trees,halo);
         int    marker_snap=fetch_treenode_snap_tree(trees,marker);
         tau=10.*((trees->t_list[halo_snap]-trees->t_list[marker_snap])/trees->t_list[halo_snap]);
      }
   }
   return(tau);
}

