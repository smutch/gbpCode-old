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

double fetch_treenode_Mpeak(tree_info *trees,tree_node_info *halo){
   tree_markers_info *markers  =fetch_treenode_precomputed_markers(trees,halo);
   tree_node_info    *peak_mass=NULL;
   if(markers!=NULL)
      peak_mass=markers->peak_mass;
   if(peak_mass==NULL)
      return(fetch_treenode_Mvir(trees,halo));
   else
      return(fetch_treenode_Mvir(trees,peak_mass));
}

