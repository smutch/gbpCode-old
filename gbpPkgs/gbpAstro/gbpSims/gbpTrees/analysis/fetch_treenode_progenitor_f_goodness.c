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

float fetch_treenode_progenitor_f_goodness(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL)
      return(fetch_treenode_descendant_f_goodness(trees,halo->progenitor_first));
   else
      return(-1.);
}

