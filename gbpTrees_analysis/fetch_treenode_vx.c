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

double fetch_treenode_vx(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL){
      halo_properties_info *properties;
      if(halo->parent==NULL){
         if(trees->group_properties!=NULL)
            properties=&(trees->group_properties[halo->snap_tree][halo->neighbour_index]);
         else
            SID_trap_error("Group properties are not defined.  They probably have not been read.",ERROR_LOGIC);
      }
      else{
         if(trees->subgroup_properties!=NULL)
            properties=&(trees->subgroup_properties[halo->snap_tree][halo->neighbour_index]);
         else
            SID_trap_error("Subgroup properties are not defined.  They probably have not been read.",ERROR_LOGIC);
      }
      return(properties->velocity_COM[0]);
   }
   return(-1.);
}

