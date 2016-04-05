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

double fetch_treenode_vy(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL){
      double vy;
      if(halo->parent==NULL){
         if(trees->group_properties!=NULL){
            halo_properties_info *properties=&(trees->group_properties[halo->snap_tree][halo->neighbour_index]);
            vy=properties->velocity_COM[1];
         }
         else if(trees->group_properties_SHORT!=NULL){
            halo_properties_SHORT_info *properties=&(trees->group_properties_SHORT[halo->snap_tree][halo->neighbour_index]);
            vy=properties->vel[1];
         }
         else if(trees->group_properties_SAGE!=NULL){
            halo_properties_SAGE_info *properties=&(trees->group_properties_SAGE[halo->snap_tree][halo->neighbour_index]);
            vy=properties->vel[1];
         }
         else
            SID_trap_error("Group properties are not defined in fetch_treenode_vy().  They probably have not been read.",ERROR_LOGIC);
      }
      else{
         if(trees->subgroup_properties!=NULL){
            halo_properties_info *properties=&(trees->subgroup_properties[halo->snap_tree][halo->neighbour_index]);
            vy=properties->velocity_COM[1];
         }
         else if(trees->subgroup_properties_SHORT!=NULL){
            halo_properties_SHORT_info *properties=&(trees->subgroup_properties_SHORT[halo->snap_tree][halo->neighbour_index]);
            vy=properties->vel[1];
         }
         else if(trees->subgroup_properties_SAGE!=NULL){
            halo_properties_SAGE_info *properties=&(trees->subgroup_properties_SAGE[halo->snap_tree][halo->neighbour_index]);
            vy=properties->vel[1];
         }
         else
            SID_trap_error("Subgroup properties are not defined in fetch_treenode_vy().  They probably have not been read.",ERROR_LOGIC);
      }
      return(vy);
   }
   return(-1.);
}

