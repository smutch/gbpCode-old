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

double fetch_treenode_x_off(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL){
      halo_properties_info *properties;
      if(halo->parent_top==NULL){
         if(trees->group_properties!=NULL)
            properties=&(trees->group_properties[halo->snap_tree][halo->neighbour_index]);
         else
            SID_trap_error("Group properties are not defined in fetch_treenode_x_off().  They probably have not been read.",ERROR_LOGIC);
      }
      else{
         if(trees->subgroup_properties!=NULL)
            properties=&(trees->subgroup_properties[halo->snap_tree][halo->neighbour_index]);
         else
            SID_trap_error("Subgroup properties are not defined in fetch_treenode_x_off().  They probably have not been read.",ERROR_LOGIC);
      }
      double expansion_factor=a_of_z(trees->z_list[halo->snap_tree]);
      return(sqrt(pow(properties->position_COM[0]-properties->position_MBP[0],2.)+
                  pow(properties->position_COM[1]-properties->position_MBP[1],2.)+
                  pow(properties->position_COM[2]-properties->position_MBP[2],2.))/(properties->R_vir/expansion_factor));
   }
   return(1.);
}

