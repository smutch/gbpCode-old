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

double fetch_treenode_R_vir(tree_info *trees,tree_node_info *halo){
   double R_vir=-1.;
   if(halo!=NULL){
      if(halo->parent_top==NULL){
         if(trees->group_properties!=NULL){
            halo_properties_info *properties=&(trees->group_properties[halo->snap_tree][halo->neighbour_index]);
            R_vir=(double)properties->R_vir;
         }
         else if(trees->group_properties_SHORT!=NULL){
            halo_properties_SHORT_info *properties=&(trees->group_properties_SHORT[halo->snap_tree][halo->neighbour_index]);
            R_vir=(double)properties->R_vir;
         }
         else
            SID_trap_error("Group properties are not defined in fetch_treenode_R_vir().  They probably have not been read.",ERROR_LOGIC);
      }
      else{
         if(trees->subgroup_properties!=NULL){
            halo_properties_info *properties=&(trees->subgroup_properties[halo->snap_tree][halo->neighbour_index]);
            R_vir=(double)properties->R_vir;
         }
         else if(trees->subgroup_properties_SHORT!=NULL){
            halo_properties_SHORT_info *properties=&(trees->subgroup_properties_SHORT[halo->snap_tree][halo->neighbour_index]);
            R_vir=(double)properties->R_vir;
         }
         else
            SID_trap_error("Subgroup properties are not defined in fetch_treenode_R_vir().  They probably have not been read.",ERROR_LOGIC);
      }
   }
   return(R_vir);
}

