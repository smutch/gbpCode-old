#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

double fetch_treenode_delta_r_MP(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL){
      tree_node_info *main_progenitor=NULL;
      find_treenode_main_progenitor(trees,halo,&main_progenitor);
      if(main_progenitor!=NULL){
         // Set the halo properties
         halo_properties_info *properties_halo=fetch_treenode_properties(trees,halo);
         halo_properties_info *properties_MP  =fetch_treenode_properties(trees,main_progenitor);
         // Compute the displacement (node that dx,dy,dz are comoving while R_vir is physical ... hence the factor of 'a')
         double dx=d_periodic((double)properties_halo->position_MBP[0]-(double)properties_MP->position_MBP[0],(double)trees->box_size);
         double dy=d_periodic((double)properties_halo->position_MBP[1]-(double)properties_MP->position_MBP[1],(double)trees->box_size);
         double dz=d_periodic((double)properties_halo->position_MBP[2]-(double)properties_MP->position_MBP[2],(double)trees->box_size);
         return(a_of_z(trees->z_list[halo->snap_tree])*sqrt(dx*dx+dy*dy+dz*dz)/properties_MP->R_vir);
      }
   }
   return(-1.);
}

