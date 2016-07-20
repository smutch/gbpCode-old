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
      int flag_found=find_treenode_main_progenitor(trees,halo,&main_progenitor);
      if(main_progenitor!=NULL){
         // Compute the displacement (node that dx,dy,dz are comoving while R_vir is physical ... hence the factor of 'a')
         double x    =fetch_treenode_x(trees,halo);
         double y    =fetch_treenode_y(trees,halo);
         double z    =fetch_treenode_z(trees,halo);
         double x_MP =fetch_treenode_x(trees,main_progenitor);
         double y_MP =fetch_treenode_y(trees,main_progenitor);
         double z_MP =fetch_treenode_z(trees,main_progenitor);
         double R_vir=fetch_treenode_R_vir(trees,main_progenitor);
         double dx   =d_periodic(x-x_MP,(double)trees->box_size);
         double dy   =d_periodic(y-y_MP,(double)trees->box_size);
         double dz   =d_periodic(z-z_MP,(double)trees->box_size);
         return(a_of_z(trees->z_list[halo->snap_tree])*sqrt(dx*dx+dy*dy+dz*dz)/R_vir);
      }
      // In this case, a halo is it's own main progenitor
      else if(flag_found)
         return(0.);
   }
   return(-1.);
}

