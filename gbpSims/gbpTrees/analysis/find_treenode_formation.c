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

int find_treenode_formation(tree_info       *trees,
                            tree_node_info  *halo,
                            double           fraction,
                            tree_node_info **peak_mass,
                            tree_node_info **fraction_of_peak_mass){
   if(halo!=NULL){
      double M_vir   =0.;
      double M_peak  =0.;
      double M_target=0.;
      tree_node_info *current_halo=halo;
      (*peak_mass)                =current_halo;
      (*fraction_of_peak_mass)    =current_halo;
      while(current_halo!=NULL){
         halo_properties_info *properties=fetch_treenode_properties(trees,current_halo);
         M_vir=properties->M_vir;
         if(M_vir>M_peak){
            (*peak_mass)=current_halo;
            M_peak      =M_vir;
            M_target    =fraction*M_peak;
         }
         if(M_vir>M_target)
            (*fraction_of_peak_mass)=current_halo;
         current_halo=current_halo->progenitor_first;
      }
      return(TRUE);
   }
   (*fraction_of_peak_mass)=NULL;
   return(FALSE);
}

