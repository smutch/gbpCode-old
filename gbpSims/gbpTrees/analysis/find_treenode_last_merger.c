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

int find_treenode_last_merger(tree_info       *trees,
                              tree_node_info  *halo,
                              double           fraction,
                              tree_node_info **remnant,
                              tree_node_info **merger_host,
                              tree_node_info **merged_halo){
   (*remnant)    =NULL;
   (*merger_host)=NULL;
   (*merged_halo)=NULL;
   if(halo!=NULL){
      tree_node_info *current_halo=halo;
      while(current_halo!=NULL && (*remnant)==NULL){
         if(current_halo->n_progenitors>1){
            double M_MP =0.; // MP =main progenitor
            double M_MMM=0.; // MMM=most massive merger
            tree_node_info *main_progenitor    =NULL;
            tree_node_info *most_massive_merger=NULL;
            tree_node_info *current_progenitor =current_halo->progenitor_first;
            while(current_progenitor!=NULL){
               halo_properties_info *properties=fetch_treenode_properties(trees,current_progenitor);
               if(properties->M_vir>M_MP){
                  M_MMM=M_MP;
                  M_MP =properties->M_vir;
                  most_massive_merger=main_progenitor;
                  main_progenitor    =current_progenitor;
               }
               else if(properties->M_vir>M_MP){
                  M_MMM              =properties->M_vir;
                  most_massive_merger=current_progenitor;
               }
               current_progenitor=current_progenitor->progenitor_next;
            }
            if((M_MMM/M_MP)>fraction){
               (*remnant)    =current_halo;
               (*merger_host)=main_progenitor;
               (*merged_halo)=most_massive_merger;
            }
         }
         current_halo=current_halo->progenitor_first;
      }
      return((*remnant)!=NULL);
   }
   return(FALSE);
}

