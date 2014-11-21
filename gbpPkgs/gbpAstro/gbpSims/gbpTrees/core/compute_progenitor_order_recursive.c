#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void compute_progenitor_order_recursive(tree_info *trees,tree_node_info *descendant,int *score_descendant,int mode){
  if(descendant==NULL)
     SID_trap_error("NULL descendant in compute_progenitor_order_recursive().",ERROR_LOGIC);
  int N_i     =descendant->n_particles;
  int max_M_iN=0;
  if(descendant->n_progenitors>1){
     tree_node_info **node_list=NULL;
     int             *score    =NULL;  
     node_list=(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*descendant->n_progenitors);
     score    =(int             *)SID_malloc(sizeof(int)             *descendant->n_progenitors);

     // Loop over each progenitor (recursively)
     tree_node_info *current_progenitor;
     size_t n_halos=0;
     current_progenitor=descendant->progenitor_first;
     while(current_progenitor!=NULL){
        if(n_halos>=descendant->n_progenitors)
           SID_trap_error("Progenitor count exceeded in compute_progenitor_order_recursive().",ERROR_LOGIC);
        node_list[n_halos]=current_progenitor;                                         // Make node  list
        compute_progenitor_order_recursive(trees,current_progenitor,&(score[n_halos]),mode); // Make score list
        if(check_mode_for_flag(mode,TREE_PROGENITOR_ORDER_DELUCIA))
           max_M_iN=MAX(max_M_iN,score[n_halos]);
        n_halos++;
        current_progenitor=current_progenitor->progenitor_next;
     }
     if(n_halos!=descendant->n_progenitors)
        SID_trap_error("Invalid halo count (%d!=%d) in compute_progenitor_order_score_recursive().",ERROR_LOGIC,
                       n_halos,descendant->n_progenitors);

     // Sort the resulting scores (ascending!) ...
     size_t *score_index=NULL;
     merge_sort(score,n_halos,&score_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

     // Reorder the pointers in *decending* order
     int             i_halo;
     tree_node_info *new_first;
     tree_node_info *new_last;
     new_first=node_list[score_index[n_halos-1]];
     new_last =node_list[score_index[0]];
     for(i_halo=0;i_halo<n_halos;i_halo++){
       if(i_halo==0)
         node_list[score_index[i_halo]]->progenitor_next=NULL;
       else
         node_list[score_index[i_halo]]->progenitor_next=node_list[score_index[i_halo-1]]; // descending
     }
     descendant->progenitor_first=new_first;
     descendant->progenitor_last =new_last;

     if(check_mode_for_flag(mode,TREE_PROGENITOR_ORDER_N_PARTICLES_PEAK))
        max_M_iN=score[score_index[n_halos-1]];

     // Clean-up 
     SID_free(SID_FARG score_index);
     SID_free(SID_FARG score);
     SID_free(SID_FARG node_list);
  }
  else if(descendant->progenitor_first!=NULL){
     compute_progenitor_order_recursive(trees,descendant->progenitor_first,score_descendant,mode);
     if(score_descendant!=NULL){
        if(check_mode_for_flag(mode,TREE_PROGENITOR_ORDER_DELUCIA))
           max_M_iN=(*score_descendant);
        else if(check_mode_for_flag(mode,TREE_PROGENITOR_ORDER_N_PARTICLES_PEAK))
           max_M_iN=(*score_descendant);
     }
  }

  // Pass order scores up the heirarchy
  if(score_descendant!=NULL){
     if(check_mode_for_flag(mode,TREE_PROGENITOR_ORDER_DELUCIA))
        // Add this progenitor's score to the descendant's sum (see De Lucia and Blaizot (2006))
        (*score_descendant)=N_i+max_M_iN;
     else if(check_mode_for_flag(mode,TREE_PROGENITOR_ORDER_N_PARTICLES))
        (*score_descendant)=N_i;
     else if(check_mode_for_flag(mode,TREE_PROGENITOR_ORDER_N_PARTICLES_PEAK))
        (*score_descendant)=MAX(N_i,max_M_iN);
     else
        SID_trap_error("Invalid mode (%d) in assign_progenitor_order_recursive().",ERROR_LOGIC,mode);
   }

}

