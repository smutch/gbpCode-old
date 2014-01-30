#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void assign_progenitor_order_recursive(tree_node_info *tree,int *M_i,int mode){
  tree_node_info  *first_new;
  tree_node_info  *last_new;
  tree_node_info  *current;
  tree_node_info **progenitors;
  int              i_progenitor;
  size_t          *M_iN_index;
  size_t          *M_iN_rank; 
  int             *M_iN,N_i,max_M_iN; // Defined in Section 2 of De Lucia and Blaizot (2006)

  N_i     =tree->halo.n_particles;
  max_M_iN=0;
  if(tree->n_progenitors>=1){
    // Sum all progenitor contributions to score
    M_iN        =(int             *)SID_malloc(sizeof(int)*tree->n_progenitors);
    progenitors =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*tree->n_progenitors);
    i_progenitor=0;
    current     =tree->progenitor_first;
    while(current!=NULL){
      M_iN[i_progenitor]=0;
      assign_progenitor_order_recursive(current,&(M_iN[i_progenitor]),mode);
      progenitors[i_progenitor]=current;
      if(check_mode_for_flag(mode,TREE_PROGENITOR_ORDER_DELUCIA))
        max_M_iN=MAX(max_M_iN,M_iN[i_progenitor]);
      i_progenitor++;
      current=current->progenitor_next;
    }
    if(i_progenitor!=tree->n_progenitors)
      SID_trap_error("There is a progenitor problem in assign_progenitor_order_recursive! (%d!=%d)",ERROR_LOGIC,i_progenitor!=tree->n_progenitors);

    // Assign progenitors by (descending) order of their score
    //   (note: sorting the sort indicies gives each progenitor's ranking in the sort)
    merge_sort(M_iN,      (size_t)tree->n_progenitors,&M_iN_index,SID_INT,   SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
    merge_sort(M_iN_index,(size_t)tree->n_progenitors,&M_iN_rank, SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
    first_new=progenitors[M_iN_index[tree->n_progenitors-1]]; 
    last_new =progenitors[M_iN_index[0]]; 
    tree->progenitor_first=first_new;
    tree->progenitor_last =last_new;
    for(i_progenitor=0;i_progenitor<tree->n_progenitors;i_progenitor++){
      if(progenitors[i_progenitor]!=last_new)
        progenitors[i_progenitor]->progenitor_next=progenitors[M_iN_index[M_iN_rank[i_progenitor]-1]]; 
      else
        progenitors[i_progenitor]->progenitor_next=NULL;
    }

    SID_free(SID_FARG progenitors);
    SID_free(SID_FARG M_iN);
    SID_free(SID_FARG M_iN_index);
    SID_free(SID_FARG M_iN_rank);
  }

  // Add this progenitor's score to the descendant's sum (see De Lucia and Blaizot (2006))
  (*M_i)=N_i+max_M_iN;
}

