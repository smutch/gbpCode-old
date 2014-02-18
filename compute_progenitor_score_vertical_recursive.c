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

void compute_progenitor_score_vertical_recursive(tree_vertical_node_info *tree,int *M_i,int mode){
  tree_vertical_node_info  *current;
  int                       i_progenitor;
  int                       M_iN,N_i,max_M_iN; // Defined in Section 2 of De Lucia and Blaizot (2006)

  N_i     =tree->halo.n_particles;
  max_M_iN=0;
  if(tree->n_progenitors>=1){
    current=tree->progenitor_first;
    i_progenitor=0;
    while(current!=NULL){
      M_iN=0;
      compute_progenitor_score_vertical_recursive(current,&M_iN,mode);
      if(check_mode_for_flag(mode,TREE_PROGENITOR_ORDER_DELUCIA))
        max_M_iN=MAX(max_M_iN,M_iN);
      i_progenitor++;
      current=current->progenitor_next;
    }
    if(i_progenitor!=tree->n_progenitors)
      SID_trap_error("There is a progenitor problem in compute_progenitor_score_recursive! (%d!=%d)",ERROR_LOGIC,i_progenitor!=tree->n_progenitors);
  }

  // Add this progenitor's score to the descendant's sum (see De Lucia and Blaizot (2006))
  (*M_i)=N_i+max_M_iN;
}

