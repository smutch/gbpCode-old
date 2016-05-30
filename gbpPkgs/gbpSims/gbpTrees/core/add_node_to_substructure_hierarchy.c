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

int add_node_to_substructure_hierarchy(tree_info       *trees,   // The tree datastructure
                                       tree_node_info  *subhalo, // Pointer to the substructure being added
                                       tree_node_info  *parent){ // Pointer to the parent being added to
  int rval=FALSE;
  // Set substructure pointers
  if(subhalo!=NULL){
     subhalo->parent=parent;
     if(parent!=NULL){
        parent->n_substructures++;
        if(parent->substructure_first==NULL)
           parent->substructure_first=subhalo;
        else
           parent->substructure_last->substructure_next=subhalo;
        parent->substructure_last=subhalo;
     }
     rval=TRUE;
  }
  return(rval);
}

