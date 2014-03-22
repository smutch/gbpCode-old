#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpTrees_build.h>

void free_trees_data(void **tree_data,void *params){
  if((*tree_data)!=NULL){
    void **tree_data_dealloc=(void **)(*tree_data);
    for(int i_snap=0;i_snap<((store_tree_data_free_parms_info *)params)->n_snaps;i_snap++)
       SID_free(SID_FARG tree_data_dealloc[i_snap]);
    SID_free(SID_FARG tree_data_dealloc);
  }
}

