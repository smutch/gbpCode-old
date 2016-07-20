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
#include <gbpTrees_analysis.h>

void add_to_treenode_list(treenode_list_info *list,tree_node_info *node){

  if(node!=NULL){
     // If this is the first group we've added, then initialize
     //   the flag for making sure that we are including only one
     //   halo type in this list
     int flag_add_a_group=(node->parent_top==NULL);
     if(list->flag_groups_list<0)
        list->flag_groups_list=flag_add_a_group;

     // Sanity check
     if(flag_add_a_group!=list->flag_groups_list)
        SID_trap_error("You are attempting to mix groups and subgroups in a treenode_list structure.  Not allowed.",ERROR_LOGIC);

     // Increment the counter
     list->n_list_local++;

     // Reallocate the array if the list has grown too big
     if(list->n_list_local>list->n_list_alloc){
       if(list->data!=NULL)
         SID_trap_error("You are attempting to grow a treenode list structure with added data past it's allocated size.  Not allowed.",ERROR_LOGIC); 
       list->n_list_alloc*=2;
       list->list=(tree_node_info **)SID_realloc(list->list,(sizeof(tree_node_info *)*list->n_list_alloc));
     }

     // Add the node
     list->list[list->n_list_local-1]=node;
  }
}

