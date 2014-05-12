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

void init_treenode_list(const char          *catalog_name,
                        int                  n_list_alloc,
                        treenode_list_info **list){
  // Allocate the data structure
  (*list)=(treenode_list_info *)SID_malloc(sizeof(treenode_list_info));
  sprintf((*list)->catalog_name,"%s",catalog_name);
  (*list)->n_list          =0;
  (*list)->n_list_alloc    =n_list_alloc;
  (*list)->list            =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_list_alloc);
  (*list)->flag_groups_list=-1;
  (*list)->data            =NULL;
}

