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

int construct_unique_tree_id(tree_node_info *tree_node,int tree_number){
  if(tree_node!=NULL)
    return(tree_node->depth_first_index);
  else
    return(-1);    
}

