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

void reset_treenode_list(treenode_list_info *list){
  for(int i_node=0;i_node<list->n_list_alloc;i_node++)
     list->list[i_node]=NULL;
  list->n_list      =0;
  list->n_list_local=0;
}

