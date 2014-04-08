#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>
#include <gbpTrees_analysis.h>
#include <assert.h>

void write_treenode_list_hist(tree_info *trees,treenode_list_info *list){
  int     n_z_bin  = MAX(1,trees->n_snaps/25);
  double  logM_min = 9.00;
//  double  logM_min = 6.50;
  double  dlogM    = 0.25;
  int     n_logM   =   24;
  double  logN_min = 1.50;
  double  dlogN    = 0.25;
  int     n_logN   =   20;
  treenode_hist_info *hist_list_M;
  treenode_hist_info *hist_list_N;
  init_treenode_hist(trees,list->catalog_name,"z","N",TREENODE_HIST_LOG_Y,&hist_list_N,n_z_bin,logN_min,dlogN,n_logN);
  init_treenode_hist(trees,list->catalog_name,"z","M",TREENODE_HIST_LOG_Y,&hist_list_M,n_z_bin,logM_min,dlogM,n_logM);
  for(int i_list=0;i_list<list->n_list;i_list++){
     tree_node_info *current_halo=list->list[i_list];
     add_to_treenode_hist(trees,hist_list_N,current_halo);
     add_to_treenode_hist(trees,hist_list_M,current_halo);
  }
  write_treenode_hist(trees,hist_list_N);
  write_treenode_hist(trees,hist_list_M);
  free_treenode_hist (&hist_list_N);
  free_treenode_hist (&hist_list_M);
}

