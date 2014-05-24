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

void write_treenode_all_hist(tree_info *trees,const char *filename_out_root){
  SID_log("Writing treenode histograms...",SID_LOG_OPEN|SID_LOG_TIMER);
  int     n_z_bin  = MAX(1,trees->n_snaps/25);
  double  logM_min = 9.00;
  //double  logM_min = 6.50;
  double  dlogM    = 0.25;
  int     n_logM   =   24;
  double  logN_min = 1.50;
  double  dlogN    = 0.25;
  int     n_logN   =   20;
  treenode_hist_info *hist_all_M;
  treenode_hist_info *hist_all_N;
  treenode_hist_info *hist_centrals_M;
  treenode_hist_info *hist_centrals_N;
  treenode_hist_info *hist_substructure_M;
  treenode_hist_info *hist_substructure_N;
  init_treenode_hist(trees,"all",         "z","M",TREENODE_HIST_LOG_Y,&hist_all_M,         n_z_bin,logM_min,dlogM,n_logM);
  init_treenode_hist(trees,"all",         "z","N",TREENODE_HIST_LOG_Y,&hist_all_N,         n_z_bin,logN_min,dlogN,n_logN);
  init_treenode_hist(trees,"centrals",    "z","M",TREENODE_HIST_LOG_Y,&hist_centrals_M,    n_z_bin,logM_min,dlogM,n_logM);
  init_treenode_hist(trees,"centrals",    "z","N",TREENODE_HIST_LOG_Y,&hist_centrals_N,    n_z_bin,logN_min,dlogN,n_logN);
  init_treenode_hist(trees,"substructure","z","M",TREENODE_HIST_LOG_Y,&hist_substructure_M,n_z_bin,logM_min,dlogM,n_logM);
  init_treenode_hist(trees,"substructure","z","N",TREENODE_HIST_LOG_Y,&hist_substructure_N,n_z_bin,logN_min,dlogN,n_logN);
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        if(check_treenode_if_central(current_halo)){
           add_to_treenode_hist(trees,hist_centrals_M,current_halo);
           add_to_treenode_hist(trees,hist_centrals_N,current_halo);
        }
        else{
           add_to_treenode_hist(trees,hist_substructure_M,current_halo);
           add_to_treenode_hist(trees,hist_substructure_N,current_halo);
        }
        add_to_treenode_hist(trees,hist_all_M,current_halo);
        add_to_treenode_hist(trees,hist_all_N,current_halo);
        current_halo=current_halo->next_neighbour;
     }
  }
  write_treenode_hist(trees,filename_out_root,hist_all_N);
  write_treenode_hist(trees,filename_out_root,hist_all_M);
  write_treenode_hist(trees,filename_out_root,hist_centrals_N);
  write_treenode_hist(trees,filename_out_root,hist_centrals_M);
  write_treenode_hist(trees,filename_out_root,hist_substructure_N);
  write_treenode_hist(trees,filename_out_root,hist_substructure_M);
  free_treenode_hist (&hist_all_N);
  free_treenode_hist (&hist_all_M);
  free_treenode_hist (&hist_centrals_N);
  free_treenode_hist (&hist_centrals_M);
  free_treenode_hist (&hist_substructure_N);
  free_treenode_hist (&hist_substructure_M);
  SID_log("Done.",SID_LOG_CLOSE);
}

