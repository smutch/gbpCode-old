#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>
#include <assert.h>

void write_treenode_all_hist(tree_info *trees,const char *filename_out_root_in,int i_type,double logM_min,double dlogM,int n_logM){
  SID_log("Writing treenode histograms...",SID_LOG_OPEN|SID_LOG_TIMER);
  int     n_z_bin  = MAX(1,trees->n_snaps/25);
  //double  logM_min = 9.00;
  //double  logM_min = 6.50;
  //double  dlogM    = 0.25;
  //int     n_logM   =   24;
  double  logN_min = 1.50;
  double  dlogN    = 0.25;
  int     n_logN   =   20;

  char group_prefix[8];
  if(i_type==0)
     sprintf(group_prefix,"");
  else
     sprintf(group_prefix,"sub");

  treenode_hist_info *hist_all_N;
  treenode_hist_info *hist_all_Npeak;
  treenode_hist_info *hist_all_M;
  treenode_hist_info *hist_all_Mpeak;
  treenode_hist_info *hist_centrals_N;
  treenode_hist_info *hist_centrals_Npeak;
  treenode_hist_info *hist_centrals_M;
  treenode_hist_info *hist_centrals_Mpeak;
  treenode_hist_info *hist_substructure_N;
  treenode_hist_info *hist_substructure_Npeak;
  treenode_hist_info *hist_substructure_M;
  treenode_hist_info *hist_substructure_Mpeak;

  char filename_out_root[MAX_FILENAME_LENGTH];
  sprintf(filename_out_root,"%s_%sgroup",filename_out_root_in,group_prefix);

  init_treenode_hist(trees,"all","z","N",     TREENODE_HIST_LOG_Y,&hist_all_N,    n_z_bin,logN_min,dlogN,n_logN);
  init_treenode_hist(trees,"all","z","N_peak",TREENODE_HIST_LOG_Y,&hist_all_Npeak,n_z_bin,logN_min,dlogN,n_logN);
  init_treenode_hist(trees,"all","z","M",     TREENODE_HIST_LOG_Y,&hist_all_M,    n_z_bin,logM_min,dlogM,n_logM);
  init_treenode_hist(trees,"all","z","M_peak",TREENODE_HIST_LOG_Y,&hist_all_Mpeak,n_z_bin,logM_min,dlogM,n_logM);
  if(i_type==1){
     init_treenode_hist(trees,"centrals",    "z","N",     TREENODE_HIST_LOG_Y,&hist_centrals_N,    n_z_bin,logN_min,dlogN,n_logN);
     init_treenode_hist(trees,"centrals",    "z","N_peak",TREENODE_HIST_LOG_Y,&hist_centrals_Npeak,n_z_bin,logN_min,dlogN,n_logN);
     init_treenode_hist(trees,"centrals",    "z","M",     TREENODE_HIST_LOG_Y,&hist_centrals_M,    n_z_bin,logM_min,dlogM,n_logM);
     init_treenode_hist(trees,"centrals",    "z","M_peak",TREENODE_HIST_LOG_Y,&hist_centrals_Mpeak,n_z_bin,logM_min,dlogM,n_logM);
     init_treenode_hist(trees,"substructure","z","N",     TREENODE_HIST_LOG_Y,&hist_substructure_N,    n_z_bin,logN_min,dlogN,n_logN);
     init_treenode_hist(trees,"substructure","z","N_peak",TREENODE_HIST_LOG_Y,&hist_substructure_Npeak,n_z_bin,logN_min,dlogN,n_logN);
     init_treenode_hist(trees,"substructure","z","M",     TREENODE_HIST_LOG_Y,&hist_substructure_M,    n_z_bin,logM_min,dlogM,n_logM);
     init_treenode_hist(trees,"substructure","z","M_peak",TREENODE_HIST_LOG_Y,&hist_substructure_Mpeak,n_z_bin,logM_min,dlogM,n_logM);
  }

  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo;
     if(i_type==0)
        current_halo=trees->first_neighbour_groups[i_snap];
     else
        current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        if(check_treenode_if_central(current_halo)){
           add_to_treenode_hist(trees,hist_centrals_N,    current_halo);
           add_to_treenode_hist(trees,hist_centrals_Npeak,current_halo);
           add_to_treenode_hist(trees,hist_centrals_M,    current_halo);
           add_to_treenode_hist(trees,hist_centrals_Mpeak,current_halo);
        }
        else if(check_treenode_if_satellite(current_halo)){
           add_to_treenode_hist(trees,hist_substructure_N,    current_halo);
           add_to_treenode_hist(trees,hist_substructure_Npeak,current_halo);
           add_to_treenode_hist(trees,hist_substructure_M,    current_halo);
           add_to_treenode_hist(trees,hist_substructure_Mpeak,current_halo);
        }
        add_to_treenode_hist(trees,hist_all_N,    current_halo);
        add_to_treenode_hist(trees,hist_all_Npeak,current_halo);
        add_to_treenode_hist(trees,hist_all_M,    current_halo);
        add_to_treenode_hist(trees,hist_all_Mpeak,current_halo);
        current_halo=current_halo->next_neighbour;
     }
  }

  SID_Allreduce(SID_IN_PLACE,(hist_all_N->array),    (hist_all_N->n_x*hist_all_N->n_y),        SID_INT,SID_SUM,SID.COMM_WORLD);
  SID_Allreduce(SID_IN_PLACE,(hist_all_Npeak->array),(hist_all_Npeak->n_x*hist_all_Npeak->n_y),SID_INT,SID_SUM,SID.COMM_WORLD);
  SID_Allreduce(SID_IN_PLACE,(hist_all_M->array),    (hist_all_M->n_x*hist_all_M->n_y),        SID_INT,SID_SUM,SID.COMM_WORLD);
  SID_Allreduce(SID_IN_PLACE,(hist_all_Mpeak->array),(hist_all_Mpeak->n_x*hist_all_Mpeak->n_y),SID_INT,SID_SUM,SID.COMM_WORLD);
  if(i_type==1){
     SID_Allreduce(SID_IN_PLACE,(hist_centrals_N->array),        (hist_centrals_N->n_x*hist_centrals_N->n_y),                SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,(hist_centrals_Npeak->array),    (hist_centrals_Npeak->n_x*hist_centrals_Npeak->n_y),        SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,(hist_centrals_M->array),        (hist_centrals_M->n_x*hist_centrals_M->n_y),                SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,(hist_centrals_Mpeak->array),    (hist_centrals_Mpeak->n_x*hist_centrals_Mpeak->n_y),        SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,(hist_substructure_N->array),    (hist_substructure_N->n_x*hist_substructure_N->n_y),        SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,(hist_substructure_Npeak->array),(hist_substructure_Npeak->n_x*hist_substructure_Npeak->n_y),SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,(hist_substructure_M->array),    (hist_substructure_M->n_x*hist_substructure_M->n_y),        SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,(hist_substructure_Mpeak->array),(hist_substructure_Mpeak->n_x*hist_substructure_Mpeak->n_y),SID_INT,SID_SUM,SID.COMM_WORLD);
  }

  write_treenode_hist(trees,filename_out_root,hist_all_N);
  write_treenode_hist(trees,filename_out_root,hist_all_Npeak);
  write_treenode_hist(trees,filename_out_root,hist_all_M);
  write_treenode_hist(trees,filename_out_root,hist_all_Mpeak);
  if(i_type==1){
     write_treenode_hist(trees,filename_out_root,hist_centrals_N);
     write_treenode_hist(trees,filename_out_root,hist_centrals_Npeak);
     write_treenode_hist(trees,filename_out_root,hist_centrals_M);
     write_treenode_hist(trees,filename_out_root,hist_centrals_Mpeak);
     write_treenode_hist(trees,filename_out_root,hist_substructure_N);
     write_treenode_hist(trees,filename_out_root,hist_substructure_Npeak);
     write_treenode_hist(trees,filename_out_root,hist_substructure_M);
     write_treenode_hist(trees,filename_out_root,hist_substructure_Mpeak);
  }

  free_treenode_hist (&hist_all_N);
  free_treenode_hist (&hist_all_Npeak);
  free_treenode_hist (&hist_all_M);
  free_treenode_hist (&hist_all_Mpeak);
  if(i_type==1){
     free_treenode_hist (&hist_centrals_N);
     free_treenode_hist (&hist_centrals_Npeak);
     free_treenode_hist (&hist_centrals_M);
     free_treenode_hist (&hist_centrals_Mpeak);
     free_treenode_hist (&hist_substructure_N);
     free_treenode_hist (&hist_substructure_Npeak);
     free_treenode_hist (&hist_substructure_M);
     free_treenode_hist (&hist_substructure_Mpeak);
  }

  SID_log("Done.",SID_LOG_CLOSE);
}

