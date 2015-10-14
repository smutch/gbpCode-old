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

double fetch_treenode_list_local_log_sigma_vx(tree_info *trees,treenode_list_info *list){
   double sigma_vx=0.;
   if(list->n_list_local>1){
      double vx_mean=0.;
      for(int i_halo=0;i_halo<list->n_list_local;i_halo++)
         vx_mean+=fetch_treenode_vx(trees,list->list[i_halo]);
      vx_mean/=(double)(list->n_list_local-1);
      for(int i_halo=0;i_halo<list->n_list_local;i_halo++)
         sigma_vx+=pow((fetch_treenode_vx(trees,list->list[i_halo])-vx_mean),2.);
      sigma_vx=sqrt(sigma_vx/(double)(list->n_list_local-1));
   }
   return(take_log10(sigma_vx));
}

