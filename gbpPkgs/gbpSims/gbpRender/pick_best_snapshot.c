#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void pick_best_snap(double a_search,double *snap_a_list,int n_snap_a_list,int *snap_best,double *snap_diff_best_return){
   double snap_diff_best=1e60;
   for(int i_snap=0;i_snap<n_snap_a_list;i_snap++){
     double snap_diff=a_search-snap_a_list[i_snap];
     if(fabs(snap_diff)<fabs(snap_diff_best)){
       snap_diff_best=snap_diff;
       (*snap_best)  =i_snap;
     }
   }
   if(snap_diff_best_return!=NULL)
     (*snap_diff_best_return)=snap_diff_best;
}

