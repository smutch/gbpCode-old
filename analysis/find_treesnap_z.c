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

int find_treesnap_z(tree_info *trees,double z_obs_exact){
  int    i_z_obs=0;
  double delta_z;
  double delta_z_min=1e6;
  for(int i_z=0;i_z<trees->n_snaps;i_z++){
     delta_z=fabs(trees->z_list[i_z]-z_obs_exact);
     if(delta_z<delta_z_min){
        i_z_obs    =i_z;
        delta_z_min=delta_z;
     }
  }
  return(i_z_obs);
}

