#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){
  char        filename_tree_in[256];
  int         n_trees;
  int         n_halos_total;
  int        *n_halos;
  int         i_tree;
  FILE       *fp;
  halo_info  *halos;
  halo_info   halo;
  int        *snap_num;
  size_t     *snap_num_index;
  int         i_snap,i_halo,j_halo,k_halo;
  int         n_halos_snap;
  int        *group_halo_first;
  int         group_halo_last;
  size_t     *group_halo_first_index;
  int        *snap_index;
  int descendant_min,descendant_max;
  int progenitor_first_min,progenitor_first_max;
  int progenitor_next_min,progenitor_next_max;
  int group_halo_first_min,group_halo_first_max;
  int group_halo_next_min,group_halo_next_max;
  int snap_num_min,snap_num_max;
  int halo_index_min,halo_index_max;
  int n_gal=0;
  int max_snap=0;
  int n_halos_max;
  int n_subtrees;
  double statistic;

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  strcpy(filename_tree_in,argv[1]);

  SID_log("Displaying tree stats for {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_tree_in);
  fp=fopen(filename_tree_in,"r");
  fread(&n_trees,      sizeof(int),1,fp);
  fread(&n_halos_total,sizeof(int),1,fp);
  SID_log("(%d trees and %d halos)...",SID_LOG_CONTINUE,n_trees,n_halos_total);
  n_halos=(int *)SID_malloc(sizeof(int)*n_trees);
  fread(n_halos,sizeof(int),n_trees,fp);
  calc_max(n_halos,&n_halos_max,n_trees,SID_INT,CALC_MODE_DEFAULT);
  SID_log("max_size=%d...",SID_LOG_CONTINUE,n_halos_max);
  halos               =(halo_info *)SID_malloc(sizeof(halo_info)*n_halos_max);
  snap_num            =(int       *)SID_malloc(sizeof(int)*n_halos_max);
  snap_index          =(int       *)SID_malloc(sizeof(int)*n_halos_max);
  group_halo_first    =(int       *)SID_malloc(sizeof(int)*n_halos_max);
  for(i_tree=0;i_tree<n_trees;i_tree++){
    fread(halos,sizeof(halo_info),n_halos[i_tree],fp);
    for(i_halo=0;i_halo<n_halos[i_tree];i_halo++) printf("%le\n",halos[i_halo].spin[0]);
/*
    printf("%5d %5d --",i_tree,n_halos[i_tree]);
    // Print tree tops
    for(i_halo=0,n_subtrees=0;i_halo<n_halos[i_tree];i_halo++){
      if(halos[i_halo].descendant<0) 
        printf(" %7d",halos[i_halo].halo_id);
    }
    // Print Mean spins
    for(i_halo=0,statistic=0.;i_halo<n_halos[i_tree];i_halo++)
      statistic+=(double)(halos[i_halo].spin[0]);
    statistic/=(double)n_halos[i_tree];
    printf(" %le",statistic);
*/
    printf("\n");
  }

  // Clean-up
  fclose(fp);
  SID_free((void **)&snap_num);
  SID_free((void **)&snap_index);
  SID_free((void **)&group_halo_first);
  SID_free((void **)&n_halos);
  SID_free((void **)&halos);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(0);
}
