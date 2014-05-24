#include <stdio.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void split_forests_n_ways(int  *n_halos_forest,
                          int   n_forests,
                          int   n_split,
                          int **halo_count_split,
                          int **forest_lo_split,
                          int **forest_hi_split){
  int i_split;
  int i_forest;
  int n_halos_used;
  int n_halos;

  // Count the total bumber of halos
  for(i_forest=0,n_halos=0;i_forest<n_forests;i_forest++) 
     n_halos+=n_halos_forest[i_forest];

  // Alloc arrays
  (*halo_count_split)=(int *)SID_calloc(sizeof(int)*n_split);
  (*forest_lo_split) =(int *)SID_malloc(sizeof(int)*n_split);  
  (*forest_hi_split) =(int *)SID_malloc(sizeof(int)*n_split);  

  // Calculate the split
  for(i_split=0,n_halos_used=0,i_forest=0;i_split<n_split;i_split++,i_forest++){
    (*forest_lo_split)[i_split] =i_forest;
    (*forest_hi_split)[i_split] =i_forest;
    (*halo_count_split)[i_split]=0;
    if(i_forest<n_forests){
      int n_halos_target;
      n_halos_target               =(n_halos-n_halos_used)/(n_split-i_split);
      (*halo_count_split)[i_split]+=n_halos_forest[i_forest];
      while((*halo_count_split)[i_split]<n_halos_target && i_forest<(n_forests-1)){
        i_forest++;
        (*forest_hi_split)[i_split]  =i_forest;
        (*halo_count_split)[i_split]+=n_halos_forest[i_forest];
      }
      n_halos_used+=(*halo_count_split)[i_split];
    }
  }
  while(i_forest<(n_forests-1)){
    i_forest++;
    (*forest_hi_split)[n_split-1]  =i_forest;
    (*halo_count_split)[n_split-1]+=n_halos_forest[i_forest];
  }

}

