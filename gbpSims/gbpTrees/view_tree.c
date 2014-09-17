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
  int         select_tree;
  int         n_trees;
  int         n_halos_total;
  int        *n_halos;
  int         i_tree;
  FILE       *fp;
  halo_properties_SAGE_info  *halos;
  halo_properties_SAGE_info   halo;
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

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  strcpy(filename_tree_in,argv[1]);
  select_tree=atoi(argv[2]);

  SID_log("Displaying tree %d from {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,select_tree,filename_tree_in);
  fp=fopen(filename_tree_in,"r");
  fread(&n_trees,      sizeof(int),1,fp);
  fread(&n_halos_total,sizeof(int),1,fp);
  SID_log("(%d trees and %d halos)...",SID_LOG_CONTINUE,n_trees,n_halos_total);
  n_halos=(int *)SID_malloc(sizeof(int)*n_trees);
  fread(n_halos,sizeof(int),n_trees,fp);
  for(i_tree=0;i_tree<select_tree;i_tree++){
    for(i_halo=0;i_halo<n_halos[i_tree];i_halo++){
      fread(&halo,sizeof(halo_properties_SAGE_info),1,fp);
      max_snap=MAX(max_snap,halo.snap_num);
    }
  }
  halos               =(halo_properties_SAGE_info *)SID_malloc(sizeof(halo_properties_SAGE_info)*n_halos[i_tree]);
  snap_num            =(int       *)SID_malloc(sizeof(int)*n_halos[i_tree]);
  snap_index          =(int       *)SID_malloc(sizeof(int)*n_halos[i_tree]);
  group_halo_first    =(int       *)SID_malloc(sizeof(int)*n_halos[i_tree]);
  fread(halos,sizeof(halo_properties_SAGE_info),n_halos[i_tree],fp);
  descendant_min      =10000;
  descendant_max      =    0;
  progenitor_first_min=10000;
  progenitor_first_max=    0;
  progenitor_next_min =10000;
  progenitor_next_max =    0;
  group_halo_first_min=10000;
  group_halo_first_max=    0;
  group_halo_next_min =10000;
  group_halo_next_max =    0;
  snap_num_min        =10000;
  snap_num_max        =    0;
  halo_index_min      =10000;
  halo_index_max      =    0; 
  for(i_halo=0;i_halo<n_halos[i_tree];i_halo++){
    snap_num[i_halo] =halos[i_halo].snap_num;

    if(halos[i_halo].descendant>=0)
      descendant_min      =MIN(descendant_min,halos[i_halo].descendant);
    if(halos[i_halo].progenitor_first>=0)
      progenitor_first_min=MIN(progenitor_first_min,halos[i_halo].progenitor_first);
    if(halos[i_halo].progenitor_next>=0)
      progenitor_next_min =MIN(progenitor_next_min,halos[i_halo].progenitor_next);
    if(halos[i_halo].group_halo_first>=0)
      group_halo_first_min=MIN(group_halo_first_min,halos[i_halo].group_halo_first);
    if(halos[i_halo].group_halo_next>=0)
      group_halo_next_min =MIN(group_halo_next_min,halos[i_halo].group_halo_next);
    if(halo.snap_num>=0)
      snap_num_min        =MIN(snap_num_min,halos[i_halo].snap_num);
    if(halos[i_halo].halo_index>=0)
      halo_index_min      =MIN(halo_index_min,halos[i_halo].halo_index);

    descendant_max      =MAX(descendant_max,halos[i_halo].descendant);
    progenitor_first_max=MAX(progenitor_first_max,halos[i_halo].progenitor_first);
    progenitor_next_max =MAX(progenitor_next_max,halos[i_halo].progenitor_next);
    group_halo_first_max=MAX(group_halo_first_max,halos[i_halo].group_halo_first);
    group_halo_next_max =MAX(group_halo_next_max,halos[i_halo].group_halo_next);
    snap_num_max        =MAX(snap_num_max,halos[i_halo].snap_num);
    halo_index_max      =MAX(halo_index_max,halos[i_halo].halo_index);
    max_snap=MAX(max_snap,halos[i_halo].snap_num);
  }

  i_tree++;
  for(;i_tree<n_trees;i_tree++){
    for(i_halo=0;i_halo<n_halos[i_tree];i_halo++){
      fread(&halo,sizeof(halo_properties_SAGE_info),1,fp);
      max_snap=MAX(max_snap,halo.snap_num);
    }
  }

  rewind(fp); 
  fread(&n_trees,      sizeof(int),1,fp);
  fread(&n_halos_total,sizeof(int),1,fp); 
  for(i_tree=0,n_gal=0;i_tree<n_trees;i_tree++){
    for(i_halo=0;i_halo<n_halos[i_tree];i_halo++){
      fread(&halo,sizeof(halo_properties_SAGE_info),1,fp);
      if(halo.snap_num==max_snap) n_gal++;
    }
  }

  SID_log("n_trees          =%d",    SID_LOG_COMMENT,n_trees);
  SID_log("n_halos[snap=%3d]=%d",    SID_LOG_COMMENT,n_gal);
  SID_log("n_halos_total    =%d",    SID_LOG_COMMENT,n_halos_total);
  SID_log("Descendants      =%d->%d",SID_LOG_COMMENT,descendant_min,      descendant_max);
  SID_log("Progenitor_first =%d->%d",SID_LOG_COMMENT,progenitor_first_min,progenitor_first_max);
  SID_log("Progenitor_next  =%d->%d",SID_LOG_COMMENT,progenitor_next_min, progenitor_next_max);
  SID_log("Group_halo_first =%d->%d",SID_LOG_COMMENT,group_halo_first_min,group_halo_first_max);
  SID_log("Group_halo_next  =%d->%d",SID_LOG_COMMENT,group_halo_next_min, group_halo_next_max);
  SID_log("Snap_num         =%d->%d",SID_LOG_COMMENT,snap_num_min,        snap_num_max);
  SID_log("Halo_index       =%d->%d",SID_LOG_COMMENT,halo_index_min,      halo_index_max);

  merge_sort((void *)snap_num,(size_t)n_halos[i_tree],&snap_num_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);
  for(i_snap=snap_num_max,i_halo=n_halos[i_tree]-1;i_snap>=snap_num_min && i_halo>=0;i_snap--){
    n_halos_snap=0;
    while(snap_num[snap_num_index[i_halo]]==i_snap && i_halo>0){
      n_halos_snap++;
      i_halo--;
    }
    if(snap_num[snap_num_index[i_halo]]==i_snap){
      n_halos_snap++;
      i_halo--;
    }
    for(j_halo=0;j_halo<n_halos_snap;j_halo++){
      group_halo_first[j_halo]=halos[snap_num_index[i_halo+j_halo+1]].group_halo_first;
      snap_index[j_halo]      =snap_num_index[i_halo+j_halo+1];
    }
    merge_sort((void *)group_halo_first,(size_t)n_halos_snap,&group_halo_first_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);
    group_halo_last=-99;
    if(n_halos_snap>0)
      printf("Snap #%3d: ",i_snap);
    for(j_halo=0;j_halo<n_halos_snap;j_halo++){
      k_halo=snap_index[group_halo_first_index[j_halo]];
      if(group_halo_last!=halos[k_halo].group_halo_first){
        if(j_halo!=0)
          printf(") ");        
        printf("(");        
      }
      else
        printf(" ");        
      // Generate output
/*
      printf("I:%5d->%5d/S:%5d/PF:%5d/PN:%5d/G:%5d->%5d",
             k_halo,
             halos[k_halo].descendant,
             halos[k_halo].n_particles,
             halos[k_halo].progenitor_first,
             halos[k_halo].progenitor_next,
             halos[k_halo].group_halo_first,
             halos[k_halo].group_halo_next);      
*/
/*
           printf("I%05d.D%05d.S%05d.F%05d.N%05d.f%05d.n%05d",
                  k_halo,
                  halos[k_halo].descendant,
                  halos[k_halo].n_particles,
                  halos[k_halo].progenitor_first,
                  halos[k_halo].progenitor_next,
                  halos[k_halo].group_halo_first,
                  halos[k_halo].group_halo_next);      
*/
/*
           printf("%05d",
                  halos[k_halo].n_particles);
*/
/*
           printf("%05d/%05d/%05d",
                  k_halo,halos[k_halo].descendant,
                  halos[k_halo].n_particles);
*/
           printf("%05d/%05d/%05d/%05d",
                  k_halo,
                  halos[k_halo].descendant,
                  halos[k_halo].group_halo_first,
                  halos[k_halo].group_halo_next);
      group_halo_last=halos[k_halo].group_halo_first;
    }
    if(n_halos_snap>0)
      printf(")\n");
    SID_free((void **)&group_halo_first_index);
  }
  SID_free((void **)&snap_num_index);

  SID_free((void **)&snap_num);
  SID_free((void **)&snap_index);
  SID_free((void **)&group_halo_first);
  SID_free((void **)&n_halos);
  SID_free((void **)&halos);
  fclose(fp);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(0);
}
