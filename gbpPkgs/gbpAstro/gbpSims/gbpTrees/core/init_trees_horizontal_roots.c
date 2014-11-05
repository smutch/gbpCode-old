#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpTrees_build.h>

void init_trees_horizontal_roots(tree_horizontal_info **groups,
                                 tree_horizontal_info **subgroups,
                                 int    *match_id,
                                 float  *match_score, 
                                 size_t *match_index,
                                 char   *match_flag_two_way,
                                 int    *n_particles_groups,
                                 int    *n_particles_subgroups, 
                                 int   **n_subgroups_group,
                                 int     n_groups_max, 
                                 int     n_subgroups_max,
                                 char   *filename_root_matches,
                                 int     i_read_start,
                                 int     i_read_stop,
                                 int     i_read_step,
                                 int     i_file_start, 
                                 int     n_wrap,
                                 int     n_halos_max,
                                 int    *max_id_group, 
                                 int    *max_tree_id_group, 
                                 int    *max_id_subgroup, 
                                 int    *max_tree_id_subgroup){
  SID_log("Initializing tree roots...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Initialize everything to a 1:1 simple match
  int     n_halos_1_matches; 
  int     n_halos_2_matches; 
  read_matches(filename_root_matches,
               i_read_stop,i_read_stop-i_read_step,
               n_halos_max,
               MATCH_GROUPS,
               &n_halos_1_matches,
               &n_halos_2_matches,
               n_particles_groups,
               NULL,
               n_subgroups_group[0],
               n_subgroups_group[1],
               match_id,
               match_score,
               match_index,
               NULL,
               F_GOODNESS_OF_MATCH);
  int i_search;
  for(i_search=1;i_search<n_wrap;i_search++)
    memcpy(n_subgroups_group[i_search],n_subgroups_group[0],n_halos_1_matches*sizeof(int));

  int i_halo;
  int j_halo;
  int k_halo;
  for(i_halo=0,(*max_id_group)=0,(*max_tree_id_group)=0;i_halo<n_groups_max;i_halo++,j_halo++){
     for(i_search=0;i_search<n_wrap;i_search++){
        groups[i_search][i_halo].file                   =   i_file_start; // The resulting file offset must be -ve for tree roots
        groups[i_search][i_halo].snap                   =   i_read_stop;
        groups[i_search][i_halo].index                  =(size_t)i_halo;
        groups[i_search][i_halo].n_back_matches         =   0;
        groups[i_search][i_halo].descendant.halo        =NULL;
        groups[i_search][i_halo].descendant.score       =  0.;
        groups[i_search][i_halo].first_progenitor.halo  =NULL;
        groups[i_search][i_halo].first_progenitor.score =  0.;
        groups[i_search][i_halo].last_progenitor.halo   =NULL;
        groups[i_search][i_halo].last_progenitor.score  =  0.;
        groups[i_search][i_halo].next_progenitor.halo   =NULL;
        groups[i_search][i_halo].next_progenitor.score  =  0.;
        groups[i_search][i_halo].forematch_first.halo   =NULL;
        groups[i_search][i_halo].forematch_first.score  =  0.;
        groups[i_search][i_halo].forematch_default.halo =NULL;
        groups[i_search][i_halo].forematch_default.score=  0.;
        groups[i_search][i_halo].bridge_backmatch.halo  =NULL;
        groups[i_search][i_halo].bridge_backmatch.score =  0.;
        groups[i_search][i_halo].back_matches           =NULL;
        groups[i_search][i_halo].type                   =TREE_CASE_INVALID;
        groups[i_search][i_halo].id                     =-1;
        groups[i_search][i_halo].main_progenitor_id     =-1;
        groups[i_search][i_halo].tree_id                =-1;
        groups[i_search][i_halo].n_particles            = 0;
        groups[i_search][i_halo].n_particles_parent     = 0;
        groups[i_search][i_halo].n_particles_largest_descendant= 0;
        groups[i_search][i_halo].n_progenitors                 = 0;
        if(i_halo<n_halos_1_matches){
           groups[i_search][i_halo].id                            =(*max_id_group);
           groups[i_search][i_halo].main_progenitor_id            =(*max_id_group);
           groups[i_search][i_halo].tree_id                       =(*max_tree_id_group);
           groups[i_search][i_halo].type                          =TREE_CASE_MAIN_PROGENITOR|TREE_CASE_NO_PROGENITORS;
           groups[i_search][i_halo].n_particles                   =n_particles_groups[i_halo];
           groups[i_search][i_halo].n_particles_parent            =n_particles_groups[i_halo];
           groups[i_search][i_halo].n_particles_largest_descendant=n_particles_groups[i_halo];
           groups[i_search][i_halo].descendant.halo               =&(subgroups[(i_search+1)%n_wrap][i_halo]);
           groups[i_search][i_halo].descendant.score              =1.;
           if(i_search!=(i_file_start%n_wrap)){
              groups[i_search][i_halo].n_progenitors=1;
              if(i_search>0){
                 groups[i_search][i_halo].first_progenitor.halo =&(groups[i_search-1][i_halo]);
                 groups[i_search][i_halo].last_progenitor.halo  =&(groups[i_search-1][i_halo]);
              }
              else{
                 groups[i_search][i_halo].first_progenitor.halo =&(groups[n_wrap-1][i_halo]);
                 groups[i_search][i_halo].last_progenitor.halo  =&(groups[n_wrap-1][i_halo]);
              }
              groups[i_search][i_halo].first_progenitor.score=1.;
              groups[i_search][i_halo].last_progenitor.score =1.;
           }
        }
     }
     if(i_halo<n_halos_1_matches){
        (*max_id_group)++;
        (*max_tree_id_group)++;
     }
  }

  // Initialize everything to a 1:1 simple match
  read_matches(filename_root_matches,
               i_read_stop,i_read_stop-i_read_step,
               n_halos_max,
               MATCH_SUBGROUPS,
               &n_halos_1_matches,
               &n_halos_2_matches,
               n_particles_subgroups,
               NULL,
               NULL,
               NULL,
               match_id,
               match_score,
               match_index,
               NULL,
               F_GOODNESS_OF_MATCH);

  for(i_halo=0,j_halo=0,k_halo=0,(*max_id_subgroup)=0,(*max_tree_id_subgroup)=0;i_halo<n_subgroups_max;i_halo++,j_halo++){
     if(j_halo>n_subgroups_group[0][k_halo] && i_halo<n_halos_1_matches){
       k_halo++;
       j_halo=0;
     }
     for(i_search=0;i_search<n_wrap;i_search++){
        subgroups[i_search][i_halo].file                   =   i_file_start; // The resulting file offset must be -ve for tree roots
        subgroups[i_search][i_halo].snap                   =   i_read_stop;
        subgroups[i_search][i_halo].index                  =(size_t)i_halo;
        subgroups[i_search][i_halo].n_back_matches         =   0;
        subgroups[i_search][i_halo].descendant.halo        =NULL;
        subgroups[i_search][i_halo].descendant.score       =  0.;
        subgroups[i_search][i_halo].first_progenitor.halo  =NULL;
        subgroups[i_search][i_halo].first_progenitor.score =  0.;
        subgroups[i_search][i_halo].last_progenitor.halo   =NULL;
        subgroups[i_search][i_halo].last_progenitor.score  =  0.;
        subgroups[i_search][i_halo].next_progenitor.halo   =NULL;
        subgroups[i_search][i_halo].next_progenitor.score  =  0.;
        subgroups[i_search][i_halo].forematch_first.halo   =NULL;
        subgroups[i_search][i_halo].forematch_first.score  =  0.;
        subgroups[i_search][i_halo].forematch_default.halo =NULL;
        subgroups[i_search][i_halo].forematch_default.score=  0.;
        subgroups[i_search][i_halo].bridge_backmatch.halo  =NULL;
        subgroups[i_search][i_halo].bridge_backmatch.score =  0.;
        subgroups[i_search][i_halo].back_matches           =NULL;
        subgroups[i_search][i_halo].type                   =TREE_CASE_INVALID;
        subgroups[i_search][i_halo].id                     =-1;
        subgroups[i_search][i_halo].main_progenitor_id     =-1;
        subgroups[i_search][i_halo].tree_id                =-1;
        subgroups[i_search][i_halo].n_particles            = 0;
        subgroups[i_search][i_halo].n_particles_parent     = 0;
        subgroups[i_search][i_halo].n_particles_largest_descendant= 0;
        subgroups[i_search][i_halo].n_progenitors                 = 0;
        if(i_halo<n_halos_1_matches){
           subgroups[i_search][i_halo].id                            =(*max_id_subgroup);
           subgroups[i_search][i_halo].main_progenitor_id            =(*max_id_subgroup);
           subgroups[i_search][i_halo].tree_id                       =(*max_tree_id_subgroup);
           subgroups[i_search][i_halo].type                          =TREE_CASE_MAIN_PROGENITOR|TREE_CASE_NO_PROGENITORS;
           subgroups[i_search][i_halo].n_particles                   =n_particles_subgroups[i_halo];
           subgroups[i_search][i_halo].n_particles_parent            =n_particles_groups[k_halo];
           subgroups[i_search][i_halo].n_particles_largest_descendant=n_particles_subgroups[i_halo];
           subgroups[i_search][i_halo].descendant.halo               =&(subgroups[(i_search+1)%n_wrap][i_halo]);
           subgroups[i_search][i_halo].descendant.score              =1.;
           if(i_search!=(i_file_start%n_wrap)){
              subgroups[i_search][i_halo].n_progenitors=1;
              if(i_search>0){
                 subgroups[i_search][i_halo].first_progenitor.halo=&(subgroups[i_search-1][i_halo]);
                 subgroups[i_search][i_halo].last_progenitor.halo =&(subgroups[i_search-1][i_halo]);
              }
              else{
                 subgroups[i_search][i_halo].first_progenitor.halo=&(subgroups[n_wrap-1][i_halo]);
                 subgroups[i_search][i_halo].last_progenitor.halo =&(subgroups[n_wrap-1][i_halo]);
              }
              subgroups[i_search][i_halo].first_progenitor.score=1.;
              subgroups[i_search][i_halo].last_progenitor.score =1.;
           }
        }
     }
     if(i_halo<n_halos_1_matches){
        (*max_id_subgroup)++;
        (*max_tree_id_subgroup)++;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

}

