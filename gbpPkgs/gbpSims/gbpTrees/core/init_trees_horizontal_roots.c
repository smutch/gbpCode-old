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
               i_read_stop,
               i_read_stop-i_read_step,
               n_halos_max,
               MATCH_GROUPS,
               &n_halos_1_matches,
               &n_halos_2_matches,
               n_particles_groups,
               NULL,
               n_subgroups_group[i_file_start%n_wrap],
               NULL,
               match_id,
               match_score,
               match_index,
               NULL,
               F_GOODNESS_OF_MATCH);

  int  n_groups=n_halos_1_matches;

  int i_group;
  for(i_group=0,(*max_id_group)=0,(*max_tree_id_group)=0;i_group<n_groups_max;i_group++){
     for(int i_search=0;i_search<n_wrap;i_search++){
        groups[i_search][i_group].file                             =   i_file_start; // The resulting file offset must be -ve for tree roots
        groups[i_search][i_group].snap                             =   i_read_stop;
        groups[i_search][i_group].index                            =(size_t)i_group;
        groups[i_search][i_group].n_back_matches                   =    0;
        groups[i_search][i_group].descendant.halo                  = NULL;
        groups[i_search][i_group].descendant.score                 =   0.;
        groups[i_search][i_group].descendant.flag_two_way          =FALSE;
        groups[i_search][i_group].descendant.flag_back_match       =FALSE;
        groups[i_search][i_group].first_progenitor.halo            = NULL;
        groups[i_search][i_group].first_progenitor.score           =   0.;
        groups[i_search][i_group].first_progenitor.flag_two_way    =FALSE;
        groups[i_search][i_group].first_progenitor.flag_back_match =FALSE;
        groups[i_search][i_group].last_progenitor.halo             = NULL;
        groups[i_search][i_group].last_progenitor.score            =   0.;
        groups[i_search][i_group].last_progenitor.flag_two_way     =FALSE;
        groups[i_search][i_group].last_progenitor.flag_back_match= FALSE;
        groups[i_search][i_group].next_progenitor.halo             = NULL;
        groups[i_search][i_group].next_progenitor.score            =   0.;
        groups[i_search][i_group].next_progenitor.flag_two_way     =FALSE;
        groups[i_search][i_group].next_progenitor.flag_back_match  =FALSE;
        groups[i_search][i_group].forematch_first.halo             = NULL;
        groups[i_search][i_group].forematch_first.score            =   0.;
        groups[i_search][i_group].forematch_first.flag_two_way     =FALSE;
        groups[i_search][i_group].forematch_first.flag_back_match  =FALSE;
        groups[i_search][i_group].forematch_default.halo           = NULL;
        groups[i_search][i_group].forematch_default.score          =   0.;
        groups[i_search][i_group].forematch_default.flag_two_way   =FALSE;
        groups[i_search][i_group].forematch_default.flag_back_match=FALSE;
        groups[i_search][i_group].forematch_best.halo              = NULL;
        groups[i_search][i_group].forematch_best.score             =   0.;
        groups[i_search][i_group].forematch_best.flag_two_way      =FALSE;
        groups[i_search][i_group].forematch_best.flag_back_match   =FALSE;
        groups[i_search][i_group].bridge_backmatch.halo            = NULL;
        groups[i_search][i_group].bridge_backmatch.score           =   0.;
        groups[i_search][i_group].bridge_backmatch.flag_two_way    =FALSE;
        groups[i_search][i_group].bridge_backmatch.flag_back_match =FALSE;
        groups[i_search][i_group].back_matches                     = NULL;
        groups[i_search][i_group].type                             =TREE_CASE_INVALID;
        groups[i_search][i_group].id                               =-1;
        groups[i_search][i_group].main_progenitor_id               =-1;
        groups[i_search][i_group].tree_id                          =-1;
        groups[i_search][i_group].n_particles                      = 0;
        groups[i_search][i_group].n_particles_parent               = 0;
        groups[i_search][i_group].n_particles_largest_descendant   = 0;
        groups[i_search][i_group].n_progenitors                    = 0;
        if(i_group<n_groups){
           groups[i_search][i_group].id                            =(*max_id_group);
           groups[i_search][i_group].main_progenitor_id            =(*max_id_group);
           groups[i_search][i_group].tree_id                       =(*max_tree_id_group);
           groups[i_search][i_group].type                          =TREE_CASE_MAIN_PROGENITOR|TREE_CASE_MOST_MASSIVE|TREE_CASE_NO_PROGENITORS;
           groups[i_search][i_group].n_particles                   =n_particles_groups[i_group];
           groups[i_search][i_group].n_particles_parent            =n_particles_groups[i_group];
           groups[i_search][i_group].n_particles_largest_descendant=n_particles_groups[i_group];
           groups[i_search][i_group].descendant.halo               =&(groups[(i_search+1)%n_wrap][i_group]);
           groups[i_search][i_group].descendant.score              =1.;
           if(i_search!=(i_file_start%n_wrap)){
              groups[i_search][i_group].n_progenitors=1;
              if(i_search>0){
                 groups[i_search][i_group].first_progenitor.halo =&(groups[i_search-1][i_group]);
                 groups[i_search][i_group].last_progenitor.halo  =&(groups[i_search-1][i_group]);
              }
              else{
                 groups[i_search][i_group].first_progenitor.halo =&(groups[n_wrap-1][i_group]);
                 groups[i_search][i_group].last_progenitor.halo  =&(groups[n_wrap-1][i_group]);
              }
              groups[i_search][i_group].first_progenitor.score=1.;
              groups[i_search][i_group].last_progenitor.score =1.;
           }
        }
     }
     if(i_group<n_groups){
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

  int n_subgroups=n_halos_1_matches;
  int i_subgroup;
  int j_subgroup;
  for(i_group=i_subgroup=j_subgroup=(*max_id_subgroup)=(*max_tree_id_subgroup)=0;i_subgroup<n_subgroups_max;i_subgroup++,j_subgroup++){  
     if(i_group<n_groups){
        while(j_subgroup==n_subgroups_group[i_file_start%n_wrap][i_group] && i_group<(n_groups-1)){
           i_group++;
           j_subgroup=0;
        }
        if(j_subgroup==n_subgroups_group[i_file_start%n_wrap][i_group]){
           i_group++;
           j_subgroup=0;
        }
     }
     for(int i_search=0;i_search<n_wrap;i_search++){
        subgroups[i_search][i_subgroup].file                             =    i_file_start; // The resulting file offset must be -ve for tree roots
        subgroups[i_search][i_subgroup].snap                             =    i_read_stop;
        subgroups[i_search][i_subgroup].index                            = (size_t)i_subgroup;
        subgroups[i_search][i_subgroup].n_back_matches                   =    0;
        subgroups[i_search][i_subgroup].descendant.halo                  = NULL;
        subgroups[i_search][i_subgroup].descendant.score                 =   0.;
        subgroups[i_search][i_subgroup].descendant.flag_two_way          =FALSE;
        subgroups[i_search][i_subgroup].descendant.flag_back_match       =FALSE;
        subgroups[i_search][i_subgroup].first_progenitor.halo            = NULL;
        subgroups[i_search][i_subgroup].first_progenitor.score           =   0.;
        subgroups[i_search][i_subgroup].first_progenitor.flag_two_way    =FALSE;
        subgroups[i_search][i_subgroup].first_progenitor.flag_back_match =FALSE;
        subgroups[i_search][i_subgroup].last_progenitor.halo             = NULL;
        subgroups[i_search][i_subgroup].last_progenitor.score            =   0.;
        subgroups[i_search][i_subgroup].last_progenitor.flag_two_way     =FALSE;
        subgroups[i_search][i_subgroup].last_progenitor.flag_back_match  =FALSE;
        subgroups[i_search][i_subgroup].next_progenitor.halo             = NULL;
        subgroups[i_search][i_subgroup].next_progenitor.score            =   0.;
        subgroups[i_search][i_subgroup].next_progenitor.flag_two_way     =FALSE;
        subgroups[i_search][i_subgroup].next_progenitor.flag_back_match  =FALSE;
        subgroups[i_search][i_subgroup].forematch_first.halo             = NULL;
        subgroups[i_search][i_subgroup].forematch_first.score            =   0.;
        subgroups[i_search][i_subgroup].forematch_first.flag_two_way     =FALSE;
        subgroups[i_search][i_subgroup].forematch_first.flag_back_match  =FALSE;
        subgroups[i_search][i_subgroup].forematch_default.halo           = NULL;
        subgroups[i_search][i_subgroup].forematch_default.score          =   0.;
        subgroups[i_search][i_subgroup].forematch_default.flag_two_way   =FALSE;
        subgroups[i_search][i_subgroup].forematch_default.flag_back_match=FALSE;
        subgroups[i_search][i_subgroup].forematch_best.halo              = NULL;
        subgroups[i_search][i_subgroup].forematch_best.score             =   0.;
        subgroups[i_search][i_subgroup].forematch_best.flag_two_way      =FALSE;
        subgroups[i_search][i_subgroup].forematch_best.flag_back_match   =FALSE;
        subgroups[i_search][i_subgroup].bridge_backmatch.halo            = NULL;
        subgroups[i_search][i_subgroup].bridge_backmatch.score           =   0.;
        subgroups[i_search][i_subgroup].bridge_backmatch.flag_two_way    =FALSE;
        subgroups[i_search][i_subgroup].bridge_backmatch.flag_back_match =FALSE;
        subgroups[i_search][i_subgroup].back_matches           =NULL;
        subgroups[i_search][i_subgroup].type                   =TREE_CASE_INVALID;
        subgroups[i_search][i_subgroup].id                     =-1;
        subgroups[i_search][i_subgroup].main_progenitor_id     =-1;
        subgroups[i_search][i_subgroup].tree_id                =-1;
        subgroups[i_search][i_subgroup].n_particles            = 0;
        subgroups[i_search][i_subgroup].n_particles_parent     = 0;
        subgroups[i_search][i_subgroup].n_particles_largest_descendant= 0;
        subgroups[i_search][i_subgroup].n_progenitors                 = 0;
        if(i_subgroup<n_subgroups){
           subgroups[i_search][i_subgroup].id                            =(*max_id_subgroup);
           subgroups[i_search][i_subgroup].main_progenitor_id            =(*max_id_subgroup);
           subgroups[i_search][i_subgroup].tree_id                       =(*max_tree_id_subgroup);
           subgroups[i_search][i_subgroup].type                          =TREE_CASE_MAIN_PROGENITOR|TREE_CASE_NO_PROGENITORS;
           if(j_subgroup==0)
              subgroups[i_search][i_subgroup].type                      |=TREE_CASE_MOST_MASSIVE;
           subgroups[i_search][i_subgroup].n_particles                   =n_particles_subgroups[i_subgroup];
           subgroups[i_search][i_subgroup].n_particles_parent            =n_particles_groups[i_group];
           subgroups[i_search][i_subgroup].n_particles_largest_descendant=n_particles_subgroups[i_subgroup];
           subgroups[i_search][i_subgroup].descendant.halo               =&(subgroups[(i_search+1)%n_wrap][i_subgroup]);
           subgroups[i_search][i_subgroup].descendant.score              =1.;
           if(i_search!=(i_file_start%n_wrap)){
              subgroups[i_search][i_subgroup].n_progenitors=1;
              if(i_search>0){
                 subgroups[i_search][i_subgroup].first_progenitor.halo=&(subgroups[i_search-1][i_subgroup]);
                 subgroups[i_search][i_subgroup].last_progenitor.halo =&(subgroups[i_search-1][i_subgroup]);
              }
              else{
                 subgroups[i_search][i_subgroup].first_progenitor.halo=&(subgroups[n_wrap-1][i_subgroup]);
                 subgroups[i_search][i_subgroup].last_progenitor.halo =&(subgroups[n_wrap-1][i_subgroup]);
              }
              subgroups[i_search][i_subgroup].first_progenitor.score=1.;
              subgroups[i_search][i_subgroup].last_progenitor.score =1.;
           }
        }
     }
     if(i_subgroup<n_subgroups){
        (*max_id_subgroup)++;
        (*max_tree_id_subgroup)++;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

}

