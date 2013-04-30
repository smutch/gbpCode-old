#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void compute_trees_horizontal(char        *filename_halo_root_in,
                              char        *filename_cat_root_in,
                              char        *filename_root_matches,
                              char        *filename_output_dir,
                              double      *a_list,
                              cosmo_info **cosmo,
                              int          i_read_start,
                              int          i_read_stop,
                              int          i_read_step,
                              int          n_search,
                              int          flag_fix_bridges,
                              int         *flag_clean){
  char        group_text_prefix[5];
  FILE       *fp;
  char       *line=NULL;
  int         line_length=0;
  int         n_strays;
  int         n_strays_drop;
  int         n_strays_bridge;
  int         i_stray;
  int         n_unlinked;
  int         n_unprocessed;
  int         n_multimatch;
  int         n_sputter;
  int         n_bridge_candidates;
  int         n_bridge_systems;
  int         n_mergers;
  int         n_mergers_drop;
  int         n_mergers_bridge;
  int         n_bridges;
  int         n_match;
  int         n_match_halos;
  int         n_back_match;
  int         i_match;
  int         j_match;
  int         k_match;
  int         n_groups_1;
  int         n_groups_2;
  int         n_groups_3;
  int         i_group;
  int         j_group;
  int         k_group;
  int         l_group;
  int         n_subgroups_1;
  int         n_subgroups_2;
  int         i_subgroup;
  int         j_subgroup;
  int         i_drop;
  int         j_drop;
  int         k_drop;
  int         i_bridge;
  int         j_bridge;
  int         n_lines;
  int         i_file;
  int         j_file;
  int         i_write;
  int         j_write;
  int         l_write;
  int         l_read;
  int         j_file_1;
  int         j_file_2;
  int         i_read;
  int         j_read;
  int         j_read_1;
  int         j_read_2;
  int         n_descendant;
  int         n_progenitor;
  int         descendant_index;
  int         progenitor_index;
  int         my_descendant_index,my_descendant_id,my_descendant_list,my_index;
  int         index;
  int         max_id;
  int         max_id_group;
  int         max_id_subgroup;
  int         my_id;
  int         my_tree;
  int         my_descendant_tree;
  int         my_idx_1;
  int         my_idx_2;
  int        *my_descendant;
  int        *n_particles;
  int        *n_particles_groups;
  int        *n_particles_subgroups;
  int         my_trunk;
  double      expansion_factor;
  double      delta_x,delta_y,delta_z,delta_vx,delta_vy,delta_vz,f_M;
  double      mass_fraction;
  double      max_mass_fraction;
  int         n_drop;
  int         n_found;
  int         n_drop_found;
  int         n_found_bridge;
  double      delta_r;
  double      delta_M;
  double      R_vir_p;
  double      R_vir_d;
  int         i_find,n_find;
  int         flag_continue;
  int         flag_drop;
  int        *match_id=NULL;
  int        *search_id=NULL;
  int         n_progenitors_max;
  int         i_search;
  int         flag_dropped;
  int         flag_first;
  int         n_particles_max;
  int         trunk_index;
  int         biggest_stray;
  int         biggest_stray_drop;
  int         biggest_stray_bridge;
  int         biggest_sputter;
  int        *n_groups=NULL;
  int        *n_subgroups=NULL;
  int         max_tree_id_group;
  int         max_tree_id_subgroup;
  int         max_tree_id;
  int       **n_subgroups_group=NULL;
  int        *n_subgroups_group_1=NULL;
  size_t    **sort_id=NULL;
  size_t    **sort_group_id=NULL;
  size_t    **sort_subgroup_id=NULL;
  size_t     *match_index=NULL;
  size_t     *bridge_index=NULL;
  size_t     *search_index=NULL;
  float      *match_score=NULL;
  int        *bridge_keep=NULL;
  int         flag_match_subgroups;
  int         flag_keep_strays=FALSE;
  int         n_k_match=2;
  int         n_snap;
  
  tree_horizontal_info **subgroups;
  tree_horizontal_info **groups;
  tree_horizontal_info **halos;
  tree_horizontal_info  *halos_i;
  bridge_info           *bridges;
  bridge_info           *bridge;

  int  n_files;
  int  n_subgroups_max;
  int  n_groups_max;
  int *n_halos;
  int  n_halos_max;
  int  n_halos_i;
  int  i_halo;
  int      n_halos_1_matches;
  int      n_halos_2_matches;
  int     j_halo;
  int     k_halo;
  int     l_halo;

  int     n_list;
  int     k_file;
  int     l_file;
  int     k_index;
  int     k_file_temp;
  int     k_index_temp;

  int     n_wrap;
  
  int     n_strayed;
  int     n_sputtered;
  int     n_bridged;
  int     n_dropped;
  int     n_bridged_systems;
  int     max_strayed_size;
  int     max_sputtered_size;
  int     max_dropped_size;
  int     max_emerged_size;
  int     n_bridge_emerged;
  int     i_file_start;
  
  tree_horizontal_stats_info stats;
 
  char  filename_output_dir_horizontal[MAX_FILENAME_LENGTH];
  char  filename_output_dir_horizontal_cases[MAX_FILENAME_LENGTH];
  char  filename_output_file_root[MAX_FILENAME_LENGTH];
  char  filename_matching_out[MAX_FILENAME_LENGTH];
  FILE *fp_matching_out;
  int   i_column;

  SID_log("Constructing horizontal merger trees for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_start,i_read_stop,i_read_step);

  if(n_search<1)
    SID_trap_error("n_search=%d but must be at least 1",ERROR_LOGIC,n_search);

  int flag_compute_fragmented=TRUE;
  int flag_compute_ghosts    =TRUE;

  if(!flag_fix_bridges)
    SID_log("Bridge-fixing is turned off.",SID_LOG_COMMENT);
  if(!flag_compute_fragmented)
    SID_log("Fragmented-halo propagation is turned off.",SID_LOG_COMMENT);
  if(!flag_compute_ghosts)
    SID_log("Ghost-populated tree construction is turned off.",SID_LOG_COMMENT);

  // Validate existing matching files &/or perfrom matching
  compute_trees_matches(filename_halo_root_in,
                        filename_root_matches,
                        i_read_start,
                        i_read_stop,
                        i_read_step,
                        &n_files,
                        &n_subgroups,
                        &n_groups,
                        n_search);

  // We need these for allocating arrays
  calc_max(n_subgroups,&n_subgroups_max,n_files,SID_INT,CALC_MODE_DEFAULT);
  calc_max(n_groups,   &n_groups_max,   n_files,SID_INT,CALC_MODE_DEFAULT);
  n_halos_max=MAX(n_subgroups_max,n_groups_max);

  // We need indices that allow us to hold-on to descendants until outputting
  //   and for the current and last i_file as well
  n_wrap=2*n_search+2;
     
  // Initialize arrays
  SID_log("Creating arrays...",SID_LOG_OPEN);
  n_particles_groups   =(int    *)SID_malloc(sizeof(int)   *n_halos_max);
  n_particles_subgroups=(int    *)SID_malloc(sizeof(int)   *n_halos_max);
  match_id             =(int    *)SID_malloc(sizeof(int)   *n_halos_max);
  match_score          =(float  *)SID_malloc(sizeof(float) *n_halos_max);
  match_index          =(size_t *)SID_malloc(sizeof(size_t)*n_halos_max);
  subgroups            =(tree_horizontal_info **)SID_malloc(sizeof(tree_horizontal_info *)*n_wrap);
  groups               =(tree_horizontal_info **)SID_malloc(sizeof(tree_horizontal_info *)*n_wrap);
  n_subgroups_group    =(int                  **)SID_malloc(sizeof(int                  *)*n_wrap);
  for(i_search=0;i_search<n_wrap;i_search++){
     subgroups[i_search]            =(tree_horizontal_info *)SID_calloc(sizeof(tree_horizontal_info)*n_subgroups_max);
     groups[i_search]               =(tree_horizontal_info *)SID_calloc(sizeof(tree_horizontal_info)*n_groups_max);       
     n_subgroups_group[i_search]    =(int                  *)SID_malloc(sizeof(int)                 *n_groups_max);       
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Process the first file separately
  //   (just give everything ids from a running index) ...
  SID_log("Initializing tree roots...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Initialize everything to a 1:1 simple match
  read_matches(filename_root_matches,
               i_read_stop,i_read_stop-i_read_step,
               MATCH_GROUPS,
               &n_halos_1_matches,
               &n_halos_2_matches,
               n_particles_groups,
               NULL,
               n_subgroups_group[0],
               n_subgroups_group[1],
               match_id,
               match_score,
               match_index);
  for(i_search=1;i_search<n_wrap;i_search++)
    memcpy(n_subgroups_group[i_search],n_subgroups_group[0],n_halos_1_matches*sizeof(int));

  i_file_start=n_files-1;

  for(i_halo=0,max_id_group=0,max_tree_id_group=0;i_halo<n_groups_max;i_halo++,j_halo++){
     for(i_search=0;i_search<n_wrap;i_search++){
        groups[i_search][i_halo].file                  =   i_file_start; // The resulting file offset must be -ve for tree roots
        groups[i_search][i_halo].snap                  =   i_read_start; 
        groups[i_search][i_halo].index                 =(size_t)i_halo;
        groups[i_search][i_halo].n_bridges             =   0;
        groups[i_search][i_halo].descendant.halo       =NULL;
        groups[i_search][i_halo].descendant.score      =  0.;
        groups[i_search][i_halo].first_progenitor.halo =NULL;
        groups[i_search][i_halo].first_progenitor.score=  0.;
        groups[i_search][i_halo].last_progenitor.halo  =NULL;
        groups[i_search][i_halo].last_progenitor.score =  0.;
        groups[i_search][i_halo].next_progenitor.halo  =NULL;
        groups[i_search][i_halo].next_progenitor.score =  0.;
        groups[i_search][i_halo].bridge_forematch.halo =NULL;
        groups[i_search][i_halo].bridge_forematch.score=  0.;
        groups[i_search][i_halo].bridge_backmatch.halo =NULL;
        groups[i_search][i_halo].bridge_backmatch.score=  0.;
        groups[i_search][i_halo].bridges               =NULL;
        groups[i_search][i_halo].type                  =TREE_CASE_INVALID;
        groups[i_search][i_halo].id                    =-1;
        groups[i_search][i_halo].main_progenitor_id    =-1;
        groups[i_search][i_halo].tree_id               =-1; 
        groups[i_search][i_halo].n_particles           = 0;
        groups[i_search][i_halo].n_particles_parent    = 0;
        groups[i_search][i_halo].n_progenitors         = 0;
        if(i_halo<n_halos_1_matches){
           groups[i_search][i_halo].id                    =max_id_group;
           groups[i_search][i_halo].main_progenitor_id    =max_id_group;
           groups[i_search][i_halo].tree_id               =max_tree_id_group;
           groups[i_search][i_halo].type                  =TREE_CASE_SIMPLE|TREE_CASE_MAIN_PROGENITOR|TREE_CASE_NO_PROGENITORS;
           groups[i_search][i_halo].n_particles           =n_particles_groups[i_halo];
           groups[i_search][i_halo].n_particles_parent    =n_particles_groups[i_halo];
           groups[i_search][i_halo].descendant.halo       =&(subgroups[(i_search+1)%n_wrap][i_halo]);
           groups[i_search][i_halo].descendant.score       =1.;
           if(i_search!=(i_file_start%n_wrap)){
              groups[i_search][i_halo].n_progenitors       =1;
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
        max_id_group++;
        max_tree_id_group++;
     }
  }

  // Initialize everything to a 1:1 simple match
  read_matches(filename_root_matches,
               i_read_stop,i_read_stop-i_read_step,
               MATCH_SUBGROUPS,
               &n_halos_1_matches,
               &n_halos_2_matches,
               n_particles_subgroups,
               NULL,
               NULL,
               NULL,
               match_id,
               match_score,
               match_index);
                     
  for(i_halo=0,j_halo=0,k_halo=0,max_id_subgroup=0,max_tree_id_subgroup=0;i_halo<n_subgroups_max;i_halo++,j_halo++){
     if(j_halo>n_subgroups_group[0][k_halo] && i_halo<n_halos_1_matches){
       k_halo++;
       j_halo=0;
     }
     for(i_search=0;i_search<n_wrap;i_search++){
        subgroups[i_search][i_halo].file                  =   i_file_start; // The resulting file offset must be -ve for tree roots
        subgroups[i_search][i_halo].snap                  =   i_read_start; 
        subgroups[i_search][i_halo].index                 =(size_t)i_halo;
        subgroups[i_search][i_halo].n_bridges             =   0;
        subgroups[i_search][i_halo].descendant.halo       =NULL;
        subgroups[i_search][i_halo].descendant.score      =  0.;
        subgroups[i_search][i_halo].first_progenitor.halo =NULL;
        subgroups[i_search][i_halo].first_progenitor.score=  0.;
        subgroups[i_search][i_halo].last_progenitor.halo  =NULL;
        subgroups[i_search][i_halo].last_progenitor.score =  0.;
        subgroups[i_search][i_halo].next_progenitor.halo  =NULL;
        subgroups[i_search][i_halo].next_progenitor.score =  0.;
        subgroups[i_search][i_halo].bridge_forematch.halo =NULL;
        subgroups[i_search][i_halo].bridge_forematch.score=  0.;
        subgroups[i_search][i_halo].bridge_backmatch.halo =NULL;
        subgroups[i_search][i_halo].bridge_backmatch.score=  0.;
        subgroups[i_search][i_halo].bridges               =NULL;
        subgroups[i_search][i_halo].type                  =TREE_CASE_INVALID;
        subgroups[i_search][i_halo].id                    =-1;
        subgroups[i_search][i_halo].main_progenitor_id    =-1;
        subgroups[i_search][i_halo].tree_id               =-1; 
        subgroups[i_search][i_halo].n_particles           = 0;
        subgroups[i_search][i_halo].n_particles_parent    = 0;
        subgroups[i_search][i_halo].n_progenitors         = 0;
        if(i_halo<n_halos_1_matches){
           subgroups[i_search][i_halo].id                    =max_id_subgroup;
           subgroups[i_search][i_halo].main_progenitor_id    =max_id_subgroup;
           subgroups[i_search][i_halo].tree_id               =max_tree_id_subgroup;
           subgroups[i_search][i_halo].type                  =TREE_CASE_SIMPLE|TREE_CASE_MAIN_PROGENITOR|TREE_CASE_NO_PROGENITORS;
           subgroups[i_search][i_halo].n_particles           =n_particles_subgroups[i_halo];
           subgroups[i_search][i_halo].n_particles_parent    =n_particles_groups[k_halo];
           subgroups[i_search][i_halo].descendant.halo       =&(subgroups[(i_search+1)%n_wrap][i_halo]);
           subgroups[i_search][i_halo].descendant.score       =1.;
           if(i_search!=(i_file_start%n_wrap)){
              subgroups[i_search][i_halo].n_progenitors       =1;
              if(i_search>0){
                 subgroups[i_search][i_halo].first_progenitor.halo =&(subgroups[i_search-1][i_halo]);
                 subgroups[i_search][i_halo].last_progenitor.halo  =&(subgroups[i_search-1][i_halo]);
              }
              else{
                 subgroups[i_search][i_halo].first_progenitor.halo =&(subgroups[n_wrap-1][i_halo]);
                 subgroups[i_search][i_halo].last_progenitor.halo  =&(subgroups[n_wrap-1][i_halo]);
              }
              subgroups[i_search][i_halo].first_progenitor.score=1.;
              subgroups[i_search][i_halo].last_progenitor.score =1.;
           }
        }
     }
     if(i_halo<n_halos_1_matches){
        max_id_subgroup++;
        max_tree_id_subgroup++;
     }
  }

  SID_log("Done.",SID_LOG_CLOSE);

  // The first snapshot is done now (set to defaults as the roots of trees) ... now loop over all other snapshots ...
  //   There are a bunch of counters at work here.  Because we aren't necessarily using every 
  //     snapshot (if i_read_step>1), we need counters to keep track of which snapshots we
  //     are working with (i_read_*,j_read_*, etc), counters to keep track of which
  //     files's we're dealing with as far as the trees indices are concerned (i_file_*,j_file_*,etc), and
  //     counters to keep track of which files are being/have been written (i_write_*,j_write_* etc).
  //     We can't write files right away because previously processed snapshots can be changed
  //     when we deal with dropped and bridged halos.
  for(i_read   =i_read_stop-i_read_step,
        i_file =i_file_start-1, 
        j_file =1,             
        i_write=i_file_start,      
        j_write=i_read_stop,
        l_write=0;      
      i_read>=i_read_start;
      i_read-=i_read_step,    
         i_file--, 
         j_file++){
    SID_log("Processing snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read);

    // Loop twice (1st to process subgroups, 2nd to process groups)
    for(k_match=0;k_match<n_k_match;k_match++){

       // Initialize a bunch of stuff depending on whether
       //   we are processing groups or subgroups
       switch(k_match){
          case 0:
          sprintf(group_text_prefix,"sub");
          flag_match_subgroups=MATCH_SUBGROUPS;
          halos               =subgroups;
          n_halos             =n_subgroups;
          n_halos_max         =n_subgroups_max;
          max_id              =max_id_subgroup;
          n_particles         =n_particles_subgroups;
          break;
          case 1:
          sprintf(group_text_prefix,"");
          flag_match_subgroups=MATCH_GROUPS;
          halos               =groups;
          n_halos             =n_groups;
          n_halos_max         =n_groups_max;
          max_id              =max_id_group;
          n_particles         =n_particles_groups;
          break;
       }
       halos_i  =halos[i_file%n_wrap];
       n_halos_i=n_halos[j_file];
       SID_log("Processing %d %sgroups...",SID_LOG_OPEN|SID_LOG_TIMER,n_halos_i,group_text_prefix);

       // Initialize tree pointer-arrays with dummy values
       for(i_halo=0;i_halo<n_halos_max;i_halo++){
          halos_i[i_halo].file                  = i_file;
          halos_i[i_halo].snap                  = i_read;
          halos_i[i_halo].index                 = (size_t)i_halo;
          halos_i[i_halo].n_bridges             =   0;
          halos_i[i_halo].descendant.halo       =NULL;
          halos_i[i_halo].descendant.score      =  0.;
          halos_i[i_halo].first_progenitor.halo =NULL;
          halos_i[i_halo].first_progenitor.score=  0.;
          halos_i[i_halo].last_progenitor.halo  =NULL;
          halos_i[i_halo].last_progenitor.score =  0.;
          halos_i[i_halo].next_progenitor.halo  =NULL;
          halos_i[i_halo].next_progenitor.score =  0.;
          halos_i[i_halo].bridge_forematch.halo =NULL;
          halos_i[i_halo].bridge_forematch.score=  0.;
          halos_i[i_halo].bridge_backmatch.halo =NULL;
          halos_i[i_halo].bridge_backmatch.score=  0.;
          SID_free(SID_FARG halos_i[i_halo].bridges);
          if(i_halo<n_halos_i)
             halos_i[i_halo].type=TREE_CASE_UNPROCESSED|TREE_CASE_NO_PROGENITORS;
          else
             halos_i[i_halo].type=TREE_CASE_INVALID;
          halos_i[i_halo].id                =-1;
          halos_i[i_halo].main_progenitor_id=-1;
          halos_i[i_halo].tree_id           =-1;
          halos_i[i_halo].n_particles       = 0;
          halos_i[i_halo].n_progenitors     = 0;
       }
       
       // Use back-matching to identify bridged halos ...
       if(flag_fix_bridges){
          SID_log("Identifying bridge candidates from back-matching...",SID_LOG_OPEN|SID_LOG_TIMER);
          SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
          //    ... first, do an initial count of matches.  This will not be a list of unique halos
          //        though, since the same halos are likely to appear in repeated snapshots.
          for(j_file_1  =i_file+1,
                j_file_2=i_file,
                j_read_1=i_read+i_read_step,
                j_read_2=i_read,
                i_search=0;
              j_read_1<=i_read_stop && i_search<n_search;
              j_file_1++,
                j_read_1+=i_read_step,
                i_search++){

             SID_log("Counting matches between files %d->%d...",SID_LOG_OPEN,j_read_1,j_read_2);

             // Read back-matching
             read_matches(filename_root_matches,
                          j_read_1,j_read_2,
                          flag_match_subgroups,
                          &n_halos_1_matches,
                          &n_halos_2_matches,
                          NULL,
                          n_particles,
                          NULL,
                          NULL,
                          match_id,
                          match_score,
                          match_index);

             // Store halo sizes
             if(i_search==0){
                for(i_halo=0;i_halo<n_halos_2_matches;i_halo++)
                   halos[i_file%n_wrap][i_halo].n_particles=n_particles[i_halo];
             }

             // Perform initial back-match count
             for(i_halo=0;i_halo<n_halos_i;i_halo++){
                j_halo=find_index_int(match_id,i_halo,n_halos_1_matches,match_index);
                while(match_id[match_index[j_halo]]==i_halo && j_halo<(n_halos_1_matches-1)){
                   halos_i[i_halo].n_bridges++;
                   j_halo++;
                }
                if(match_id[match_index[j_halo]]==i_halo && j_halo==(n_halos_1_matches-1)){
                   halos_i[i_halo].n_bridges++;
                   j_halo++;
                }
             }
             SID_log("Done.",SID_LOG_CLOSE);
          }
       
          //    ... second, do a conservative allocation using the non-unique counts and reset the counter.
          for(i_halo=0;i_halo<n_halos_i;i_halo++){
             if((halos_i[i_halo].n_bridges)>0)
                (halos_i[i_halo].bridges)=(bridge_info *)SID_calloc(sizeof(bridge_info)*(halos_i[i_halo].n_bridges));
             else
                (halos_i[i_halo].bridges)=NULL;
             halos_i[i_halo].n_bridges =0;
          }

          //    ... third, assemble the list of unique back-matched halos.
          for(j_file_1  =i_file+1,
                j_file_2=i_file,
                j_read_1=i_read+i_read_step,
                j_read_2=i_read,
                i_search=0;
              j_read_1<=i_read_stop && i_search<n_search;
              j_file_1++,
                j_read_1+=i_read_step,
                i_search++){

             SID_log("Finding unique matches between files %d->%d...",SID_LOG_OPEN,j_read_1,j_read_2);
          
             // Read back-matching
             read_matches(filename_root_matches,
                          j_read_1,j_read_2,
                          flag_match_subgroups,
                          &n_halos_1_matches,
                          &n_halos_2_matches,
                          NULL,
                          NULL,
                          NULL,
                          NULL,
                          match_id,
                          match_score,
                          match_index);
                             
             // For all the halos in i_file_1 with back-matches ...
             for(i_halo=0;i_halo<n_halos_i;i_halo++){
                if((halos_i[i_halo].bridges)!=NULL){
                   // Scan over the list of halos from snapshot=j_read_1 
                   //   that match this halo in j_read_2 ...
                   bridges=halos_i[i_halo].bridges;
                   j_halo =find_index_int(match_id,i_halo,n_halos_1_matches,match_index);
                   // Loop over all but the last halo in the list ...
                   while(match_id[match_index[j_halo]]==i_halo && j_halo<(n_halos_1_matches-1)){
                      // Check to see if this halo is already in the bridge list (keep all sputtered/strayed halos as well) ...
                      flag_continue=TRUE;
                      if(halos[j_file_1%n_wrap][match_index[j_halo]].id>=0){
                         for(k_halo=0,flag_continue=TRUE;k_halo<halos_i[i_halo].n_bridges && flag_continue;k_halo++){
                            if(bridges[k_halo].halo->id==halos[j_file_1%n_wrap][match_index[j_halo]].id)
                               flag_continue=FALSE;
                         }
                      }
                      // ... if not, add it
                      if(flag_continue){
                         bridges[halos_i[i_halo].n_bridges].score=match_score[match_index[j_halo]];
                         bridges[halos_i[i_halo].n_bridges].halo =&(halos[j_file_1%n_wrap][match_index[j_halo]]);
                         (halos_i[i_halo].n_bridges)++;
                      }
                      j_halo++;
                   }
                   // ... then do the last halo in the list ...
                   if(match_id[match_index[j_halo]]==i_halo && j_halo==(n_halos_1_matches-1)){
                      // Check to see if this halo is already in the list (keep all sputtered/strayed halos as well) ...
                      flag_continue=TRUE;
                      if(halos[j_file_1%n_wrap][match_index[j_halo]].id>=0){
                         for(k_halo=0;k_halo<halos_i[i_halo].n_bridges && flag_continue;k_halo++){
                            if(bridges[k_halo].halo->id==halos[j_file_1%n_wrap][match_index[j_halo]].id)
                               flag_continue=FALSE;
                         }
                      }
                      // ... if not, add it
                      if(flag_continue){
                         bridges[halos_i[i_halo].n_bridges].score=match_score[match_index[j_halo]];
                         bridges[halos_i[i_halo].n_bridges].halo =&(halos[j_file_1%n_wrap][match_index[j_halo]]);
                         (halos_i[i_halo].n_bridges)++;
                      }
                      j_halo++;
                   }
                }
             }
             SID_log("Done.",SID_LOG_CLOSE);
          }

          // ... lastly, reorder the emerged halos by score, keep only the most immediate bridge descendants and finalize the list ...
          SID_log("Re-ordering bridges...",SID_LOG_OPEN);
          for(i_halo=0;i_halo<n_halos_i;i_halo++){
             if((halos_i[i_halo].n_bridges)>1){ 

                // We may need to remove several halos from the list.  This array will keep track of this.
                bridge_keep=(int *)SID_malloc(sizeof(int)*halos_i[i_halo].n_bridges);

                // Reorder the bridges by their score.  We make a temporary copy of the list 
                //   to do this and initially set all bridges as halos to keep..
                bridges=(bridge_info *)SID_calloc(sizeof(bridge_info)*(halos_i[i_halo].n_bridges));
                for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;j_halo++){
                   bridge=&(halos_i[i_halo].bridges[j_halo]);
                   memcpy(&(bridges[j_halo]),bridge,sizeof(bridge_info));
                   match_score[j_halo]=bridge->score;
                   bridge_keep[j_halo]=TRUE;
                }
                merge_sort((void *)match_score,(size_t)(halos_i[i_halo].n_bridges),&bridge_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

                // Remove any mutual descendants from the list
                //   (since they have their own IDs, this is 
                //    needed to avoid calling them emerged halos)
                for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;j_halo++){
                   bridge = &(bridges[bridge_index[j_halo]]);
                   tree_horizontal_info *current;
                   // ... walk the tree upwards ...
                   current=bridge->halo->descendant.halo;
                   if(current!=NULL)
                      k_file=current->file;
                   l_file=k_file;
                   while(current!=NULL && k_file>=l_file && k_file<MIN(n_files,i_file+(n_search+1))){
                      for(k_halo=0;k_halo<halos_i[i_halo].n_bridges;k_halo++){
                         bridge = &(bridges[bridge_index[k_halo]]);
                         if(bridge->halo==current)
                            bridge_keep[k_halo]=FALSE;
                      }
                      current=current->descendant.halo;
                      l_file=k_file;
                      if(current!=NULL)
                         k_file =current->file;
                   }
                }

                // Since we may have trimmed the list, recount the number remaining
                n_list=halos_i[i_halo].n_bridges;
                for(j_halo=n_list-1,halos_i[i_halo].n_bridges=0;j_halo>=0;j_halo--){
                  if(bridge_keep[j_halo])
                     halos_i[i_halo].n_bridges++;
                }

                // We've removed some halos and may not actually be a bridged halo anymore.  Clean-up if so.
                if(halos_i[i_halo].n_bridges<1){
                   halos_i[i_halo].type&=(~TREE_CASE_BRIDGED);
                   SID_free(SID_FARG halos_i[i_halo].bridges);
                   halos_i[i_halo].n_bridges=0;
                }
                else{
                   halos_i[i_halo].type|=TREE_CASE_BRIDGED;
                   n_bridged++;

                   // Because we've overallocated previously, reallocate the bridge list here to save RAM.
                   SID_free(SID_FARG halos_i[i_halo].bridges);
                   (halos_i[i_halo].bridges)=(bridge_info *)SID_calloc(sizeof(bridge_info)*(halos_i[i_halo].n_bridges));

                   // Copy the sorted temporary list to the permanent list.
                   for(j_halo=n_list-1,l_halo=0;j_halo>=0;j_halo--){
                      if(bridge_keep[j_halo]){
                         memcpy(&(halos_i[i_halo].bridges[l_halo]),&(bridges[bridge_index[j_halo]]),sizeof(bridge_info));
                         if(halos[(bridges[bridge_index[j_halo]].halo->file)%n_wrap][bridges[bridge_index[j_halo]].halo->index].bridge_backmatch.halo==NULL){
                            halos[(bridges[bridge_index[j_halo]].halo->file)%n_wrap][bridges[bridge_index[j_halo]].halo->index].bridge_backmatch.halo =&(halos_i[i_halo]);
                            halos[(bridges[bridge_index[j_halo]].halo->file)%n_wrap][bridges[bridge_index[j_halo]].halo->index].bridge_backmatch.score=match_score[bridge_index[j_halo]];
                         }
                         l_halo++;
                      }
                   }
                }

                // Clean-up
                SID_free(SID_FARG bridge_keep);
                SID_free(SID_FARG bridge_index);
                SID_free(SID_FARG bridges);
             }
             // This halo is not a bridge.  Perform cleaning.
             else{
                halos_i[i_halo].type&=(~TREE_CASE_BRIDGED);
                SID_free(SID_FARG halos_i[i_halo].bridges);
                halos_i[i_halo].n_bridges=0;
             }
          }
          SID_log("Done.",SID_LOG_CLOSE);
          SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
          SID_log("Done.",SID_LOG_CLOSE);
       }

       // Perform forward-matching
       SID_log("Constructing progenitors from forward-matching...",SID_LOG_OPEN|SID_LOG_TIMER);
       for(j_file_1  =i_file,
             j_file_2=i_file+1,
             j_read_1=i_read,
             j_read_2=i_read+i_read_step,
             i_search=0;
           j_read_2<=i_read_stop && i_search<n_search;
           j_file_2++,
             j_read_2+=i_read_step,
             i_search++){

          // Read forward-matching
          read_matches(filename_root_matches,
                       j_read_1,j_read_2,
                       flag_match_subgroups,
                       &n_halos_1_matches,
                       &n_halos_2_matches,
                       n_particles,
                       NULL,
                       n_subgroups_group[j_file_1%n_wrap],
                       n_subgroups_group[j_file_2%n_wrap],
                       match_id,
                       match_score,
                       match_index);

          // Store halo sizes
          if(!flag_fix_bridges){
             if(i_search==0){
                for(i_halo=0;i_halo<n_halos_1_matches;i_halo++)
                   halos[i_file%n_wrap][i_halo].n_particles=n_particles[i_halo];
             }
          }
 
          // Perform matching for all the halos in i_file_1.  This loop should deal completely with
          //   all simple matches and dropped halos.  It also identifies matches to bridges, which
          //   require special treatment in the loop that follows (to look for matches to emergent halos)
          //   and at the end of the loop over j_read/i_search to finalize those not matched to emerged halos.
          for(i_halo=0;i_halo<n_halos_1_matches;i_halo++){
             // If this halo hasn't been processed during earlier searches ...
             if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_UNPROCESSED)){
                my_descendant_index=match_id[i_halo];
                // If this halo has been matched to something in i_file_2 ...
                if(my_descendant_index>=0){
                   if(my_descendant_index<n_halos_2_matches)
                      set_halo_and_descendant(halos,
                                              i_file,
                                              i_halo,
                                              j_file_2,
                                              my_descendant_index,
                                              match_score[i_halo],
                                              &max_id,
                                              n_wrap);
                   else
                      SID_log_warning("descendant ID out of bounds (ie. %d>%d) in snapshot %03d -> snapshot %03d %sgroup matching for i_halo=%d.",
                                      SID_WARNING_DEFAULT,my_descendant_index,n_halos_2_matches-1,j_read_1,j_read_2,group_text_prefix,i_halo);
                }
             }
          }

          // Try to match halos-matched-to-bridges to the candidate emergent halos identified in the bridge-lists
          for(i_halo=0,n_drop=0;i_halo<n_halos_1_matches;i_halo++){
             if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED)){
                // Loop over all the emergent halos identified with the bridge that i_halo has been matched to
                if(halos_i[i_halo].bridge_forematch.halo==NULL)
                  SID_trap_error("Bridge match not defined during emerged halo search.",ERROR_LOGIC);
                bridges=halos[(halos_i[i_halo].bridge_forematch.halo->file)%n_wrap][halos_i[i_halo].bridge_forematch.halo->index].bridges;
                if(bridges==NULL)
                  SID_trap_error("Bridges not defined during emerged halo search.",ERROR_LOGIC);
                n_list=halos[(halos_i[i_halo].bridge_forematch.halo->file)%n_wrap][halos_i[i_halo].bridge_forematch.halo->index].n_bridges;
                for(k_halo=0;k_halo<n_list;k_halo++){
                   if(bridges[k_halo].halo->file==j_file_2 && match_id[i_halo]==bridges[k_halo].halo->index){
                      set_halo_and_descendant(halos,
                                              i_file,
                                              i_halo,
                                              bridges[k_halo].halo->file,
                                              bridges[k_halo].halo->index,
                                              match_score[i_halo],
                                              &max_id,
                                              n_wrap);
                   }
                }
             }
          }
       }
       SID_log("Done.",SID_LOG_CLOSE);
      
       // Finalize matches to unprocessed halos ...
       SID_log("Applying defaults to unprocessed halos...",SID_LOG_OPEN|SID_LOG_TIMER);
       // ... first unprocessed matches to bridged halos (apply default behavior)
       for(i_halo=0,n_drop=0;i_halo<n_halos_1_matches;i_halo++){
          if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED)){
             // The descendant info is already set.  Just increment it's progenitor counter and set this halo's info.
             //   Leave the target flag untouched so we can later identify which BRIDGE_PROGENITORS were found
             halos_i[i_halo].type|=(TREE_CASE_BRIDGE_FINALIZE|TREE_CASE_BRIDGE_DEFAULT);
             set_halo_and_descendant(halos,
                                     i_file,
                                     i_halo,
                                     halos_i[i_halo].bridge_forematch.halo->file,
                                     halos_i[i_halo].bridge_forematch.halo->index,
                                     halos_i[i_halo].bridge_forematch.score,
                                     &max_id,
                                     n_wrap);
          }
       }
       // ... then assign flags for halos not successfully processed.  They must be strays.
       for(i_halo=0;i_halo<n_halos_i;i_halo++){
          if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_UNPROCESSED)){
             halos_i[i_halo].type|=TREE_CASE_STRAYED;
             halos_i[i_halo].type&=(~TREE_CASE_UNPROCESSED);
          }
       }
       SID_log("Done.",SID_LOG_CLOSE);
 
       // Now that we have assigned all the IDs for the halos in the active snapshot,
       //   we need to remove all descendants of the bridge from the lists of candidate emerged halos 
       //   to avoid incorrectly matching to the main progenitor's descendants later-on when when 
       //   we are scaning emerged halo candidates.  Real matches to bridges are dealt-with
       //   when halos marked TREE_CASE_BRIDGE_DEFAULT are processed. 
       SID_log("Removing main progenitors from candidate emerged halo lists...",SID_LOG_OPEN|SID_LOG_TIMER);
       for(i_halo=0;i_halo<n_halos_i;i_halo++){
          // Check all bridged halos ...
          if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGED)){
             n_list=halos_i[i_halo].n_bridges;
             // ... and check all of their descendants ...
             tree_horizontal_info *current;
             current=halos_i[i_halo].descendant.halo;
             if(current!=NULL)
                k_file=current->file;
             l_file=k_file;
             while(current!=NULL && k_file>=l_file && k_file<MIN(n_files,i_file+(n_search+1))){
                for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;){
                   bridge=&(halos_i[i_halo].bridges[j_halo]);
                   if(bridge->halo==current){
                      // Remove halos by decrementing the counter ...
                      halos_i[i_halo].n_bridges--;
                      // ... and sliding all the halos down ...
                      for(k_halo=j_halo;k_halo<halos_i[i_halo].n_bridges;k_halo++)
                         memcpy(&(halos_i[i_halo].bridges[k_halo]),&(halos_i[i_halo].bridges[k_halo+1]),sizeof(bridge_info));
                   }
                   else
                     j_halo++; // We only need to increment the counter if we don't find a match
                }
                current=current->descendant.halo;
                l_file=k_file;
                if(current!=NULL)
                   k_file =current->file;
             }

             // Since we may have removed items, we might not have a bridged halo any more.
             if(halos_i[i_halo].n_bridges<1){
                halos_i[i_halo].type&=(~TREE_CASE_BRIDGED);
                SID_free(SID_FARG halos_i[i_halo].bridges);
                halos_i[i_halo].n_bridges=0;
             }
             // If this halo is still bridged, label the halos in the remaining back-match list as candidate emerged halos
             else{
                for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;j_halo++){
                   bridge = &(halos_i[i_halo].bridges[j_halo]);
                   bridge->halo->type|=TREE_CASE_EMERGED_CANDIDATE;
                }
             }
          }
       }
       SID_log("Done.",SID_LOG_CLOSE);

       // This has to be written right after a snapshot is read and processed (because it needs all forward scan information), 
       //    so it is separate from the rest of the log output code.
       SID_log("Writing candidate emerged halo information...",SID_LOG_OPEN|SID_LOG_TIMER);
       sprintf(filename_output_dir_horizontal,      "%s/horizontal",   filename_output_dir);
       sprintf(filename_output_dir_horizontal_cases,"%s/special_cases",filename_output_dir_horizontal);
       mkdir(filename_output_dir,                   02755);
       mkdir(filename_output_dir_horizontal,        02755);
       mkdir(filename_output_dir_horizontal_cases,  02755);
       strcpy(filename_output_file_root,filename_output_dir);
       strip_path(filename_output_file_root);
       sprintf(filename_matching_out,"%s/%s.%sgroups_emerged_candidates_%d",filename_output_dir_horizontal_cases,filename_output_file_root,group_text_prefix,i_read);
       fp_matching_out=fopen(filename_matching_out,"w");
       i_column=1;
       fprintf(fp_matching_out,"# (%02d): Bridge number\n",                    i_column++);
       fprintf(fp_matching_out,"# (%02d): Halo file\n",                        i_column++);
       fprintf(fp_matching_out,"# (%02d): Halo index\n",                       i_column++);
       fprintf(fp_matching_out,"# (%02d): Halo ID\n",                          i_column++);
       fprintf(fp_matching_out,"# (%02d): Tree ID\n",                          i_column++);
       fprintf(fp_matching_out,"# (%02d): Descendant file\n",                  i_column++);
       fprintf(fp_matching_out,"# (%02d): Descendant index\n",                 i_column++);
       fprintf(fp_matching_out,"# (%02d): Descendant ID\n",                    i_column++);
       fprintf(fp_matching_out,"# (%02d): Bridge file\n",                      i_column++);
       fprintf(fp_matching_out,"# (%02d): Bridge index\n",                     i_column++);
       fprintf(fp_matching_out,"# (%02d): Bridge ID\n",                        i_column++);
       fprintf(fp_matching_out,"# (%02d): Descendant match score\n",           i_column++);
       fprintf(fp_matching_out,"# (%02d): Bridge     match score\n",           i_column++);
       fprintf(fp_matching_out,"# (%02d): Number of particles in halo\n",      i_column++);
       fprintf(fp_matching_out,"# (%02d): Number of particles in descendant\n",i_column++);
       fprintf(fp_matching_out,"# (%02d): Number of particles in bridge\n",    i_column++);
       for(i_halo=0;i_halo<n_halos_i;i_halo++){
          // Print bridge info
          if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGED)){
             for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;j_halo++){
                bridge=&(halos_i[i_halo].bridges[j_halo]);
                fprintf(fp_matching_out,"%3d %7d   %4d %7d %7d   %4d %7d %7d   %4d %7d %7d   %10.4f %10.4f   %7d %7d %7d\n",
                        j_halo,
                        i_read,
                        halos_i[i_halo].tree_id,
                        i_halo,
                        halos_i[i_halo].id,
                        set_match_snapshot(&(halos_i[i_halo].descendant)),
                        set_match_index(&(halos_i[i_halo].descendant)),
                        set_match_id(&(halos_i[i_halo].descendant)),
                        set_match_snapshot(bridge),
                        set_match_index(bridge),
                        set_match_id(bridge),
                        set_match_score(&(halos_i[i_halo].descendant)),
                        set_match_score(bridge),
                        halos_i[i_halo].n_particles,
                        set_match_n_particles(&(halos_i[i_halo].descendant)),
                        set_match_n_particles(bridge));
             }
          }
       }
       fclose(fp_matching_out);
       SID_log("Done.",SID_LOG_CLOSE);

       // Report some statistics
       //   n.b.: This is only an estimate in some cases, since subsequent snapshots may alter this snapshot.  
       //         See the written log file for the most accurate numbers.
       compute_trees_horizontal_stats(halos_i,n_halos_i,n_halos_max,&stats,TRUE);
       SID_log("Results (estimates which may change with continued processing):",SID_LOG_OPEN);
       SID_log("# of halos                  =%-8d",SID_LOG_COMMENT,stats.n_halos);
       SID_log("# of simple matches         =%-8d (%d mergers)",SID_LOG_COMMENT,stats.n_simple,stats.n_mergers);
       if(stats.n_strayed>0)
          SID_log("# of strayed halos          =%-8d (largest=%d particles)",SID_LOG_COMMENT,stats.n_strayed,stats.max_strayed_size);
       else
          SID_log("# of strayed halos          =%-8d",SID_LOG_COMMENT,stats.n_strayed);
       if(stats.n_sputtered>0)
          SID_log("# of sputtering halos       =%-8d (largest=%d particles)",SID_LOG_COMMENT,stats.n_sputtered,stats.max_sputtered_size);
       else
          SID_log("# of sputtering halos       =%-8d",SID_LOG_COMMENT,stats.n_sputtered);
       if(stats.n_dropped>0)
          SID_log("# of dropped halos          =%-8d (largest=%d particles)",SID_LOG_COMMENT,stats.n_dropped,stats.max_dropped_size);
       else
          SID_log("# of dropped halos          =%-8d",SID_LOG_COMMENT,stats.n_dropped);
       SID_log("# of bridged halos          =%-8d",SID_LOG_COMMENT,stats.n_bridged);
       SID_log("# of bridge progenitors     =%-8d",SID_LOG_COMMENT,stats.n_bridge_progenitors);
       if(stats.n_emerged_progenitors>0)
          SID_log("# of emerged progenitors    =%-8d (largest=%d particles)",SID_LOG_COMMENT,stats.n_emerged_progenitors,stats.max_emerged_progenitor_size);
       else
          SID_log("# of emerged progenitors    =%-8d",SID_LOG_COMMENT,stats.n_emerged_progenitors);      
       SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);
      
       // Update things changed in this k_match iteration
       switch(k_match){
          case 0:
            max_id_subgroup=max_id;
            break;
          case 1:
            max_id_group=max_id;
            break;
       }
       SID_log("Done.",SID_LOG_CLOSE);
    } // k_match
    
    // Write tree info once a few files have been processed
    //   and no more dropped groups need to be given ids
    if(j_file>n_search){
       check_for_fragmented_halos(0,subgroups,n_subgroups[l_write],i_write,j_write,n_wrap);
       check_for_fragmented_halos(1,groups,   n_groups[l_write],   i_write,j_write,n_wrap);
       write_trees_horizontal((void **)groups, 
                              (void **)subgroups,
                              n_groups[l_write],   n_groups_max,   0,
                              n_subgroups[l_write],n_subgroups_max,0,
                              n_subgroups_group,
                              max_tree_id_subgroup,
                              max_tree_id_group,
                              i_write,
                              j_write,
                              l_write,
                              n_wrap,
                              i_file_start,
                              filename_cat_root_in,
                              filename_output_dir,
                              a_list,
                              cosmo,
                              n_k_match,
                              TREE_HORIZONTAL_WRITE_EXTENDED|TREE_HORIZONTAL_WRITE_ALLCASES);
       i_write--;
       l_write++;
       j_write-=i_read_step;
    }
    SID_log("Done.",SID_LOG_CLOSE);
  } // loop over snaps

  // Write the remaining snapshots
  for(;j_write>=i_read_start;i_write--,j_write-=i_read_step,l_write++){
     check_for_fragmented_halos(0,subgroups,n_subgroups[l_write],i_write,j_write,n_wrap);
     check_for_fragmented_halos(1,groups,   n_groups[l_write],   i_write,j_write,n_wrap);
     write_trees_horizontal((void **)groups,   
                            (void **)subgroups,
                            n_groups[l_write],   n_groups_max,   0,
                            n_subgroups[l_write],n_subgroups_max,0,
                            n_subgroups_group,
                            max_tree_id_subgroup,
                            max_tree_id_group,
                            i_write,
                            j_write,
                            l_write,
                            n_wrap,
                            i_file_start,
                            filename_cat_root_in,
                            filename_output_dir,
                            a_list,
                            cosmo,
                            n_k_match,
                            TREE_HORIZONTAL_WRITE_EXTENDED|TREE_HORIZONTAL_WRITE_ALLCASES);
  }
  int i_write_last;
  int l_write_last;
  int j_write_last;
  i_write_last=i_write+1;
  j_write_last=j_write+i_read_step;
  l_write_last=l_write-1;

  // Clean-up
  SID_log("Freeing arrays...",SID_LOG_OPEN);
  for(i_search=0;i_search<n_wrap;i_search++){
     // Free subgroup information
     for(i_halo=0;i_halo<n_subgroups_max;i_halo++)
        SID_free(SID_FARG subgroups[i_search][i_halo].bridges);
     SID_free(SID_FARG subgroups[i_search]);

     // Free group information
     for(i_halo=0;i_halo<n_groups_max;i_halo++)
        SID_free(SID_FARG groups[i_search][i_halo].bridges);
     SID_free(SID_FARG groups[i_search]);
  }
  SID_free(SID_FARG subgroups);
  SID_free(SID_FARG groups);
  SID_free(SID_FARG match_id);
  SID_free(SID_FARG match_score);
  SID_free(SID_FARG match_index);
  SID_free(SID_FARG n_particles_groups);
  SID_free(SID_FARG n_particles_subgroups);
  SID_log("Done.",SID_LOG_CLOSE);

  // At this point, fragmented halos are only labeled when they appear.
  //    This will propagate the fragmented halo flags forward in time.
  if(flag_compute_fragmented)
     compute_trees_horizontal_fragmented(n_groups,
                                         n_subgroups,
                                         n_subgroups_group,
                                         i_file_start,
                                         i_write_last,
                                         j_write_last,
                                         l_write_last,
                                         i_read_stop,
                                         i_read_step,
                                         max_tree_id_subgroup,
                                         max_tree_id_group,
                                         n_subgroups_max,
                                         n_groups_max,
                                         n_search,
                                         n_files,
                                         n_wrap,
                                         n_k_match,
                                         a_list,
                                         cosmo,
                                         filename_output_dir);

  // Compute ghost-populated trees if we're asked to
  if(flag_compute_ghosts)
     compute_trees_horizontal_ghosts(n_groups,
                                     n_subgroups,
                                     n_subgroups_group,
                                     i_read_start,
                                     i_file_start,
                                     i_write_last,
                                     j_write_last,
                                     l_write_last,
                                     i_read_stop,
                                     i_read_step,
                                     max_tree_id_subgroup,
                                     max_tree_id_group,
                                     n_subgroups_max,
                                     n_groups_max,
                                     n_search,
                                     n_files,
                                     n_wrap,
                                     n_k_match,
                                     a_list,
                                     cosmo,
                                     filename_cat_root_in,
                                     filename_output_dir);

  // Some final clean-up
  SID_free(SID_FARG n_groups);
  SID_free(SID_FARG n_subgroups);
  for(i_search=0;i_search<n_wrap;i_search++)
     SID_free(SID_FARG n_subgroups_group[i_search]);
  SID_free(SID_FARG n_subgroups_group);

  // Set flag_clean=TRUE so that the output files generated here
  //  will be treated as temporary in the next step (if called)
  //(*flag_clean)=TRUE;

  SID_log("Done.",SID_LOG_CLOSE);
}

