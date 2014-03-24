#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void compute_trees_horizontal(char        *filename_halo_root_in,
                              char        *filename_cat_root_in,
                              char        *filename_snap_list_in,
                              char        *filename_root_matches,
                              char        *filename_output_dir,
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
  int         max_id         =0;
  int         max_id_group   =0;
  int         max_id_subgroup=0;
  int        *my_descendant;
  int        *n_particles;
  int        *n_particles_groups;
  int        *n_particles_subgroups;
  int         my_trunk;
  double      expansion_factor;
  int         n_found;
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
  int     n_halos_1_matches;
  int     n_halos_2_matches;
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
  int     i_file_start;
  
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
  int flag_compute_ghosts    =FALSE;

  if(!flag_fix_bridges)
    SID_log("Bridge-fixing is turned off.",SID_LOG_COMMENT);
  if(!flag_compute_fragmented)
    SID_log("Fragmented-halo propagation is turned off.",SID_LOG_COMMENT);
  if(!flag_compute_ghosts)
    SID_log("Ghost-populated tree construction is turned off.",SID_LOG_COMMENT);

  // Create the output directory
  mkdir(filename_output_dir,02755);

  // Create snapshot expansion factor list
  double *a_list=NULL;
  int     n_a_list_in;
  write_a_list(filename_snap_list_in,
               filename_output_dir,
               i_read_start,
               i_read_stop,
               i_read_step);
  read_a_list(filename_output_dir,
              &a_list,
              &n_a_list_in);

  write_tree_run_parameters(filename_output_dir,
                            i_read_start,
                            i_read_stop,
                            i_read_step,
                            n_search,
                            flag_fix_bridges,
                            flag_compute_fragmented,
                            flag_compute_ghosts);

  // Validate existing matching files &/or perfrom matching
  //if(!compute_trees_matches(filename_halo_root_in,
  //                          filename_root_matches,
  //                          i_read_start,
  //                          i_read_stop,
  //                          i_read_step,
  //                          n_search,
  //                          WRITE_MATCHES_MODE_TREES|WRITE_MATCHES_PERFORM_CHECK))
  //   SID_trap_error("Matching could not be completed.  Terminating.",ERROR_LOGIC);
  read_matches_header(filename_root_matches,
                      i_read_start,
                      i_read_stop,
                      i_read_step,
                      &n_files,
                      &n_subgroups,
                      &n_groups,
                      &n_subgroups_max,
                      &n_groups_max,
                      &n_halos_max);

  // We need these for allocating arrays
  calc_max(n_subgroups,&n_subgroups_max,n_files,SID_INT,CALC_MODE_DEFAULT);
  calc_max(n_groups,   &n_groups_max,   n_files,SID_INT,CALC_MODE_DEFAULT);
  n_halos_max=MAX(n_subgroups_max,n_groups_max);

  // We need enough indices to allow us to hold-on to descendants until outputting
  //   and for the current and last i_file as well
  n_wrap      =2*n_search+2;
  i_file_start=n_files-1;
     
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
  init_trees_horizontal_roots(groups,
                              subgroups,
                              match_id,
                              match_score,
                              match_index,
                              n_particles_groups,
                              n_particles_subgroups,
                              n_subgroups_group,
                              n_groups_max,
                              n_subgroups_max,
                              filename_root_matches,
                              i_read_start,
                              i_read_stop,
                              i_read_step,
                              i_file_start,
                              n_wrap,
                              n_halos_max,
                              &max_id_group,
                              &max_tree_id_group,
                              &max_id_subgroup,
                              &max_tree_id_subgroup);

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

       // Initialize a bunch of stuff which depends on whether
       //   we are processing groups or subgroups
       switch(k_match){
          case 0:
          sprintf(group_text_prefix,"sub");
          flag_match_subgroups=MATCH_SUBGROUPS;
          halos               =subgroups;
          n_halos             =n_subgroups;
          n_halos_max         =n_subgroups_max;
          max_id              =max_id_subgroup;
          max_tree_id         =max_tree_id_subgroup;
          n_particles         =n_particles_subgroups;
          break;
          case 1:
          sprintf(group_text_prefix,"");
          flag_match_subgroups=MATCH_GROUPS;
          halos               =groups;
          n_halos             =n_groups;
          n_halos_max         =n_groups_max;
          max_id              =max_id_group;
          max_tree_id         =max_tree_id_group;
          n_particles         =n_particles_groups;
          break;
       }
       halos_i  =halos[i_file%n_wrap];
       n_halos_i=n_halos[j_file];

       SID_log("Processing %d %sgroups...",SID_LOG_OPEN|SID_LOG_TIMER,n_halos_i,group_text_prefix);

       // Initialize tree pointer-arrays with dummy values
       init_trees_horizontal_snapshot(halos_i,n_halos_i,i_read,i_file,n_halos_max);
       
       // Use back-matching to identify bridged halos ...
       if(flag_fix_bridges)
          identify_bridges(halos,
                           halos_i,
                           n_halos_i,
                           match_id,
                           match_score,
                           match_index,
                           n_particles,
                           i_file,
                           i_read,
                           i_read_start,
                           i_read_stop,
                           i_read_step,
                           n_search,
                           n_wrap,
                           n_halos_max,
                           n_files,
                           filename_root_matches,
                           flag_match_subgroups);

       // Perform forward-matching
       construct_progenitors(halos,
                             halos_i,
                             n_subgroups_group,
                             n_halos_i,
                             match_id,
                             match_score,
                             match_index,
                             n_particles,
                             i_file,
                             i_read,
                             i_read_start,
                             i_read_stop,
                             i_read_step,
                             n_search,
                             n_wrap,
                             n_halos_max,
                             n_files,
                             flag_fix_bridges,
                             &max_id,
                             &n_halos_1_matches,
                             &n_halos_2_matches,
                             filename_root_matches,
                             group_text_prefix,
                             flag_match_subgroups);
      
       // Finalize matches to unprocessed halos ...
       apply_horizontal_tree_defaults(n_halos_1_matches,
                                      n_halos_i,
                                      halos,
                                      halos_i,
                                      i_file,
                                      n_wrap,
                                      &max_id,
                                      &max_tree_id);
 
       // Now that we have assigned all the IDs for the halos in the active snapshot,
       //   we need to remove all descendants of bridged halos from the lists of candidate emerged halos 
       //   to avoid incorrectly matching to the main progenitor's descendants later-on when when 
       //   we are scaning emerged halo candidates.  Real matches to bridges are dealt-with
       //   when halos marked TREE_CASE_BRIDGE_DEFAULT are processed. 
       clean_emerged_halo_list(halos_i,
                               n_halos_i,
                               i_file,
                               n_search,
                               n_files);

       // This has to be written right after a snapshot is read and processed (because it needs all forward scan information), 
       //    so it is separate from the rest of the log output code.
       write_trees_horizontal_emerged_candidates(i_read,n_halos_i,halos_i,group_text_prefix,filename_output_dir,j_file==1);

       // Report some statistics
       //   n.b.: This is only an estimate in some cases, since subsequent snapshots may alter this snapshot.  
       //         See the final written log.txt file for accurate numbers.
       write_trees_horizontal_report(n_halos_i,n_halos_max,halos_i);

       // Update the max_id variables
       switch(k_match){
          case 0:
            max_id_subgroup     =max_id;
            max_tree_id_subgroup=max_tree_id;
            break;
          case 1:
            max_id_group     =max_id;
            max_tree_id_group=max_tree_id;
            break;
       }
       SID_log("Done.",SID_LOG_CLOSE);
    } // k_match
 
    // Write trees once a few files have been processed
    //   and no more dropped groups etc. need to be given ids
    if(j_file>n_search){
       int mode_write;
       if(flag_compute_ghosts || flag_compute_fragmented)
          mode_write=TREE_HORIZONTAL_WRITE_EXTENDED|TREE_HORIZONTAL_WRITE_ALLCASES|TREE_HORIZONTAL_WRITE_CHECK_FRAGMENTED;
       else
          mode_write=TREE_HORIZONTAL_WRITE_ALLCASES|TREE_HORIZONTAL_WRITE_CHECK_FRAGMENTED;
       write_trees_horizontal((void **)groups, 
                              (void **)subgroups,
                              n_groups[l_write],   n_groups_max,   
                              n_subgroups[l_write],n_subgroups_max,
                              n_subgroups_group,
                              max_tree_id_subgroup,
                              max_tree_id_group,
                              i_write,
                              j_write,
                              l_write,
                              i_read_step,
                              n_search,
                              n_wrap,
                              i_file_start,
                              filename_cat_root_in,
                              filename_output_dir,
                              a_list,
                              cosmo,
                              n_k_match,
                              l_write==0,
                              mode_write);
       i_write--;
       l_write++;
       j_write-=i_read_step;
    }
    SID_log("Done.",SID_LOG_CLOSE);

  } // loop over snaps

  // Write the remaining snapshots
  for(;j_write>=i_read_start;i_write--,j_write-=i_read_step,l_write++){
     int mode_write;
     if(flag_compute_ghosts || flag_compute_fragmented)
        mode_write=TREE_HORIZONTAL_WRITE_EXTENDED|TREE_HORIZONTAL_WRITE_ALLCASES|TREE_HORIZONTAL_WRITE_CHECK_FRAGMENTED;
     else
        mode_write=TREE_HORIZONTAL_WRITE_ALLCASES|TREE_HORIZONTAL_WRITE_CHECK_FRAGMENTED;
     write_trees_horizontal((void **)groups,   
                            (void **)subgroups,
                            n_groups[l_write],   n_groups_max,   
                            n_subgroups[l_write],n_subgroups_max,
                            n_subgroups_group,
                            max_tree_id_subgroup,
                            max_tree_id_group,
                            i_write,
                            j_write,
                            l_write,
                            i_read_step,
                            n_search,
                            n_wrap,
                            i_file_start,
                            filename_cat_root_in,
                            filename_output_dir,
                            a_list,
                            cosmo,
                            n_k_match,
                            l_write==0,
                            mode_write);
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

  // If extended horizontal tree files were written for fragmented
  //    halo propagation or ghost tree construction, remove them.
  if(flag_compute_ghosts || flag_compute_fragmented){
     SID_log("Removing temporary tree files...",SID_LOG_OPEN);
     for(j_write=i_read_stop;j_write>=i_read_start;j_write-=i_read_step){
        char filename_output_dir_horizontal[MAX_FILENAME_LENGTH];
        char filename_output_dir_horizontal_trees[MAX_FILENAME_LENGTH];
        char filename_remove[MAX_FILENAME_LENGTH];
        sprintf(filename_output_dir_horizontal,      "%s/horizontal",filename_output_dir);
        sprintf(filename_output_dir_horizontal_trees,"%s/trees",     filename_output_dir_horizontal);
        sprintf(filename_remove,"%s/horizontal_trees_tmp_%03d.dat",filename_output_dir_horizontal_trees,j_write);
        remove(filename_remove);
     }
     SID_log("Done.",SID_LOG_CLOSE);
  }

  // Some final clean-up
  SID_free(SID_FARG n_groups);
  SID_free(SID_FARG n_subgroups);
  for(i_search=0;i_search<n_wrap;i_search++)
     SID_free(SID_FARG n_subgroups_group[i_search]);
  SID_free(SID_FARG n_subgroups_group);
  SID_free(SID_FARG a_list);

  // Force the forest construction to use all snapshots
  int n_search_forests=i_read_stop;

  // Construct tree->forest mappings
  compute_forests(filename_output_dir,n_search_forests);

  SID_log("Done.",SID_LOG_CLOSE);
}

