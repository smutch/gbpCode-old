#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void read_trees_match_scores(tree_info *trees,char *filename_SSimPL_dir,int mode){
  int i_read;
  int j_read;
  int i_snap;
  int j_snap;
  int k_snap;
  int i_neighbour;

  SID_log("Reading tree match scores...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Create the array we will return
  if(check_mode_for_flag(mode,READ_TREES_MATCH_SCORES_GROUPS)){
     trees->group_match_scores=(float **)SID_malloc(sizeof(float *)*trees->n_snaps);
     for(int i_snap=0;i_snap<trees->n_snaps;i_snap++)
        trees->group_match_scores[i_snap]=(float *)SID_calloc(sizeof(float)*trees->n_groups_snap_local[i_snap]);
  }
  if(check_mode_for_flag(mode,READ_TREES_MATCH_SCORES_SUBGROUPS)){
     trees->subgroup_match_scores=(float **)SID_malloc(sizeof(float *)*trees->n_snaps);
     for(int i_snap=0;i_snap<trees->n_snaps;i_snap++)
        trees->subgroup_match_scores[i_snap]=(float *)SID_calloc(sizeof(float)*trees->n_subgroups_snap_local[i_snap]);
  }
  float **match_score_groups_local   =trees->group_match_scores;
  float **match_score_subgroups_local=trees->subgroup_match_scores;

  // Process each snapshot in turn
  int    *nebr_idx_list_local;
  int    *file_idx_list_local;
  int    *list_init_local;
  size_t *file_idx_list_local_index;
  int     n_list_local;
  file_idx_list_local=(int *)SID_malloc(sizeof(int)*trees->max_n_subgroups_snap_local);
  nebr_idx_list_local=(int *)SID_malloc(sizeof(int)*trees->max_n_subgroups_snap_local);
  list_init_local    =(int *)SID_malloc(sizeof(int)*trees->max_n_subgroups_snap_local);
  for(i_snap=0;i_snap<(trees->n_snaps-1);i_snap++){ // Don't process the last snapshot since it will not have any descendants
     i_read=trees->snap_list[i_snap];
     SID_log("Processing snapshot %03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read);

     // Loop twice; once for subgroups and then once for groups
     for(int i_type=0;i_type<2;i_type++){
        // Set some group/subgroup specific things
        int             n_halos;
        tree_node_info *first_neighbour;
        char            group_text_prefix[5];
        int             flag_proceed=TRUE;
        float          *match_score_local;
        switch(i_type){
           case 0:
              flag_proceed=check_mode_for_flag(mode,READ_TREES_MATCH_SCORES_SUBGROUPS);
              if(flag_proceed){
                 sprintf(group_text_prefix,"sub");
                 n_halos            =trees->n_subgroups_snap_local[i_snap];
                 first_neighbour    =trees->first_neighbour_subgroups[i_snap];
                 match_score_local  =match_score_subgroups_local[i_snap];
              }
              break;
           case 1:
              flag_proceed=check_mode_for_flag(mode,READ_TREES_MATCH_SCORES_SUBGROUPS);
              if(flag_proceed){
                 sprintf(group_text_prefix,"");
                 n_halos            =trees->n_groups_snap_local[i_snap];
                 first_neighbour    =trees->first_neighbour_groups[i_snap];
                 match_score_local  =match_score_groups_local[i_snap];
              }
              break;
        }

        if(flag_proceed){

           // Initialize a validation array
           for(i_neighbour=0;i_neighbour<n_halos;i_neighbour++) 
              list_init_local[i_neighbour]=0;

           // Scan over the matching range
           for(j_snap=i_snap+1,k_snap=1;j_snap<trees->n_snaps && k_snap<=trees->n_search;j_snap++,k_snap++){
              j_read=trees->snap_list[j_snap];

              // Create a sorted list of locally stored progenitor-descendant pairs between these snapshots
              tree_node_info *current;
              n_list_local=0;
              current     =first_neighbour;
              while(current!=NULL){
                 if(current->descendant!=NULL){
                    int file_offset = current->descendant->snap_tree-i_snap;
                    if(file_offset==k_snap){
                       file_idx_list_local[n_list_local]=current->file_index;
                       nebr_idx_list_local[n_list_local]=current->neighbour_index;
                       n_list_local++;
                    }
                 }
                 else
                    list_init_local[current->neighbour_index]=-1;
                 current=current->next_neighbour;
              }
              merge_sort(file_idx_list_local,(size_t)n_list_local,&file_idx_list_local_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);

              // Set filenames
              char   filename_in_dir_snap[256];
              char   filename_in[256];
              sprintf(filename_in_dir_snap,"%s/trees/matches/%03d",           filename_SSimPL_dir,i_read);
              sprintf(filename_in,         "%s/%sgroup_matches_%03d_%03d.dat",filename_in_dir_snap,group_text_prefix,i_read,j_read);

              // Read header
              SID_fp fp_in;
              int    i_read_file;
              int    j_read_file;
              int    n_groups_i;
              int    n_groups_j;
              SID_fopen(filename_in,"r",&fp_in);
              SID_fread_all(&i_read_file,sizeof(int),1,&fp_in);
              SID_fread_all(&j_read_file,sizeof(int),1,&fp_in);
              SID_fread_all(&n_groups_i, sizeof(int),1,&fp_in);
              SID_fread_all(&n_groups_j, sizeof(int),1,&fp_in);
              if(i_read_file!=i_read || j_read_file!=j_read)
                 SID_trap_error("Invalid file numbers in an input match file (ie %d!=%d or %d!=%d in {%s})",ERROR_LOGIC,
                                i_read_file,i_read,j_read_file,j_read,filename_in);
              if(i_read_file!=trees->snap_list[i_snap] || j_read_file!=trees->snap_list[j_snap])
                 SID_trap_error("The wrong match file is being read (ie %d!=%d or %d!=%d)",ERROR_LOGIC,
                                i_read_file,trees->snap_list[i_snap],j_read_file,trees->snap_list[j_snap]);

              // Read matching data
              int k_read;
              int l_read;
              SID_fskip(sizeof(int),   n_groups_i,&fp_in);
              SID_fskip(sizeof(size_t),n_groups_i,&fp_in);
              for(k_read=0,l_read=0;k_read<n_groups_i;k_read++){
                 float match_score;
                 SID_fread_all(&match_score,sizeof(float),1,&fp_in);
                 if(l_read<n_list_local){
                    if(k_read==file_idx_list_local[file_idx_list_local_index[l_read]]){
                       match_score_local[nebr_idx_list_local[file_idx_list_local_index[l_read]]]=match_score;
                       list_init_local[nebr_idx_list_local[file_idx_list_local_index[l_read]]]++;
                       l_read++;
                    }
                 }
              }
              if(l_read!=n_list_local)
                 SID_trap_error("An incorrect number of matches were read (ie. %d!=%d) for snaphot %d->%d",ERROR_LOGIC,l_read,n_list_local,i_read,j_read);
              SID_fclose(&fp_in);
              SID_free(SID_FARG file_idx_list_local_index);
           }
        }

        // Check that all neighbours have been processed correctly
        for(i_neighbour=0;i_neighbour<n_halos;i_neighbour++)
           if(list_init_local[i_neighbour]!=1 && list_init_local[i_neighbour]>=0)
              SID_trap_error("Neighbour %d of %d was not processed correctly (init=%d) for snapshot %d.",ERROR_LOGIC,
                             i_neighbour,n_halos,list_init_local[i_neighbour],i_read);
     }
     SID_log("Done.",SID_LOG_CLOSE);
  }
  SID_free(SID_FARG nebr_idx_list_local);
  SID_free(SID_FARG file_idx_list_local);
  SID_free(SID_FARG list_init_local);

  SID_log("Done.",SID_LOG_CLOSE);
}

