#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpTrees_build.h>

void identify_bridges(tree_horizontal_info **halos,
                      tree_horizontal_info  *halos_i,
                      int     n_halos_i,
                      int    *match_id,
                      float  *match_score,
                      size_t *match_index,
                      char   *match_flag_two_way,
                      int    *n_particles,
                      int     i_file,
                      int     i_read,
                      int     i_read_start,
                      int     i_read_stop,
                      int     i_read_step,
                      int     n_search,
                      int     n_wrap,
                      int     n_halos_max,
                      int     n_files,
                      char   *filename_root_matches,
                      int     flag_match_subgroups){
    SID_log("Identifying bridge candidates from back-matching...",SID_LOG_OPEN|SID_LOG_TIMER);
    SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);

    // First, do an initial count of matches.  This will not be a list of unique halos
    //    though, since the same halos are likely to appear in repeated snapshots.
    int j_file_1;
    int j_file_2;
    int j_read_1;
    int j_read_2;
    int i_search;
    int i_halo;
    int j_halo;
    int k_halo;
    int l_halo;
    int n_halos_1_matches;
    int n_halos_2_matches;
    back_match_info *back_matches;
    back_match_info *back_match;
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
                    j_read_1,j_read_2,n_halos_max,
                    flag_match_subgroups,
                    &n_halos_1_matches,
                    &n_halos_2_matches,
                    NULL,
                    n_particles,
                    NULL,
                    NULL,
                    match_id,
                    match_score,
                    match_index,
                    NULL,
                    F_GOODNESS_OF_MATCH);

       // Store halo sizes for the current snapshot's halos
       if(i_search==0){
          for(i_halo=0;i_halo<n_halos_2_matches;i_halo++)
             halos[i_file%n_wrap][i_halo].n_particles=n_particles[i_halo];
       }

       // Perform initial back-match count
       for(i_halo=0;i_halo<n_halos_i;i_halo++){
          j_halo=find_index_int(match_id,i_halo,n_halos_1_matches,match_index);
          while(match_id[match_index[j_halo]]==i_halo && j_halo<(n_halos_1_matches-1)){
             halos_i[i_halo].n_back_matches++;
             j_halo++;
          }
          if(match_id[match_index[j_halo]]==i_halo && j_halo==(n_halos_1_matches-1)){
             halos_i[i_halo].n_back_matches++;
             j_halo++;
          }
       }
       SID_log("Done.",SID_LOG_CLOSE);
    } // loop over snaps

    //    ... second, do a conservative allocation using the non-unique counts and reset the counter.
    for(i_halo=0;i_halo<n_halos_i;i_halo++){
       if((halos_i[i_halo].n_back_matches)>0)
          (halos_i[i_halo].back_matches)=(back_match_info *)SID_calloc(sizeof(back_match_info)*(halos_i[i_halo].n_back_matches));
       else
          (halos_i[i_halo].back_matches)=NULL;
       halos_i[i_halo].n_back_matches =0;
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
                    j_read_1,j_read_2,n_halos_max,
                    flag_match_subgroups,
                    &n_halos_1_matches,
                    &n_halos_2_matches,
                    NULL,
                    NULL,
                    NULL,
                    NULL,
                    match_id,
                    match_score,
                    match_index,
                    NULL,
                    F_GOODNESS_OF_MATCH);

       // For all the halos in i_file_1 with back-matches ...
       int flag_continue;
       for(i_halo=0;i_halo<n_halos_i;i_halo++){
          if((halos_i[i_halo].back_matches)!=NULL){
             // Scan over the list of halos from snapshot=j_read_1 
             //   that match this halo in j_read_2 ...
             back_matches=halos_i[i_halo].back_matches;
             j_halo =find_index_int(match_id,i_halo,n_halos_1_matches,match_index);
             // Loop over all but the last halo in the list ...
             while(match_id[match_index[j_halo]]==i_halo && j_halo<(n_halos_1_matches-1)){
                // Check to see if this halo ID is already in the back match list (keep all -ve ID halos as well;
                //   a tree-walk descendant cull will be done below as well) ...
                flag_continue=TRUE;
                if(halos[j_file_1%n_wrap][match_index[j_halo]].id>=0){
                   for(k_halo=0;k_halo<halos_i[i_halo].n_back_matches && flag_continue;k_halo++){
                      if(back_matches[k_halo].halo->id==halos[j_file_1%n_wrap][match_index[j_halo]].id)
                         flag_continue=FALSE;
                   }
                }
                // ... if not, add it
                if(flag_continue){
                   back_matches[halos_i[i_halo].n_back_matches].halo =&(halos[j_file_1%n_wrap][match_index[j_halo]]);
                   back_matches[halos_i[i_halo].n_back_matches].file =j_file_1;
                   back_matches[halos_i[i_halo].n_back_matches].score=match_score[match_index[j_halo]];
                   (halos_i[i_halo].n_back_matches)++;
                }
                j_halo++;
             }
             // ... then do the last halo in the list ...
             if(match_id[match_index[j_halo]]==i_halo && j_halo==(n_halos_1_matches-1)){
                // Check to see if this halo ID is already in the back match list (keep all -ve ID halos as well;
                //   a tree-walk descendant cull will be done below as well) ...
                flag_continue=TRUE;
                if(halos[j_file_1%n_wrap][match_index[j_halo]].id>=0){
                   for(k_halo=0;k_halo<halos_i[i_halo].n_back_matches && flag_continue;k_halo++){
                      if(back_matches[k_halo].halo->id==halos[j_file_1%n_wrap][match_index[j_halo]].id)
                         flag_continue=FALSE;
                   }
                }
                // ... if not, add it
                if(flag_continue){
                   back_matches[halos_i[i_halo].n_back_matches].halo =&(halos[j_file_1%n_wrap][match_index[j_halo]]);
                   back_matches[halos_i[i_halo].n_back_matches].file =j_file_1;
                   back_matches[halos_i[i_halo].n_back_matches].score=match_score[match_index[j_halo]];
                   (halos_i[i_halo].n_back_matches)++;
                }
                j_halo++;
             }
          }
       }
       SID_log("Done.",SID_LOG_CLOSE);
    }

    // ... lastly, reorder the back-matched halos by the sizes of their most massive descendants.  Keep only 
    //        the most immediate bridge descendants and finalize the list ...
    SID_log("Re-ordering back_matches...",SID_LOG_OPEN|SID_LOG_TIMER);
    for(i_halo=0;i_halo<n_halos_i;i_halo++){
       // We may need to remove several halos from the list.  This array will keep track of this.
       size_t *back_match_index=NULL;
       int    *backmatch_keep  =(int *)SID_malloc(sizeof(int)*halos_i[i_halo].n_back_matches);

       // Reorder the back_matches by the size of their most massive descendant.  We make a temporary copy of the list 
       //   to do this and initially set all back_matches as halos to keep..
       back_matches=(back_match_info *)SID_calloc(sizeof(back_match_info)*(halos_i[i_halo].n_back_matches));
       for(j_halo=0;j_halo<halos_i[i_halo].n_back_matches;j_halo++){
          back_match=&(halos_i[i_halo].back_matches[j_halo]);
          memcpy(&(back_matches[j_halo]),back_match,sizeof(back_match_info));
          match_score[j_halo]   =(float)(back_match->halo->n_particles_largest_descendant);
          backmatch_keep[j_halo]=TRUE;
       }
       merge_sort((void *)match_score,(size_t)(halos_i[i_halo].n_back_matches),&back_match_index,SID_FLOAT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

       // Remove any mutual descendants from the list
       //   (since they may have different IDs and passed the test above)
       for(j_halo=0;j_halo<halos_i[i_halo].n_back_matches;j_halo++){
          int k_file;
          int l_file;
          // ... walk the tree upwards for each back matched halo...
          tree_horizontal_info *current;
          back_match = &(back_matches[j_halo]);
          current=back_match->halo->descendant.halo;
          if(current!=NULL)
             k_file=current->file;
          l_file=k_file;
          while(current!=NULL && k_file>=l_file && k_file<MIN(n_files,i_file+(n_search+1))){
             for(k_halo=0;k_halo<halos_i[i_halo].n_back_matches;k_halo++){
                if(j_halo!=k_halo){ // Don't waste time checking a halo against itself
                   back_match = &(back_matches[k_halo]);
                   if(back_match->halo==current)
                      backmatch_keep[k_halo]=FALSE;
                }
             }
             current=current->descendant.halo;
             l_file=k_file;
             if(current!=NULL)
                k_file =current->file;
          }
       }

       // Remove any back matches which have already been assigned to a halo.  This makes
       //    sure that the backmatch is uniquely set to the most immediate backmatched halo.
       for(j_halo=0;j_halo<halos_i[i_halo].n_back_matches;j_halo++){
          if(backmatch_keep[j_halo]){
             int backmatch_file =back_matches[j_halo].halo->file;
             int backmatch_index=back_matches[j_halo].halo->index;
             backmatch_keep[j_halo]=((halos[backmatch_file%n_wrap][backmatch_index].bridge_backmatch.halo)==NULL); 
          }
       }

       // Since we may have trimmed the list, recount the number remaining
       int n_list=halos_i[i_halo].n_back_matches;
       for(j_halo=n_list-1,halos_i[i_halo].n_back_matches=0;j_halo>=0;j_halo--){
          if(backmatch_keep[j_halo])
             halos_i[i_halo].n_back_matches++;
       }

       // Set bridge halo flag
       if(halos_i[i_halo].n_back_matches>1)
          halos_i[i_halo].type|=TREE_CASE_BRIDGED;
       else
          halos_i[i_halo].type&=(~TREE_CASE_BRIDGED);

       // Because we've overallocated previously, reallocate the back matched list here to save RAM ...
       SID_free(SID_FARG halos_i[i_halo].back_matches);
       (halos_i[i_halo].back_matches)=(back_match_info *)SID_calloc(sizeof(back_match_info)*(halos_i[i_halo].n_back_matches));

       // ... and copy the sorted temporary list to the permanent list (DESCENDING ORDER!)
       for(j_halo=n_list-1,l_halo=0;j_halo>=0;j_halo--){
          int j_halo_sorted=back_match_index[j_halo];
          if(backmatch_keep[j_halo_sorted]){
             memcpy(&(halos_i[i_halo].back_matches[l_halo]),&(back_matches[j_halo_sorted]),sizeof(back_match_info));
             halos[(back_matches[j_halo_sorted].halo->file)%n_wrap][back_matches[j_halo_sorted].halo->index].bridge_backmatch.halo =
                &(halos_i[i_halo]);
             halos[(back_matches[j_halo_sorted].halo->file)%n_wrap][back_matches[j_halo_sorted].halo->index].bridge_backmatch.score=
                back_matches[j_halo_sorted].score;
             l_halo++;
          }
       }

       // Clean-up
       SID_free(SID_FARG backmatch_keep);
       SID_free(SID_FARG back_match_index);
       SID_free(SID_FARG back_matches);
    }
    SID_log("Done.",SID_LOG_CLOSE);

    SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
    SID_log("Done.",SID_LOG_CLOSE);
}

