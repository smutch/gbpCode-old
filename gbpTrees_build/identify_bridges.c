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
    bridge_info *bridges;
    bridge_info *bridge;
    int         *bridge_keep=NULL;
    size_t      *bridge_index=NULL;
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
                    match_index);

       // For all the halos in i_file_1 with back-matches ...
       int flag_continue;
       for(i_halo=0;i_halo<n_halos_i;i_halo++){
          if((halos_i[i_halo].bridges)!=NULL){
             // Scan over the list of halos from snapshot=j_read_1 
             //   that match this halo in j_read_2 ...
             
             bridges=halos_i[i_halo].bridges;
             j_halo =find_index_int(match_id,i_halo,n_halos_1_matches,match_index);
             // Loop over all but the last halo in the list ...
             while(match_id[match_index[j_halo]]==i_halo && j_halo<(n_halos_1_matches-1)){
                // Check to see if this halo is already in the bridge list (keep all strayed halos as well) ...
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
                // Check to see if this halo is already in the list (keep all strayed halos as well) ...
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
             int k_file;
             int l_file;
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
          int n_list;
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

