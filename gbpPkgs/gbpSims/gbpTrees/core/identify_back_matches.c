#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpTrees_build.h>

void identify_back_matches(tree_horizontal_info **halos,
                           tree_horizontal_info  *halos_i,
                           match_info           **back_matches_i,
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
    SID_log("Identifying back-matches...",SID_LOG_OPEN|SID_LOG_TIMER);
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
             halos_i[i_halo].n_particles=n_particles[i_halo];
       }

       // Perform initial back-match count
       for(i_halo=0;i_halo<n_halos_2_matches;i_halo++){
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
    //        Add one to the allocation for a possible forematch at well.  Re-zero the count so that
    //        we can perform a more accurate count below.  We use a temporary array for this since
    //        we will have to reallocate after trimming the list.
    int n_allocate_backmatch=0;
    for(i_halo=0;i_halo<n_halos_2_matches;i_halo++)
       n_allocate_backmatch+=(1+halos_i[i_halo].n_back_matches);
    match_info *back_matches_i_temp=(match_info *)SID_calloc(sizeof(match_info)*n_allocate_backmatch);
    for(i_halo=0,n_allocate_backmatch=0;i_halo<n_halos_2_matches;i_halo++){
       halos_i[i_halo].back_matches  =&(back_matches_i_temp[n_allocate_backmatch]);
       n_allocate_backmatch         +=(1+halos_i[i_halo].n_back_matches);
       halos_i[i_halo].n_back_matches=0;
    }

    //    ... third, assemble the lists of first fore and back matches.
    for(j_file_1  =i_file+1,
          j_file_2=i_file,
          j_read_1=i_read+i_read_step,
          j_read_2=i_read,
          i_search=0;
        j_read_1<=i_read_stop && i_search<n_search;
        j_file_1++,
          j_read_1+=i_read_step,
          i_search++){

       SID_log("Finding matches between files %d->%d...",SID_LOG_OPEN,j_read_1,j_read_2);

       // Read fore matches
       read_matches(filename_root_matches,
                    j_read_2,j_read_1,n_halos_max,
                    flag_match_subgroups,
                    &n_halos_2_matches,
                    &n_halos_1_matches,
                    NULL,
                    NULL,
                    NULL,
                    NULL,
                    match_id,
                    match_score,
                    match_index,
                    match_flag_two_way,
                    F_GOODNESS_OF_MATCH);

       if(n_halos_1_matches>0 && n_halos_2_matches>0){
          // Assemble first fore matches for each halo
          for(i_halo=0;i_halo<n_halos_2_matches;i_halo++){
             j_halo=match_id[i_halo];
             if(j_halo>=0){
                tree_horizontal_info *halo_i=&(halos[j_file_2%n_wrap][i_halo]);
                if(halo_i->forematch_first.halo==NULL){
                   tree_horizontal_info *halo_j           =&(halos[j_file_1%n_wrap][j_halo]);
                   halo_i->forematch_first.halo           =halo_j;
                   halo_i->forematch_first.score          =match_score[i_halo];
                   halo_i->forematch_first.flag_two_way   =match_flag_two_way[i_halo];
                   halo_i->forematch_first.flag_back_match=FALSE;
                }
             }
          }

          // Read back matches
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
                       match_flag_two_way,
                       F_GOODNESS_OF_MATCH);

          // Scan over all the current snapshot's halos ...
          for(j_halo=0;j_halo<n_halos_1_matches;j_halo++){
             // Count this halo as a back match ...
             //    if there is a good match and this halo has not 
             //    previously been assigned as a backmatch to an
             //    earlier halo.  This makes sure that the backmatch
             //    is uniquely set to the most immediate backmatched 
             //    halo.  This is important (for example) in cases where 
             //    a halo is repeatedly emerging from a bridge or (even 
             //    more importantly) when an emerged halo is lost for 
             //    several snapshots.  We only want to use this as an 
             //    emerged candidate for one (the most immediate) instance,
             //    just when the halo actually emerges, and not before.
             tree_horizontal_info *halo_j=&(halos[j_file_1%n_wrap][j_halo]);
             i_halo=match_id[j_halo];
             if(i_halo>=0 && (halo_j->bridge_backmatch.halo)==NULL){
                tree_horizontal_info *halo_i        =&(halos_i[i_halo]);
                match_info           *back_matches  =halo_i->back_matches;
                int                   n_back_matches=halo_i->n_back_matches;
                // Check to see if this halo's ID has already been added
                int halo_j_id=halo_j->id;
                int flag_add =TRUE;
                if(halo_j_id>=0){
                   for(k_halo=0;k_halo<n_back_matches && flag_add;k_halo++){
                      if(back_matches[k_halo].halo->id==halo_j_id)
                         flag_add=FALSE;
                   }
                }
                // Add this halo to the preliminary list if it passed the above test(s)
                if(flag_add){
                   back_matches[n_back_matches].halo           =halo_j;
                   back_matches[n_back_matches].score          =match_score[j_halo];
                   back_matches[n_back_matches].flag_two_way   =match_flag_two_way[j_halo];
                   back_matches[n_back_matches].flag_back_match=TRUE;
                   halo_i->n_back_matches++;
                }
             }
          }
       }
       SID_log("Done.",SID_LOG_CLOSE);
    }

    // Insert first forematches into the backmatch list
    SID_log("Inserting first forematches into the back match list...",SID_LOG_OPEN|SID_LOG_TIMER);
    for(i_halo=0;i_halo<n_halos_2_matches;i_halo++){
       // If we have found a forematch for this halo ...
       tree_horizontal_info *halo_i=&(halos_i[i_halo]);
       if(halo_i->forematch_first.halo!=NULL){
          // Look to see if it is present in the current list of backmatches
          match_info *back_matches =halo_i->back_matches;
          for(k_halo=0;k_halo<halo_i->n_back_matches;k_halo++){
             if(halo_i->forematch_first.halo==back_matches[k_halo].halo) break;
          }
          // We are either overwritting a back match with a forematch or
          //    adding it to the end of the list here
          memcpy(&(back_matches[k_halo]),&(halo_i->forematch_first),sizeof(match_info));
          if(k_halo==halo_i->n_back_matches)
             halo_i->n_back_matches++;
       }
    }
    SID_log("Done.",SID_LOG_CLOSE);

    // Find the largest temporary back match count.  This will
    //    be used for allocating a couple of temporary arrays.
    int n_allocate_max=0;
    for(i_halo=0;i_halo<n_halos_2_matches;i_halo++)
       n_allocate_max=MAX(n_allocate_max,halos_i[i_halo].n_back_matches);

    // ... lastly, finalize the list, keeping only the most immediate descendants.
    SID_log("Finalizing list...",SID_LOG_OPEN|SID_LOG_TIMER);
    int        *backmatch_keep=(int        *)SID_malloc(sizeof(int)       *n_allocate_max);
    match_info *back_matches  =(match_info *)SID_calloc(sizeof(match_info)*n_allocate_max);
    match_info *empty_match   =(match_info *)SID_calloc(sizeof(match_info));
    for(i_halo=0;i_halo<n_halos_2_matches;i_halo++){
       // Zero temporary arrays
       for(j_halo=0;j_halo<halos_i[i_halo].n_back_matches;j_halo++){
          backmatch_keep[j_halo]=FALSE;
          back_matches[j_halo]  =empty_match[0];
       }

       // Reorder the back_matches.  We make a temporary copy of the list 
       //   to do this and initially set all back_matches as halos to keep..
       size_t *back_match_index=NULL;
       for(j_halo=0;j_halo<halos_i[i_halo].n_back_matches;j_halo++){
          match_info *back_match=&(halos_i[i_halo].back_matches[j_halo]);
          memcpy(&(back_matches[j_halo]),back_match,sizeof(match_info));
          // Set the criteria by which the halos will be sorted here
          match_score[j_halo]   =(float)(back_match->halo->n_particles);
          backmatch_keep[j_halo]=TRUE;
       }
       merge_sort((void *)match_score,(size_t)(halos_i[i_halo].n_back_matches),&back_match_index,SID_FLOAT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

       // IMPORTANT: Remove any mutual descendants from the list
       //   (since they may have different IDs and passed the test above).
       //   Keep the most immediate instance.
       // n.b.: As it stands, this is a order-N^2 process.  Sorting by snapshot
       //       first (and then resorting by size later) could turn this into an 
       //       order-NlogN process.
       for(j_halo=0;j_halo<halos_i[i_halo].n_back_matches;j_halo++){
          int k_file;
          int l_file;
          // ... start with each back match's descendant...
          match_info           *back_match=&(back_matches[j_halo]);
          tree_horizontal_info *current   =back_match->halo->descendant.halo;
          if(current!=NULL)
             k_file=current->file;
          l_file=k_file;
          // ... walk the tree upwards for each back matched halo...
          while(current!=NULL && k_file>=l_file && k_file<MIN(n_files,i_file+(n_search+1))){
             // ... loop over the other back matches ...
             for(k_halo=0;k_halo<halos_i[i_halo].n_back_matches;k_halo++){
                // ... don't waste time checking a halo against itself or against one already removed
                if(j_halo!=k_halo && backmatch_keep[k_halo]){ 
                   // ... if we've walked into another back match, remove it.  This
                   //     ensures that we are removing the least immediate instances.
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

       // Since we may have trimmed the list, recount the number remaining
       int n_list=halos_i[i_halo].n_back_matches;
       for(j_halo=n_list-1,halos_i[i_halo].n_back_matches=0;j_halo>=0;j_halo--){
          if(backmatch_keep[j_halo])
             halos_i[i_halo].n_back_matches++;
       }

       // Slide the sorted temporary list to the permanent list (DESCENDING ORDER!)
       for(j_halo=n_list-1,l_halo=0;j_halo>=0;j_halo--){
          int j_halo_sorted=back_match_index[j_halo];
          if(backmatch_keep[j_halo_sorted]){
             match_info *back_match=&(back_matches[j_halo_sorted]);
             memcpy(&(halos_i[i_halo].back_matches[l_halo]),back_match,sizeof(match_info));
             // This check is true if the match is a back match or is 2-way ... we want to
             //    ensure that only the most immediate back match is used for each halo.
             //    This is achieved if we set bridge_backmatch.halo!=NULL
             if(back_match->flag_two_way || back_match->flag_back_match){
                halos[(back_match->halo->file)%n_wrap][back_match->halo->index].bridge_backmatch.halo =
                   &(halos_i[i_halo]);
                halos[(back_match->halo->file)%n_wrap][back_match->halo->index].bridge_backmatch.score=
                   back_match->score;
                halos[(back_match->halo->file)%n_wrap][back_match->halo->index].bridge_backmatch.flag_two_way=
                   back_match->flag_two_way;
                halos[(back_match->halo->file)%n_wrap][back_match->halo->index].bridge_backmatch.flag_back_match=
                   back_match->flag_back_match;
             }
             l_halo++;
          }
       }

       // Clean-up
       SID_free(SID_FARG back_match_index);
    }
    // Clean-up
    SID_free(SID_FARG backmatch_keep);
    SID_free(SID_FARG back_matches);
    SID_free(SID_FARG empty_match);

    // Because we've (likely) overallocated previously, reallocate and repopulate the back match array here to save RAM ...
    n_allocate_backmatch=0;
    for(i_halo=0;i_halo<n_halos_2_matches;i_halo++)
       n_allocate_backmatch+=halos_i[i_halo].n_back_matches;
    (*back_matches_i)=(match_info *)SID_calloc(sizeof(match_info)*n_allocate_backmatch);
    for(i_halo=0,k_halo=0;i_halo<n_halos_2_matches;i_halo++){
       l_halo=k_halo;
       tree_horizontal_info *halo_i=&(halos_i[i_halo]);
       for(j_halo=0;j_halo<halos_i[i_halo].n_back_matches;j_halo++,k_halo++)
          memcpy(&((*back_matches_i)[k_halo]),&(halo_i->back_matches[j_halo]),sizeof(match_info));
       halo_i->back_matches=&((*back_matches_i)[l_halo]);
    }
    SID_free(SID_FARG back_matches_i_temp);

    SID_log("Done.",SID_LOG_CLOSE);

    // !! n.b.: The TREE_CASE_BRIDGED and TREE_CASE_EMERGED_CANDIDATE flags will be set later in identify_bridges(), 
    //          once all main progenitor decisions have been made for this snapshot.  The chosen main progenitor 
    //          will not be given this flag and only subsequent matches to the others will be considered as 
    //          possible emerged halos.

    SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
    SID_log("Done.",SID_LOG_CLOSE);
}
