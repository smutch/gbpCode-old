#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void match_halos(plist_info  *plist_1_in,
                 int          i_file_1_in,
                 int         *mark_list_1,
                 int          n_mark_1_local,
                 plist_info  *plist_2_in,
                 int          i_file_2_in,
                 int         *mark_list_2_local,
                 int          n_mark_2_local,
                 char        *catalog_1to2,
                 int          mode){
  plist_info  *plist_store;
  plist_info  *plist_1;
  int          i_file_1;
  plist_info  *plist_2;
  int          i_file_2;
  int     i_group;
  int     i_mark;
  size_t  i_particle;
  size_t  j_particle;
  int     k_particle;
  int     i_hist;
  int     j_hist;
  size_t  n_particles_1;
  size_t  n_particles_1_all;
  size_t  n_particles_2;
  size_t  n_particles_2_local;
  size_t  n_particles_2_all;
  size_t  n_groups_1_dim;
  size_t  n_mark_1_dim;
  int     n_mark_1_all;
  int     n_groups_1;
  int     n_groups_1_all;
  size_t  n_groups_2_s;
  size_t  n_groups_2_local_s;
  int     n_groups_2;
  int     n_groups_2_local;
  int     n_groups_2_all;
  int     n_match;
  int     n_match_all;
  int     n_match_1;
  int     n_match_1_all;
  char    parm_txt[256];
  int    *match;
  int    *n_particles_group_1;
  int    *n_particles_group_2;
  int    *n_particles_group_2_local;
  int    *group_offset_1;
  int     rank_offset_1;
  int    *group_offset_2;
  int    *group_offset_2_local;
  int     rank_offset_2;
  int     rank_offset_2_local;
  size_t *mark_list_index_1      =NULL;
  size_t *mark_list_index_2      =NULL;
  size_t *mark_list_index_2_local=NULL;
  size_t *index_1;
  size_t *index_2;
  size_t *id_1;
  size_t *id_2;
  size_t *id_2_local;
  size_t  id_1_i;
  int    *group_index_1;
  int    *group_index_2;
  int    *group_index_2_local;
  int     n_hist;
  int    *hist_list;
  float  *hist_score;
  float  *match_score;
  int    *match_rank;
  float   rank;
  int     i_rank;
  int     flag;
  size_t  group;
  size_t  hist_size=100;
  int     counter;
  int     n_p_match;
  double  f_min=0.9;
  int     flag_store_score;
  int     flag_read_marked;
  int     flag_read_this;
  char    catalog_1[5];
  char    catalog_2[5];
  int     flag_match_subgroups;
  int     flag_match_substructure;
  int     flag_match_continue;
  int     rank_offset;

  //SID_profile_start("match_groups",SID_PROFILE_NOTMPIENABLED);
  //SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

  // Check if we are back-matching.  If we are, then switch the inputs.
  if(check_mode_for_flag(mode,MATCH_BACK)){
    plist_1 =plist_2_in;
    i_file_1=i_file_2_in;
    plist_2 =plist_1_in;
    i_file_2=i_file_1_in; 
  }
  else{
    plist_1 =plist_1_in;
    i_file_1=i_file_1_in;
    plist_2 =plist_2_in;
    i_file_2=i_file_2_in;     
  }

  // Create string versions of the snapshot numbers
  //   (used for setting data labels in storage)
  sprintf(catalog_1,"%03d",i_file_1);
  sprintf(catalog_2,"%03d",i_file_2);

  // Determine if we are matching groups or subgroups (default is subgroups)
  if(check_mode_for_flag(mode,MATCH_GROUPS) && check_mode_for_flag(mode,MATCH_SUBGROUPS))
    SID_trap_error("match_groups can't match both groups and subgroups (yet).",ERROR_LOGIC);
  else if(check_mode_for_flag(mode,MATCH_GROUPS))
    flag_match_subgroups=FALSE;
  else if(check_mode_for_flag(mode,MATCH_SUBGROUPS))
    flag_match_subgroups=TRUE;
  else
    flag_match_subgroups=TRUE;

  // Are we matching substructure (we will ignore self matches and
  //   select the largest of multiple matches if so)
  if(check_mode_for_flag(mode,MATCH_SUBSTRUCTURE)){
    flag_match_subgroups   =FALSE;
    flag_match_substructure=TRUE;
  }
  else
    flag_match_substructure=FALSE;

  // Fetch some info and print log message
  if(flag_match_subgroups){
    SID_log("Matching subgroups from catalog {%s} to catalog {%s}...",
            SID_LOG_OPEN,
            catalog_1,
            catalog_2);
    if(ADaPS_exist(plist_1->data,"n_subgroups_%s",catalog_1))
      n_groups_1=((int *)ADaPS_fetch(plist_1->data,"n_subgroups_%s",catalog_1))[0];
    else
      n_groups_1=0;
    if(ADaPS_exist(plist_1->data,"n_subgroups_all_%s",catalog_1))
      n_groups_1_all=((int *)ADaPS_fetch(plist_1->data,"n_subgroups_all_%s",catalog_1))[0];
    else
      n_groups_1_all=0;
    if(ADaPS_exist(plist_2->data,"n_subgroups_%s",catalog_2))
      n_groups_2_local=((int *)ADaPS_fetch(plist_2->data,"n_subgroups_%s",catalog_2))[0];
    else
      n_groups_2_local=0;
    if(ADaPS_exist(plist_2->data,"n_subgroups_all_%s",catalog_2))
      n_groups_2_all=((int *)ADaPS_fetch(plist_2->data,"n_subgroups_all_%s",catalog_2))[0];
    else
      n_groups_2_all=0;
  }
  else{
    if(check_mode_for_flag(mode,MATCH_SUBSTRUCTURE))
      SID_log("Finding substructure from catalog {%s} in catalog {%s}...",
              SID_LOG_OPEN,
              catalog_1,
              catalog_2);
    else
      SID_log("Matching groups catalog {%s} to catalog {%s}...",
              SID_LOG_OPEN,
              catalog_1,
              catalog_2);
    if(ADaPS_exist(plist_1->data,"n_groups_%s",catalog_1))
      n_groups_1=((int *)ADaPS_fetch(plist_1->data,"n_groups_%s",catalog_1))[0];
    else
      n_groups_1=0;
    if(ADaPS_exist(plist_1->data,"n_groups_all_%s",catalog_1))
      n_groups_1_all=((int *)ADaPS_fetch(plist_1->data,"n_groups_all_%s",catalog_1))[0];
    else
      n_groups_1_all=0;
    if(ADaPS_exist(plist_2->data,"n_groups_%s",catalog_2))
      n_groups_2_local=((int *)ADaPS_fetch(plist_2->data,"n_groups_%s",catalog_2))[0];
    else
      n_groups_2_local=0;
    if(ADaPS_exist(plist_2->data,"n_groups_all_%s",catalog_2))
      n_groups_2_all=((int *)ADaPS_fetch(plist_2->data,"n_groups_all_%s",catalog_2))[0];
    else
      n_groups_2_all=0;
  }
  
  // Allow a mark list to be passed if we are
  //   only interested in certain objects
if(SID.I_am_Master) fprintf(stderr,"test1\n"); SID_Barrier(SID.COMM_WORLD);
  calc_sum_global(&n_mark_1_local,&n_mark_1_all,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
if(SID.I_am_Master) fprintf(stderr,"test2\n"); SID_Barrier(SID.COMM_WORLD);
  if(n_mark_1_all>0){
    flag_read_marked=TRUE;
    n_match_1       =n_mark_1_local;
    n_match_1_all   =n_mark_1_all;
  }
  else{
    flag_read_marked=FALSE;
    n_match_1       =n_groups_1;
    n_match_1_all   =n_groups_1_all;
  }

  // Perform matching if there are groups to match
if(SID.I_am_Master) fprintf(stderr,"test3\n"); SID_Barrier(SID.COMM_WORLD);
  if(n_groups_1_all>0 && n_groups_2_all>0){                

if(SID.I_am_Master) fprintf(stderr,"test4\n"); SID_Barrier(SID.COMM_WORLD);
    // Fetch needed info for catalog_1
    n_particles_1_all  =((size_t *)ADaPS_fetch(plist_1->data,"n_particles_all_%s",catalog_1))[0];
    n_particles_1      =((size_t *)ADaPS_fetch(plist_1->data,"n_particles_%s",    catalog_1))[0];
    id_1               = (size_t *)ADaPS_fetch(plist_1->data,"particle_ids_%s",   catalog_1);
    if(flag_match_subgroups){
      n_particles_group_1= (int *)ADaPS_fetch(plist_1->data,"n_particles_subgroup_%s",    catalog_1);
      group_offset_1     = (int *)ADaPS_fetch(plist_1->data,"particle_offset_subgroup_%s",catalog_1);
      rank_offset_1      =((int *)ADaPS_fetch(plist_1->data,"rank_offset_subgroup_%s",    catalog_1))[0];
    }
    else{
      n_particles_group_1= (int *)ADaPS_fetch(plist_1->data,"n_particles_group_%s",    catalog_1);
      group_offset_1     = (int *)ADaPS_fetch(plist_1->data,"particle_offset_group_%s",catalog_1);
      rank_offset_1      =((int *)ADaPS_fetch(plist_1->data,"rank_offset_group_%s",    catalog_1))[0];
    }
if(SID.I_am_Master) fprintf(stderr,"test5\n"); SID_Barrier(SID.COMM_WORLD);

    // Fetch needed info for catalog_2
    n_particles_2_all  =((size_t *)ADaPS_fetch(plist_2->data,"n_particles_all_%s",catalog_2))[0];
    n_particles_2_local=((size_t *)ADaPS_fetch(plist_2->data,"n_particles_%s",    catalog_2))[0];
    id_2_local         = (size_t *)ADaPS_fetch(plist_2->data,"particle_ids_%s",   catalog_2);
    if(flag_match_subgroups){
      n_particles_group_2_local= (int *)ADaPS_fetch(plist_2->data,"n_particles_subgroup_%s",    catalog_2);
      group_offset_2_local     = (int *)ADaPS_fetch(plist_2->data,"particle_offset_subgroup_%s",catalog_2);
      rank_offset_2_local      =((int *)ADaPS_fetch(plist_2->data,"rank_offset_subgroup_%s",    catalog_2))[0];
    }
    else{
      n_particles_group_2_local= (int *)ADaPS_fetch(plist_2->data,"n_particles_group_%s",    catalog_2);
      group_offset_2_local     = (int *)ADaPS_fetch(plist_2->data,"particle_offset_group_%s",catalog_2);
      rank_offset_2_local      =((int *)ADaPS_fetch(plist_2->data,"rank_offset_group_%s",    catalog_2))[0];
    }
if(SID.I_am_Master) fprintf(stderr,"test6\n"); SID_Barrier(SID.COMM_WORLD);

    // Sort the mark list (if there is one)
    if(mark_list_1==NULL)
      mark_list_index_1=NULL;
    else
      merge_sort((void *)mark_list_1,(size_t)n_mark_1_local,&mark_list_index_1,SID_INT,SORT_COMPUTE_INDEX,FALSE);
if(SID.I_am_Master) fprintf(stderr,"test7\n"); SID_Barrier(SID.COMM_WORLD);

    // We need group ids of the 1st catalog if matching substructure
    if(flag_match_substructure){
      group_index_1=(int *)SID_malloc(sizeof(int)*n_particles_1);
      for(i_particle=0;i_particle<n_particles_1;i_particle++)
        group_index_1[i_particle]=-1;
      if(mark_list_1==NULL){
        for(i_group=0;i_group<n_groups_1;i_group++){
          for(i_particle=0;i_particle<n_particles_group_1[i_group];i_particle++)
            group_index_1[group_offset_1[i_group]+i_particle]=i_group;
        }
      }
      else{
        for(i_group=0,i_mark=0;i_group<n_groups_1 && i_mark<n_mark_1_local;i_group++){
          if(mark_list_1[mark_list_index_1[i_mark]]==i_group){
            for(i_particle=0;i_particle<n_particles_group_1[i_group];i_particle++)
              group_index_1[group_offset_1[i_group]+i_particle]=i_group;
            i_mark++;
          }
        }
      }
    }
if(SID.I_am_Master) fprintf(stderr,"test8\n"); SID_Barrier(SID.COMM_WORLD);

    // Generate group ids for the second catalog
    group_index_2_local=(int *)SID_malloc(sizeof(int)*n_particles_2_local);
    for(i_particle=0;i_particle<n_particles_2_local;i_particle++)
      group_index_2_local[i_particle]=-1;
    if(mark_list_2_local==NULL){
      mark_list_index_2_local=NULL;
      for(i_group=0;i_group<n_groups_2_local;i_group++){
        for(i_particle=0;i_particle<n_particles_group_2_local[i_group];i_particle++)
          group_index_2_local[group_offset_2_local[i_group]+i_particle]=i_group;
      }
    }
    else{
      merge_sort((void *)mark_list_2_local,(size_t)n_mark_2_local,&mark_list_index_2_local,SID_INT,SORT_COMPUTE_INDEX,FALSE);
      for(i_group=0,i_mark=0;i_group<n_groups_2_local && i_mark<n_mark_2_local;i_group++){
        if(mark_list_2_local[mark_list_index_2_local[i_mark]]==i_group){
          for(i_particle=0;i_particle<n_particles_group_2_local[i_group];i_particle++)
            group_index_2_local[group_offset_2_local[i_group]+i_particle]=i_group;
          i_mark++;
        }
      }
    }
if(SID.I_am_Master) fprintf(stderr,"test9\n"); SID_Barrier(SID.COMM_WORLD);

    // Create the arrays that will hold the results
    match     =(int   *)SID_malloc(sizeof(int)*n_match_1);
    if(SID.n_proc>1){
      match_rank=(int *)SID_malloc(sizeof(int)*n_match_1);
      for(i_mark=0;i_mark<n_match_1;i_mark++)
        match_rank[i_mark]=SID.My_rank;
    }
    else
      match_rank=NULL;
    hist_score=(float *)SID_calloc(sizeof(float)*hist_size);
    hist_list =(int   *)SID_calloc(sizeof(int)  *hist_size);
if(SID.I_am_Master) fprintf(stderr,"test10\n"); SID_Barrier(SID.COMM_WORLD);

    // Initialize the matches to -1, the default if no match is found
    //   (for substructure matching, -1 means a group is matched only to itself)
    for(i_mark=0;i_mark<n_match_1;i_mark++)
      match[i_mark]=-1;

    // Check if we are storing the match score.  Allocate array and set flag=TRUE if so.
    if(check_mode_for_flag(mode,MATCH_STORE_SCORE) || SID.n_proc>1){
      match_score=(float *)SID_calloc(sizeof(float)*n_match_1);
      flag_store_score=TRUE;
    }
    else
      flag_store_score=FALSE;
  
    // Perform matching
if(SID.I_am_Master) fprintf(stderr,"test11\n"); SID_Barrier(SID.COMM_WORLD);
    for(i_rank=0,n_match=0;i_rank<SID.n_proc;i_rank++){

       // Point arrays at themselves when matching a rank against itself
       if(i_rank==0){
         n_particles_2 =n_particles_2_local;
         group_index_2 =group_index_2_local;
         group_offset_2=group_offset_2_local;
         rank_offset_2 =rank_offset_2_local;
         id_2          =id_2_local;
       }
       // Create buffer arrays and perform exchanges if matching against another rank
       else{
         // Determine how many particles and groups need to be echanged
if(SID.I_am_Master) fprintf(stderr,"test12\n"); SID_Barrier(SID.COMM_WORLD);
         exchange_ring_buffer(&n_particles_2_local,
                              sizeof(size_t),
                              1,
                              &n_particles_2,
                              NULL,
                              i_rank);
if(SID.I_am_Master) fprintf(stderr,"test13\n"); SID_Barrier(SID.COMM_WORLD);
         exchange_ring_buffer(&n_groups_2_local,
                              sizeof(size_t),
                              1,
                              &n_groups_2,
                              NULL,
                              i_rank);
if(SID.I_am_Master) fprintf(stderr,"test14\n"); SID_Barrier(SID.COMM_WORLD);

         // Allocate arrays for the exchange
         group_index_2 =(int    *)SID_calloc(n_particles_2*sizeof(int));
         group_offset_2=(int    *)SID_calloc(n_groups_2   *sizeof(int));
         id_2          =(size_t *)SID_calloc(n_particles_2*sizeof(size_t));

         // Perform exchange
         exchange_ring_buffer(group_index_2_local,
                              sizeof(int),
                              n_particles_2_local,
                              group_index_2,
                              &n_particles_2,
                              i_rank);
         n_groups_2_local_s=(size_t)n_groups_2_local;
         exchange_ring_buffer(group_offset_2_local,
                              sizeof(int),
                              n_groups_2_local_s,
                              group_offset_2,
                              &n_groups_2_s,
                              i_rank);
         exchange_ring_buffer(&rank_offset_2_local,
                              sizeof(int),
                              1,
                              &rank_offset_2,
                              NULL,
                              i_rank);
         exchange_ring_buffer(id_2_local,
                              sizeof(size_t),
                              n_particles_2_local,
                              id_2,
                              &n_particles_2,
                              i_rank);
if(SID.I_am_Master) fprintf(stderr,"test15\n"); SID_Barrier(SID.COMM_WORLD);
       }

       // Sort the target (ie second) list of ids
if(SID.I_am_Master) fprintf(stderr,"test16\n"); SID_Barrier(SID.COMM_WORLD);
       merge_sort((void *)id_2,(size_t)n_particles_2,&index_2,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
if(SID.I_am_Master) fprintf(stderr,"test17\n"); SID_Barrier(SID.COMM_WORLD);

       for(i_group=0,i_particle=0,i_mark=0;
           i_group<n_groups_1 && i_mark<n_match_1;
           i_group++){
         // Decide if we are trying to match this halo or not
         if(flag_read_marked){
           if(mark_list_index_1!=NULL){
             if(mark_list_1[mark_list_index_1[i_mark]]==i_group){
               flag_read_this=TRUE;
             }
             else{
               flag_read_this=FALSE;
             }
           }
           else{
             if(mark_list_1[i_mark]==i_group)
               flag_read_this=TRUE;
             else
               flag_read_this=FALSE;          
           }
         }
         else
           flag_read_this=TRUE;

         // If we are trying to match this halo, then ...
         if(flag_read_this){
           // Loop over particles in i_group'th group in catalog_1 -- generate match scores to 
           //   overlaping halos in catalog_2 from the rank that we are currently analyzing
           n_p_match=n_particles_group_1[i_group];
           merge_sort((void *)(&(id_1[group_offset_1[i_group]])),(size_t)(n_p_match),&index_1,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
           // ... initialize the histogram ...
           for(i_hist=0;i_hist<hist_size;i_hist++){
             hist_score[i_hist]=0;
             hist_list[i_hist] =0;
           }
           // ... perform loop over particles  ...
           for(i_particle=group_offset_1[i_group],k_particle=0,n_hist=0;
               k_particle<n_p_match;
               i_particle++,k_particle++){
             // Check if this particle matches the current catalog_1 halo ...
             id_1_i    =id_1[group_offset_1[i_group]+index_1[k_particle]];
             j_particle=find_index(id_2,id_1_i,n_particles_2,index_2);
             while(id_1_i>id_2[index_2[j_particle]] && j_particle<(n_particles_2-2)) j_particle++;
             if(id_1_i>id_2[index_2[j_particle]])                                    j_particle++;
             if(id_1_i==id_2[index_2[j_particle]])
               flag_match_continue=TRUE;
             else
               flag_match_continue=FALSE;
             // ... if it does ...
             while(flag_match_continue){
               // ...and if the corresponding group_id in catalog_2 has already been 
               //    involved in a match, then add to its histogram bin 
               for(i_hist=0,flag=TRUE;i_hist<n_hist && flag;i_hist++){
                 if(hist_list[i_hist]==group_index_2[index_2[j_particle]]){
                   switch(flag_match_substructure){
                   case TRUE:
                     rank               =1.;
                     hist_score[i_hist]+=1.;
                     break;
                   default:
                     rank               =(float)(1+index_2[j_particle]-group_offset_2[hist_list[i_hist]]);
                     hist_score[i_hist]+=(float)pow(rank,-TWO_THIRDS);
                     break;
                   }
                   flag=FALSE;
                 }
               }
               // ... else, add the new group id to the histogram (if it is a valid halo; ie. +ve and included in the matching process)
               if(flag && group_index_2[index_2[j_particle]]>=0){
                 // Reallocate the arrays if we have too many matches ...
                 if(n_hist>=hist_size){
                   hist_size*=2;
                   hist_list =(int   *)realloc(hist_list, hist_size*sizeof(int));
                   hist_score=(float *)realloc(hist_score,hist_size*sizeof(float));
                   for(i_hist=n_hist;i_hist<hist_size;i_hist++){
                      hist_score[i_hist]=0.;
                      hist_list[i_hist] =0;
                   }
                 }
                 switch(flag_match_substructure){
                 case TRUE:
                   if(n_particles_group_1[i_group]<
                      n_particles_group_2[group_index_2[index_2[j_particle]]]){
                     hist_list[n_hist] =group_index_2[index_2[j_particle]];
                     rank              =1.;
                     hist_score[n_hist]=1.;
                     n_hist++;
                   }
                   break;
                 default:
                   hist_list[n_hist] =group_index_2[index_2[j_particle]];
                   rank              =(float)(1+index_2[j_particle]-group_offset_2[hist_list[n_hist]]);
                   hist_score[n_hist]=(float)pow(rank,-TWO_THIRDS);
                   n_hist++;
                   break;
                 }
               }
               if(flag_match_substructure && j_particle<(n_particles_2-1)){
                 j_particle++;
                 if(id_1_i!=id_2[index_2[j_particle]]) flag_match_continue=FALSE;
               }
               else
                 flag_match_continue=FALSE;
             }
           }
if(SID.I_am_Master) fprintf(stderr,"test18\n"); SID_Barrier(SID.COMM_WORLD);
           SID_free(SID_FARG index_1);
           // Find the best match ON THE CURRENT RANK; assign a value of -1 for unmatched groups
           j_hist=0;
           if(n_hist>0){
             switch(flag_match_substructure){
               // Keep the smallest system if we are matching substructure
             case TRUE:
               for(i_hist=1,j_hist=0;i_hist<n_hist;i_hist++){
                 if(n_particles_group_2[hist_list[i_hist]]<
                    n_particles_group_2[hist_list[j_hist]])
                   j_hist=i_hist;
               }
               break;
               // Keep the system with the highest score otherwise
             default:
               for(i_hist=1,j_hist=0;i_hist<n_hist;i_hist++){
                 if(hist_score[i_hist]>hist_score[j_hist])
                   j_hist=i_hist;
               }
               break;
             }
             // Set match results here
             if(match_score!=NULL){ // Always true if n_proc>1
               if(hist_score[j_hist]>match_score[i_mark]){
                 match_score[i_mark]=hist_score[j_hist];
                 match[i_mark]      =hist_list[j_hist]+rank_offset_2;
                 if(match_rank!=NULL) // Always true if n_proc>1
                   match_rank[i_mark]=i_rank;
               }
             }
             else
                match[i_mark]=hist_list[j_hist];
if(SID.I_am_Master) fprintf(stderr,"test19\n"); SID_Barrier(SID.COMM_WORLD);
             n_match++;
           }
           i_mark++;
         }
       } // Loop over local groups
if(SID.I_am_Master) fprintf(stderr,"test20\n"); SID_Barrier(SID.COMM_WORLD);

       // Clean-up
       SID_free((void **)&index_2);
       if(i_rank!=0){
          SID_free(SID_FARG group_index_2);
          SID_free(SID_FARG group_offset_2);
          SID_free(SID_FARG id_2);
       }

    } // Loop over ranks
    calc_sum_global(&n_match,&n_match_all,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
    n_match=n_match_all;

    // Clean-up
    SID_free((void **)&group_index_2);
    SID_free((void **)&hist_list);
    SID_free((void **)&hist_score);
    if(mark_list_index_1!=NULL)
      SID_free((void **)&mark_list_index_1);
    if(mark_list_index_2!=NULL)
      SID_free((void **)&mark_list_index_2);
    if(!check_mode_for_flag(mode,MATCH_STORE_SCORE) && SID.n_proc>1)
      SID_free(SID_FARG match_score);

    if(flag_match_subgroups)
      SID_log("%d of %d subgroups matched to %d subgroups...",SID_LOG_CONTINUE,
              n_match,
              n_match_1_all,
              n_groups_2_all);
    else
      SID_log("%d of %d groups matched to %d groups...",SID_LOG_CONTINUE,
              n_match,
              n_match_1_all,
              n_groups_2_all);

    // Store matches
    if(check_mode_for_flag(mode,MATCH_STORE_2))
      plist_store=plist_2;
    else
      plist_store=plist_1;
    if(check_mode_for_flag(mode,MATCH_BACK)){
      ADaPS_store(&(plist_store->data),(void *)(&n_match),"n_back_match_%s",ADaPS_SCALAR_INT,catalog_1to2);
      ADaPS_store(&(plist_store->data),(void *)(match),   "back_match_%s",  ADaPS_DEFAULT,   catalog_1to2);
      if(flag_store_score)
        ADaPS_store(&(plist_store->data),(void *)(match_score),"back_match_score_%s",ADaPS_DEFAULT,catalog_1to2);
      if(SID.n_proc>1)
        ADaPS_store(&(plist_store->data),(void *)(match_rank),"back_match_rank_%s",ADaPS_DEFAULT,catalog_1to2);
    }
    else{
      ADaPS_store(&(plist_store->data),(void *)(&n_match),"n_match_%s",ADaPS_SCALAR_INT,catalog_1to2);
      ADaPS_store(&(plist_store->data),(void *)(match),   "match_%s",  ADaPS_DEFAULT,   catalog_1to2);
      if(flag_store_score)
        ADaPS_store(&(plist_store->data),(void *)(match_score),"match_score_%s",ADaPS_DEFAULT,catalog_1to2);
      if(SID.n_proc>1)
        ADaPS_store(&(plist_store->data),(void *)(match_rank),"match_rank_%s",ADaPS_DEFAULT,catalog_1to2);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
  else{
    if(check_mode_for_flag(mode,MATCH_STORE_2))
      plist_store=plist_2;
    else
      plist_store=plist_1;
    n_match=0;
    if(check_mode_for_flag(mode,MATCH_BACK)){
      ADaPS_store(&(plist_store->data),(void *)(&n_match),"n_back_match_%s",ADaPS_SCALAR_INT,catalog_1to2);
    }
    else{
      ADaPS_store(&(plist_store->data),(void *)(&n_match),"n_match_%s",ADaPS_SCALAR_INT,catalog_1to2);
    }
    SID_log("NO GROUPS TO MATCH!",SID_LOG_CLOSE);
  }
  //SID_profile_stop(SID_PROFILE_DEFAULT);
  //SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
}

