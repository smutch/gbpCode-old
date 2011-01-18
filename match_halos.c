#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpHalos.h>

void match_halos(plist_info  *plist_1_in,
                 int          i_file_1_in,
                 int         *mark_list_1,
                 int          n_mark_1,
                 plist_info  *plist_2_in,
                 int          i_file_2_in,
                 int         *mark_list_2,
                 int          n_mark_2,
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
  size_t  n_particles_2;
  int     n_groups_1;
  size_t  n_groups_1_dim;
  size_t  n_mark_1_dim;
  int     n_groups_2;
  int     n_match;
  int     n_match_1;
  char    parm_txt[256];
  int    *match;
  int    *n_particles_group_1;
  int    *n_particles_group_2;
  int    *group_offset_1;
  int    *group_offset_2;
  size_t *mark_list_index_1=NULL;
  size_t *mark_list_index_2=NULL;
  size_t *index_1;
  size_t *index_2;
  size_t *id_1;
  size_t *id_2;
  size_t  id_1_i;
  int    *group_index_1;
  int    *group_index_2;
  size_t  n_hist;
  int    *hist_list;
  float  *hist_score;
  float  *match_score;
  float   rank;
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

  SID_profile_start("match_groups",SID_PROFILE_NOTMPIENABLED);

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
    if(ADaPS_exist(plist_2->data,"n_subgroups_%s",catalog_2))
      n_groups_2=((int *)ADaPS_fetch(plist_2->data,"n_subgroups_%s",catalog_2))[0];
    else
      n_groups_2=0;
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
    if(ADaPS_exist(plist_2->data,"n_groups_%s",catalog_2))
      n_groups_2=((int *)ADaPS_fetch(plist_2->data,"n_groups_%s",catalog_2))[0];
    else
      n_groups_2=0;
  }
  
  // Allow a mark list to be passed if we are
  //   only interested in certain objects
  if(mark_list_1!=NULL && n_mark_1>0){
    flag_read_marked=TRUE;
    n_match_1       =n_mark_1;
  }
  else{
    flag_read_marked=FALSE;
    n_match_1       =n_groups_1;
  }

  // Perform matching if there are groups to match
  if(n_groups_1>0 && n_groups_2>0){                

    // Fetch needed info for catalog_1
    n_particles_1=((size_t *)ADaPS_fetch(plist_1->data,"n_particles_%s",catalog_1))[0];
    id_1         = (size_t *)ADaPS_fetch(plist_1->data,"particle_ids_%s",catalog_1);
    if(flag_match_subgroups){
      n_particles_group_1=(int *)ADaPS_fetch(plist_1->data,"n_particles_subgroup_%s",    catalog_1);
      group_offset_1     =(int *)ADaPS_fetch(plist_1->data,"particle_offset_subgroup_%s",catalog_1);
    }
    else{
      n_particles_group_1=(int *)ADaPS_fetch(plist_1->data,"n_particles_group_%s",    catalog_1);
      group_offset_1     =(int *)ADaPS_fetch(plist_1->data,"particle_offset_group_%s",catalog_1);
    }

    // Fetch needed info for catalog_2
    n_particles_2=((size_t *)ADaPS_fetch(plist_2->data,"n_particles_%s",catalog_2))[0];
    id_2         = (size_t *)ADaPS_fetch(plist_2->data,"particle_ids_%s",catalog_2);
    if(flag_match_subgroups){
      n_particles_group_2=(int *)ADaPS_fetch(plist_2->data,"n_particles_subgroup_%s",    catalog_2);
      group_offset_2     =(int *)ADaPS_fetch(plist_2->data,"particle_offset_subgroup_%s",catalog_2);
    }
    else{
      n_particles_group_2=(int *)ADaPS_fetch(plist_2->data,"n_particles_group_%s",    catalog_2);
      group_offset_2     =(int *)ADaPS_fetch(plist_2->data,"particle_offset_group_%s",catalog_2);
    }

    // Initialize a bunch of arrays
    //   (We need group ids of 1st catalog if matching substructure)
    if(mark_list_1==NULL)
      mark_list_index_1=NULL;
    else
      merge_sort((void *)mark_list_1,(size_t)n_mark_1,&mark_list_index_1,SID_INT,SORT_COMPUTE_INDEX,FALSE);
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
        for(i_group=0,i_mark=0;i_group<n_groups_1 && i_mark<n_mark_1;i_group++){
          if(mark_list_1[mark_list_index_1[i_mark]]==i_group){
            for(i_particle=0;i_particle<n_particles_group_1[i_group];i_particle++)
              group_index_1[group_offset_1[i_group]+i_particle]=i_group;
            i_mark++;
          }
        }
      }
    }
    group_index_2=(int *)SID_malloc(sizeof(int)*n_particles_2);
    for(i_particle=0;i_particle<n_particles_2;i_particle++)
      group_index_2[i_particle]=-1;
    if(mark_list_2==NULL){
      mark_list_index_2=NULL;
      for(i_group=0;i_group<n_groups_2;i_group++){
        for(i_particle=0;i_particle<n_particles_group_2[i_group];i_particle++)
          group_index_2[group_offset_2[i_group]+i_particle]=i_group;
      }
    }
    else{
      merge_sort((void *)mark_list_2,(size_t)n_mark_2,&mark_list_index_2,SID_INT,SORT_COMPUTE_INDEX,FALSE);
      for(i_group=0,i_mark=0;i_group<n_groups_2 && i_mark<n_mark_2;i_group++){
        if(mark_list_2[mark_list_index_2[i_mark]]==i_group){
          for(i_particle=0;i_particle<n_particles_group_2[i_group];i_particle++)
            group_index_2[group_offset_2[i_group]+i_particle]=i_group;
          i_mark++;
        }
      }
    }
    merge_sort((void *)id_2,(size_t)n_particles_2,&index_2,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
    match     =(int   *)SID_malloc(sizeof(int)  *n_match_1);
    hist_score=(float *)SID_malloc(sizeof(float)*hist_size);
    hist_list =(int   *)SID_malloc(sizeof(int)  *hist_size);
    for(i_mark=0;i_mark<n_match_1;i_mark++) 
      match[i_mark]=-1;

    // Check if we are storing the match score.  Allocate array and set flag=TRUE if so.
    if(check_mode_for_flag(mode,MATCH_STORE_SCORE)){
      match_score=(float *)SID_malloc(sizeof(float)*n_match_1);
      for(i_mark=0;i_mark<n_match_1;i_mark++) 
        match_score[i_mark]=0.;
      flag_store_score=TRUE;
    }
    else
      flag_store_score=FALSE;
  
    // Perform matching
    for(i_group=0,i_particle=0,i_mark=0,n_match=0;
        i_group<n_groups_1 && i_mark<n_match_1;
        i_group++){
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
      if(flag_read_this){
        // Loop over particles in i_group'th group in catalog_1 -- generate match scores to overlaping halos in catalog_2
        n_p_match=n_particles_group_1[i_group];
        merge_sort((void *)(&(id_1[group_offset_1[i_group]])),(size_t)(n_p_match),&index_1,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
        for(i_hist=0;i_hist<hist_size;i_hist++){
          hist_score[i_hist]=0.;
          hist_list[i_hist] =0;
        }
        for(i_particle=group_offset_1[i_group],k_particle=0,n_hist=0;
            k_particle<n_p_match;
            i_particle++,k_particle++){
          id_1_i    =id_1[group_offset_1[i_group]+index_1[k_particle]];
          j_particle=find_index(id_2,id_1_i,n_particles_2,index_2);
          while(id_1_i>id_2[index_2[j_particle]] && j_particle<(n_particles_2-2)) j_particle++;
          if(id_1_i>id_2[index_2[j_particle]])                                    j_particle++;
          // If we have matched particles ...
          if(id_1_i==id_2[index_2[j_particle]])
            flag_match_continue=TRUE;
          else
            flag_match_continue=FALSE;
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
              if(n_hist>=hist_size){
                hist_size*=2;
                hist_list =(int   *)realloc(hist_list, hist_size*sizeof(int));
                hist_score=(float *)realloc(hist_score,hist_size*sizeof(float));
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
        SID_free((void **)&index_1);
        // Find the best match; assign a value of -1 for unmatched groups
        //   (for substructure matching, -1 means a group is matched only to itself)
        j_hist=0;
        switch(n_hist){
        case 0:
          switch(flag_match_substructure){
          case TRUE:
            //match[i_mark]=i_mark;
            match[i_mark]=-1;
            break;
          default:
            match[i_mark]=-1;
            break;
          }
          break;
        default:
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
          match[i_mark]=hist_list[j_hist];
          if(flag_store_score)
            match_score[i_mark]=hist_score[j_hist];
          n_match++;
          break;
        }
        i_mark++;
      }
    }

    // Clean-up
    SID_free((void **)&group_index_2);
    SID_free((void **)&index_2);
    SID_free((void **)&hist_list);
    SID_free((void **)&hist_score);
    if(mark_list_index_1!=NULL)
      SID_free((void **)&mark_list_index_1);
    if(mark_list_index_2!=NULL)
      SID_free((void **)&mark_list_index_2);

    if(flag_match_subgroups)
      SID_log("%d of %d subgroups matched to %d subgroups...",SID_LOG_CONTINUE,
              n_match,
              n_match_1,
              n_groups_2);
    else
      SID_log("%d of %d groups matched to %d groups...",SID_LOG_CONTINUE,
              n_match,
              n_match_1,
              n_groups_2);

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
    }
    else{
      ADaPS_store(&(plist_store->data),(void *)(&n_match),"n_match_%s",ADaPS_SCALAR_INT,catalog_1to2);
      ADaPS_store(&(plist_store->data),(void *)(match),   "match_%s",  ADaPS_DEFAULT,   catalog_1to2);
      if(flag_store_score)
        ADaPS_store(&(plist_store->data),(void *)(match_score),"match_score_%s",ADaPS_DEFAULT,catalog_1to2);
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
  SID_profile_stop(SID_PROFILE_DEFAULT);
}
