#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpHalos.h>

void match_halos(plist_info  *plist_1_in,
                 int          i_file_1_in,
                 int         *mark_list_1,
                 int          n_mark_1_local,
                 plist_info  *plist_2_in,
                 int          i_file_2_in,
                 int         *mark_list_2_local,
                 int          n_mark_2_local,
                 const char  *catalog_1to2,
                 int          mode){
  plist_info  *plist_store;
  plist_info  *plist_1;
  int          i_file_1;
  plist_info  *plist_2;
  int          i_file_2;
  int     i_group;
  int     j_group;
  int     i_mark;
  size_t  i_particle;
  size_t  j_particle;
  size_t  n_particles_1;
  size_t  n_particles_1_local;
  size_t  n_particles_1_all;
  size_t  n_particles_2;
  size_t  n_particles_2_local;
  size_t  n_particles_2_all;
  size_t  n_groups_1_dim;
  size_t  n_mark_1_dim;
  int     n_mark_1_all;
  int     n_groups_1;
  int     n_groups_1_local;
  int     n_groups_1_all;
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
  int    *file_index_1;
  size_t *group_offset_1;
  size_t *group_offset_2;
  size_t *group_offset_2_local;
  int    *file_index_2;
  int    *file_index_2_local;
  size_t *mark_list_index_1      =NULL;
  size_t *mark_list_index_2      =NULL;
  size_t *mark_list_index_2_local=NULL;
  size_t *index_2      =NULL;
  size_t *id_1;
  size_t *id_2;
  size_t *id_2_local;
  size_t  id_1_i;
  int    *group_index_1;
  int    *group_index_2;
  int    *group_index_2_local;
  int    *particle_rank_1;
  int    *hist_list  =NULL;
  float  *hist_score =NULL;
  float  *match_score=NULL;
  float   rank;
  int     i_rank;
  int     i_lookup;
  int     flag;
  size_t  group;
  int     counter;
  double  f_min=0.9;
  int     flag_store_score;
  int     flag_read_marked;
  int     flag_match_this;
  char    catalog_1[5];
  char    catalog_2[5];
  int     flag_match_subgroups;
  int     flag_match_substructure;
  int     flag_PHK_decomp;
  size_t  n_particles_exchange_1;
  int     n_groups_exchange_1;
  size_t  n_particles_exchange_2;
  int     n_groups_exchange_2;
  int     n_score_lookup_table;
  float  *score_lookup_table;

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

  // matching of substructure is broken right now.
  if(flag_match_substructure)
    SID_trap_error("Matching of substructure is broken in match_halos() right now.",ERROR_LOGIC);

  // Check if we are using Peano-Hilbert Key (PHK) decomposition
  int n_groups_boundary_1;
  int n_particles_boundary_1;
  int n_groups_boundary_2_local;
  int n_particles_boundary_2_local;
  int n_bits_PHK_1;
  int n_bits_PHK_2;
  if(ADaPS_exist(plist_1->data,"n_bits_PHK_%s",catalog_1)){
     flag_PHK_decomp=TRUE;
     if(!ADaPS_exist(plist_2->data,"n_bits_PHK_%s",catalog_2))
        SID_trap_error("Both catalogs must be loaded with PHK decompositions.",ERROR_LOGIC);
     n_bits_PHK_1=((int *)ADaPS_fetch(plist_1->data,"n_bits_PHK_%s",catalog_1))[0];
     n_bits_PHK_2=((int *)ADaPS_fetch(plist_2->data,"n_bits_PHK_%s",catalog_2))[0];
     n_particles_boundary_1      =((int *)ADaPS_fetch(plist_1->data,"n_particles_boundary_%s",catalog_1))[0];
     n_particles_boundary_2_local=((int *)ADaPS_fetch(plist_2->data,"n_particles_boundary_%s",catalog_2))[0];
     if(flag_match_subgroups){
        n_groups_boundary_1      =((int *)ADaPS_fetch(plist_1->data,"n_subgroups_boundary_%s",catalog_1))[0];
        n_groups_boundary_2_local=((int *)ADaPS_fetch(plist_2->data,"n_subgroups_boundary_%s",catalog_2))[0];
     }
     else{
        n_groups_boundary_1      =((int *)ADaPS_fetch(plist_1->data,"n_groups_boundary_%s",   catalog_1))[0];
        n_groups_boundary_2_local=((int *)ADaPS_fetch(plist_2->data,"n_groups_boundary_%s",   catalog_2))[0];
     }
     if(n_bits_PHK_1!=n_bits_PHK_2 && (n_bits_PHK_1!=0 && n_bits_PHK_2!=0)) // n_bits=0 if there are no groups
        SID_trap_error("PHK decomposition of the two catalogs is incompatible (ie. %d!=%d).",ERROR_LOGIC,n_bits_PHK_1,n_bits_PHK_2);
  }
  else
     flag_PHK_decomp     =FALSE;

  // Fetch some info and print log message
  if(flag_match_subgroups){
    SID_log("Matching subgroups from catalog {%s} to catalog {%s}...",
            SID_LOG_OPEN|SID_LOG_TIMER,
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
              SID_LOG_OPEN|SID_LOG_TIMER,
              catalog_1,
              catalog_2);
    else
      SID_log("Matching groups catalog {%s} to catalog {%s}...",
              SID_LOG_OPEN|SID_LOG_TIMER,
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
  n_groups_1_local=n_groups_1;

  // A status message if using PHKs
  if(flag_PHK_decomp)
     SID_log("Groups are PH decomposed using %d-bit per dimension keys.",SID_LOG_COMMENT,n_bits_PHK_1);
  else
     SID_log("Groups are NOT PH decomposed.",SID_LOG_COMMENT);

  // Allow a mark list to be passed if we are
  //   only interested in certain objects
  calc_sum_global(&n_mark_1_local,&n_mark_1_all,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
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
  if(n_groups_1_all>0 && n_groups_2_all>0){                

    // Create a lookup table to save time computing match scores
    n_score_lookup_table=256;
    score_lookup_table  =(float *)SID_malloc(sizeof(float)*n_score_lookup_table);
    for(i_lookup=0;i_lookup<n_score_lookup_table;i_lookup++)
      score_lookup_table[i_lookup]=(float)pow((double)i_lookup,-TWO_THIRDS);

    // Fetch needed info for catalog_1
    n_particles_1_all  =((size_t *)ADaPS_fetch(plist_1->data,"n_particles_all_%s",catalog_1))[0];
    n_particles_1_local=((size_t *)ADaPS_fetch(plist_1->data,"n_particles_%s",    catalog_1))[0];
    n_particles_1      =((size_t *)ADaPS_fetch(plist_1->data,"n_particles_%s",    catalog_1))[0];
    id_1               = (size_t *)ADaPS_fetch(plist_1->data,"particle_ids_%s",   catalog_1);
    if(flag_match_subgroups){
      n_particles_group_1= (int    *)ADaPS_fetch(plist_1->data,"n_particles_subgroup_%s",    catalog_1);
      group_offset_1     = (size_t *)ADaPS_fetch(plist_1->data,"particle_offset_subgroup_%s",catalog_1);
      file_index_1       = (int    *)ADaPS_fetch(plist_1->data,"file_index_subgroups_%s",    catalog_1);
    }
    else{
      n_particles_group_1= (int    *)ADaPS_fetch(plist_1->data,"n_particles_group_%s",    catalog_1);
      group_offset_1     = (size_t *)ADaPS_fetch(plist_1->data,"particle_offset_group_%s",catalog_1);
      file_index_1       = (int    *)ADaPS_fetch(plist_1->data,"file_index_groups_%s",    catalog_1);
    }
//size_t count;
//size_t id_max=(((size_t)2160)*((size_t)2160)*((size_t)2160))-1;
//fprintf(stderr,"[%d]cat1: {%s}\n",SID.My_rank,catalog_1); 
//if(SID.I_am_Master)
//   ADaPS_status(plist_1->data);
//SID_Barrier(SID.COMM_WORLD);
//for(size_t i_id=0,count=0;i_id<n_particles_1_local;i_id++){
//   if(id_1[i_id]>id_max) count++;
//}
//fprintf(stderr,"[%d]count1=%zd\n",SID.My_rank,count);

    // Fetch needed info for catalog_2
    n_particles_2_all  =((size_t *)ADaPS_fetch(plist_2->data,"n_particles_all_%s",catalog_2))[0];
    n_particles_2_local=((size_t *)ADaPS_fetch(plist_2->data,"n_particles_%s",    catalog_2))[0];
    id_2_local         = (size_t *)ADaPS_fetch(plist_2->data,"particle_ids_%s",   catalog_2);
    if(flag_match_subgroups){
      n_particles_group_2_local= (int    *)ADaPS_fetch(plist_2->data,"n_particles_subgroup_%s",    catalog_2);
      group_offset_2_local     = (size_t *)ADaPS_fetch(plist_2->data,"particle_offset_subgroup_%s",catalog_2);
      file_index_2_local       = (int    *)ADaPS_fetch(plist_2->data,"file_index_subgroups_%s",    catalog_2);
    }
    else{
      n_particles_group_2_local= (int    *)ADaPS_fetch(plist_2->data,"n_particles_group_%s",    catalog_2);
      group_offset_2_local     = (size_t *)ADaPS_fetch(plist_2->data,"particle_offset_group_%s",catalog_2);
      file_index_2_local       = (int    *)ADaPS_fetch(plist_2->data,"file_index_groups_%s",    catalog_2);
    }
fprintf(stderr,"cat2: {%s}\n",catalog_2); 
if(SID.I_am_Master)
   ADaPS_status(plist_2->data);
SID_Barrier(SID.COMM_WORLD);
for(size_t i_id=0,count=0;i_id<n_particles_2_local;i_id++){
   if(id_2_local[i_id]>id_max) count++;
}   
fprintf(stderr,"[%d]count2=%zd\n",SID.My_rank,count);
SID_exit(ERROR_NONE);

    // Set the number of particles that need to be checked by exchanges
    if(flag_PHK_decomp){
       n_groups_exchange_1   =n_groups_boundary_1;
       n_particles_exchange_1=n_particles_boundary_1;
       n_groups_exchange_2   =n_groups_boundary_2_local;
       n_particles_exchange_2=n_particles_boundary_2_local;
    }
    else{
       n_groups_exchange_1   =n_groups_1;
       n_particles_exchange_1=n_particles_1;
       n_groups_exchange_2   =n_groups_2_local;
       n_particles_exchange_2=n_particles_2_local;
    }

    // Sort the mark list (if there is one)
    if(mark_list_1==NULL)
      mark_list_index_1=NULL;
    else
      merge_sort((void *)mark_list_1,(size_t)n_mark_1_local,&mark_list_index_1,SID_INT,SORT_COMPUTE_INDEX,FALSE);

    // We need group ids of the 1st catalog if matching substructure
    group_index_1  =(int *)SID_malloc(sizeof(int)*n_particles_1);
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

    // We need the rank of each particle
    if(!flag_match_substructure){
       particle_rank_1=(int *)SID_malloc(sizeof(int)*n_particles_1);
       for(i_particle=0;i_particle<n_particles_1;i_particle++)
          particle_rank_1[i_particle]=-1;
       for(i_group=0;i_group<n_groups_1;i_group++){
          for(i_particle=0;i_particle<n_particles_group_1[i_group];i_particle++)
             particle_rank_1[group_offset_1[i_group]+i_particle]=i_particle+1;
       }
    }

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
      SID_free(SID_FARG mark_list_index_2_local);
    }

    // Create the array that will hold the results and set values to a 
    //   default of -1 for groups that don't get matched
    //   (for substructure matching, -1 means a group is matched only to itself)
    match =(int *)SID_malloc(sizeof(int)*n_match_1);
    for(i_mark=0;i_mark<n_match_1;i_mark++)
      match[i_mark]=-1;

    // Check if we are storing the match score.  Allocate array and set flag=TRUE if so.
    match_score=(float *)SID_calloc(sizeof(float)*n_match_1);
    if(check_mode_for_flag(mode,MATCH_STORE_SCORE))
      flag_store_score=TRUE;
    else
      flag_store_score=FALSE;

    // Allocate some buffers for rank exchanges
    size_t  n_particles_2_max;
    int     n_groups_2_max;
    SID_Allreduce(&n_particles_2_local,&n_particles_2_max,1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
    SID_Allreduce(&n_groups_2_local,   &n_groups_2_max,   1,SID_INT,   SID_MAX,SID.COMM_WORLD);
    if(SID.n_proc>1){
       SID_log("Allocate exchange buffers...",SID_LOG_OPEN);
       SID_log("Group    buffer size=%d groups",   SID_LOG_COMMENT,n_groups_2_max);
       SID_log("Particle buffer size=%d particles",SID_LOG_COMMENT,n_particles_2_max);
       n_particles_group_2=(int    *)SID_calloc(n_groups_2_max   *sizeof(int));
       group_index_2      =(int    *)SID_calloc(n_particles_2_max*sizeof(int));
       file_index_2       =(int    *)SID_calloc(n_groups_2_max   *sizeof(int));
       id_2               =(size_t *)SID_calloc(n_particles_2_max*sizeof(size_t));
       index_2            =(size_t *)SID_calloc(n_particles_2_max*sizeof(size_t));
       SID_log("Done.",SID_LOG_CLOSE|SID_LOG_TIMER);
    }

    // Create array of histograms for the matching
    short int  *n_hist_array;
    short int  *hist_size_array;
    int       **hist_list_array;
    float     **hist_score_array;
    SID_log("Allocate histogram arrays...",SID_LOG_OPEN|SID_LOG_TIMER);
    n_hist_array    =(short int  *)SID_calloc(sizeof(short int)*n_groups_1);
    hist_size_array =(short int  *)SID_calloc(sizeof(short int)*n_groups_1);
    hist_list_array =(int       **)SID_malloc(sizeof(int *)    *n_groups_1);
    hist_score_array=(float     **)SID_malloc(sizeof(float *)  *n_groups_1);
    for(i_group=0;i_group<n_groups_1;i_group++){ 
       // Decide if we are trying to match this halo or not
       if(flag_read_marked){
          if(mark_list_index_1!=NULL){
             if(mark_list_1[mark_list_index_1[i_mark]]==i_group){
               flag_match_this=TRUE;
             }
             else{
               flag_match_this=FALSE;
             }
          }
          else{
             if(mark_list_1[i_mark]==i_group)
               flag_match_this=TRUE;
             else
               flag_match_this=FALSE;          
          }
       }
       else
          flag_match_this=TRUE;
       n_hist_array[i_group]=0;
       if(flag_match_this){
          hist_size_array[i_group] =4;
          hist_list_array[i_group] =(int   *)SID_calloc(sizeof(int)  *hist_size_array[i_group]);
          hist_score_array[i_group]=(float *)SID_calloc(sizeof(float)*hist_size_array[i_group]);
       }
       else{
          hist_size_array[i_group] =0;
          hist_list_array[i_group] =NULL;
          hist_score_array[i_group]=NULL;
       }
    }
    int hist_size_max=0;
    for(i_group=0;i_group<n_groups_1;i_group++){
       if(hist_size_array[i_group]>hist_size_max) 
          hist_size_max=hist_size_array[i_group]; 
    }
    SID_log("Done.",SID_LOG_CLOSE|SID_LOG_TIMER);

    // Perform matching
    for(i_rank=0,n_match=0;i_rank<SID.n_proc;i_rank++){
       if(SID.n_proc>1)
          SID_log("Processing rank %4d of %4d...",SID_LOG_OPEN|SID_LOG_TIMER,i_rank+1,SID.n_proc);

       // Sort the IDs of the first catalog.  Only needs to be done for i_rank==0 for
       //   all particles and then for i_rank==1 for exchanged particles.
       size_t *index_1      =NULL;
       size_t *index_2_local=NULL;
       if(i_rank==0){
          SID_log("Sorting all IDs...",SID_LOG_OPEN|SID_LOG_TIMER);
          merge_sort(id_1,      (size_t)(n_particles_1),      &index_1,      SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE); //*
          merge_sort(id_2_local,(size_t)(n_particles_2_local),&index_2_local,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
          SID_log("Done.",SID_LOG_CLOSE);
       }
       // We just need to sort the boundary particles once
       else if(i_rank==1){
          SID_log("Sorting boundary IDs...",SID_LOG_OPEN|SID_LOG_TIMER);
          SID_free(SID_FARG index_1);
          SID_free(SID_FARG index_2_local);
          merge_sort(id_1,      (size_t)(n_particles_exchange_1),&index_1,      SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
          merge_sort(id_2_local,(size_t)(n_particles_exchange_2),&index_2_local,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
          SID_log("Done.",SID_LOG_CLOSE);
       }

       // Point arrays at themselves when matching a rank against itself
       if(i_rank==0){
         n_particles_2=n_particles_2_local;
         n_groups_2   =n_groups_2_local;
         if(SID.n_proc>1){
            memcpy(n_particles_group_2,n_particles_group_2_local,sizeof(int)   *n_groups_2);
            memcpy(group_index_2,      group_index_2_local,      sizeof(int)   *n_particles_2);
            memcpy(file_index_2,       file_index_2_local,       sizeof(int)   *n_groups_2);
            memcpy(id_2,               id_2_local,               sizeof(size_t)*n_particles_2);
            memcpy(index_2,            index_2_local,            sizeof(size_t)*n_particles_2);
         }
         else{
            n_particles_group_2=n_particles_group_2_local;
            group_index_2      =group_index_2_local;
            file_index_2       =file_index_2_local;
            id_2               =id_2_local;
            index_2            =index_2_local;
         }
       }
       // Create buffer arrays and perform exchanges if matching against another rank
       else{
         // During exchanges, we only need to work with the n_exchange halos/particles in each catalog
         n_particles_1      =n_particles_exchange_1;
         n_groups_1         =n_groups_exchange_1;
         n_particles_2_local=n_particles_exchange_2;
         n_groups_2_local   =n_groups_exchange_2;

         // Determine how many particles and groups need to be exchanged
         SID_log("Performing exchange...",SID_LOG_OPEN|SID_LOG_TIMER);
         SID_log("particle count...",SID_LOG_COMMENT);
         exchange_ring_buffer(&n_particles_2_local,
                              sizeof(size_t),
                              1,
                              &n_particles_2,
                              NULL,
                              i_rank);
         SID_log("group count...",SID_LOG_COMMENT);SID_Barrier(SID.COMM_WORLD);
         exchange_ring_buffer(&n_groups_2_local, // send
                              sizeof(int),       // send/recv
                              1,                 // send
                              &n_groups_2,       // recv
                              NULL,              // recv
                              i_rank); 

         // Perform exchange
         SID_log("group sizes...",SID_LOG_COMMENT);
         exchange_ring_buffer(n_particles_group_2_local,
                              sizeof(int),
                              (size_t)n_groups_2_local,
                              n_particles_group_2,
                              NULL,
                              i_rank);
         SID_log("group indices...",SID_LOG_COMMENT);
         exchange_ring_buffer(group_index_2_local,
                              sizeof(int),
                              (size_t)n_particles_2_local,
                              group_index_2,
                              NULL,
                              i_rank);
         SID_log("file indices...",SID_LOG_COMMENT);
         exchange_ring_buffer(file_index_2_local,
                              sizeof(int),
                              (size_t)n_groups_2_local,
                              file_index_2,
                              NULL,
                              i_rank);
         SID_log("ids...",SID_LOG_COMMENT);
         exchange_ring_buffer(id_2_local,
                              sizeof(size_t),
                              (size_t)n_particles_2_local,
                              id_2,
                              NULL,
                              i_rank);
         SID_log("sort indices...",SID_LOG_COMMENT);
         exchange_ring_buffer(index_2_local,
                              sizeof(size_t),
                              (size_t)n_particles_2_local,
                              index_2,
                              NULL,
                              i_rank);
         SID_log("Done.",SID_LOG_CLOSE);
       }

       // Perform matching
       int       max_id_gap     = 100; // If two consecuative IDs are more than this appart, use bisection to find next index
       int       flag_use_bisect=TRUE; // This will skip all the zeros at the beginning of the list due to unused particles and jump big id gaps
       size_t    id_1_last      =   0;
       short int hist_size;
       size_t    idx_1;
       size_t    idx_2;
       SID_log("Performing matching...",SID_LOG_OPEN);
       for(i_particle=0,j_particle=0;
           i_particle<n_particles_1 && j_particle<n_particles_2;
           i_particle++,flag_use_bisect=FALSE){
          idx_1  =index_1[i_particle];
          i_group=group_index_1[idx_1];

          // What group does this particle belong to and are we trying to match it?
          if(i_group>=0)
            hist_size=hist_size_array[i_group];
          else
            hist_size=0;
          if(hist_size>0){

             // Set the ID we are looking for and decide how we're going to find it
             id_1_i=id_1[idx_1];
             if((id_1_i-id_1_last)>max_id_gap)
                flag_use_bisect=TRUE;
             id_1_last=id_1_i;

             // Look for the current particle in the second catalog ...
             if(flag_use_bisect)
                j_particle=find_index(id_2,id_1_i,n_particles_2,index_2);
             if(j_particle<(n_particles_2-1)){
                while(id_1_i>id_2[index_2[j_particle]] && j_particle<(n_particles_2-2)) j_particle++;
                if(id_1_i>id_2[index_2[j_particle]])                                    j_particle++;
             }
             idx_2  =index_2[j_particle];
             j_group=group_index_2[idx_2];

             // ... if we found it, update the appropriate histogram ...
             if(id_1_i==id_2[idx_2] && j_group>=0){
                int file_index_2_j;
                hist_list     =hist_list_array[i_group];
                hist_score    =hist_score_array[i_group];
                file_index_2_j=file_index_2[j_group];
                // ... if the corresponding group_id in catalog_2 has already been 
                //     involved in a match, then add to its score
                short int i_hist;
                short int j_hist;
                switch(flag_match_substructure){
                   case TRUE:
                      for(i_hist=0,flag=TRUE;i_hist<n_hist_array[i_group] && flag;i_hist++){
                         if(hist_list[i_hist]==file_index_2_j) 
                            hist_score[i_hist]+=1.;
                         flag=FALSE;
                      }
                      break;
                   default:
                      for(i_hist=0,flag=TRUE;i_hist<n_hist_array[i_group] && flag;i_hist++){
                         if(hist_list[i_hist]==file_index_2_j){
                            switch(particle_rank_1[idx_1]<n_score_lookup_table){
                              case TRUE:
                                 hist_score[i_hist]+=score_lookup_table[particle_rank_1[idx_1]];
                                 break;
                              default:
                                 hist_score[i_hist]+=(float)pow((double)(particle_rank_1[idx_1]),-TWO_THIRDS);
                                 break;
                            }
                            flag=FALSE;
                         }
                      }
                      break;
                }

                // ... else, add the new group id to the histogram (if it is a valid halo; ie. +ve and included in the matching process)
                if(flag){
                   // Reallocate the arrays if we have too many matches ...
                   if(n_hist_array[i_group]>=hist_size_array[i_group]){
                      short int n_old;
                      n_old=hist_size_array[i_group];
                      if(hist_size_array[i_group]>=hist_size_max)
                         hist_size_max=2*(int)hist_size_array[i_group];
                      hist_size_array[i_group]*=2;
                      /*
                      hist_list_array[i_group] =(int   *)SID_realloc(hist_list_array[i_group], (size_t)(hist_size_array[i_group])*sizeof(int));
                      hist_score_array[i_group]=(float *)SID_realloc(hist_score_array[i_group],(size_t)(hist_size_array[i_group])*sizeof(float));
                      hist_list                =hist_list_array[i_group];
                      hist_score               =hist_score_array[i_group];
                      */
                      hist_list_array[i_group] =(int   *)SID_calloc((size_t)(hist_size_array[i_group])*sizeof(int));
                      hist_score_array[i_group]=(float *)SID_calloc((size_t)(hist_size_array[i_group])*sizeof(float));
                      memcpy(hist_list_array[i_group], hist_list, (size_t)n_old*sizeof(int));
                      memcpy(hist_score_array[i_group],hist_score,(size_t)n_old*sizeof(float));
                      SID_free(SID_FARG hist_list);
                      SID_free(SID_FARG hist_score);
                      hist_list =hist_list_array[i_group];
                      hist_score=hist_score_array[i_group];
                   }
                   switch(flag_match_substructure){
                      case TRUE:
                         if(n_particles_group_1[i_group]<
                            n_particles_group_2[j_group]){
                            hist_list[0]         =file_index_2_j;
                            hist_score[0]        =1.;
                            n_hist_array[i_group]=1;
                         }
                         break;
                      default:
                         hist_list[n_hist_array[i_group]] =file_index_2_j;
                         switch(particle_rank_1[idx_1]<n_score_lookup_table){
                           case TRUE:
                              hist_score[n_hist_array[i_group]]=score_lookup_table[particle_rank_1[idx_1]];
                              break;
                           default:
                              hist_score[n_hist_array[i_group]]=(float)pow((double)(particle_rank_1[idx_1]),-TWO_THIRDS);
                              break;
                         }
                         n_hist_array[i_group]++;
                         break;
                   }
                } // Add a new match to the histogram for this group
             } // If we found this particle and it's involved in the matching
          } // If we are matching this group and particle
       } // Loop over local particles in catalog 1
       SID_log("Done.",SID_LOG_CLOSE);
       SID_free(SID_FARG index_1);
       SID_free(SID_FARG index_2_local);
       if(SID.n_proc>1)
          SID_log("Done.",SID_LOG_CLOSE);
    } // Loop over ranks

    // Determine the best matches
    for(i_group=0,i_mark=0;i_group<n_groups_1_local;i_group++){
       if(hist_size_array[i_group]>0){
          // Find the best match here
          if(n_hist_array[i_group]>0){
             short int i_hist;
             short int j_hist;
             j_hist    =0;
             hist_list =hist_list_array[i_group];
             hist_score=hist_score_array[i_group];
             switch(flag_match_substructure){
                // Keep the smallest system if we are matching substructure
                case TRUE:
                   for(i_hist=1;i_hist<n_hist_array[i_group];i_hist++){
                     if(n_particles_group_2[hist_list[i_hist]]<
                        n_particles_group_2[hist_list[j_hist]])
                       j_hist=i_hist;
                   }
                   break;
                // Keep the system with the highest score otherwise
                default:
                   for(i_hist=1;i_hist<n_hist_array[i_group];i_hist++){
                     if(hist_score[i_hist]>hist_score[j_hist])
                       j_hist=i_hist;
                   }
                   break;
             }
             // Set match results here; the default (set above) is -1 for unmatched groups
             if(hist_score[j_hist]>0.){
                match_score[i_mark]=hist_score[j_hist];
                match[i_mark]      =hist_list[j_hist];
                n_match++;
             }
             if(match[i_mark]<(-1) || match[i_mark]>=n_groups_2_all){
                SID_log_warning("Invalid match_id (%d) for i_mark/i_group=%d/%d.  There are %d objects in the target catalog. n_hist/j_hist=%d/%d",
                                ERROR_LOGIC,match[i_mark],i_mark,i_group,n_groups_2_all,n_hist_array[i_group],j_hist);
                SID_trap_error("Invalid match_id (%d) for i_mark/i_group=%d/%d.  There are %d objects in the target catalog.",
                               ERROR_LOGIC,match[i_mark],i_mark,i_group,n_groups_2_all);
             }
          }
          i_mark++;
       }
    }

    // Check that the hist_size array didn't over-flow
    int hist_size_max_global;
    int hist_size_limit=32767; // Largest number supported by 16-bit signed int
    calc_max_global(&hist_size_max,&hist_size_max_global,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
    SID_log("Largest hist size=%d",SID_LOG_COMMENT,hist_size_max_global);
    if(hist_size_max_global>hist_size_limit)
       SID_trap_error("The histogram size array has overflowed (ie. %d>%d)",SID_LOG_COMMENT,hist_size_max_global,hist_size_limit);

    // How many matches across all ranks?
    calc_sum_global(&n_match,&n_match_all,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
    n_match=n_match_all;

    // Clean-up
    SID_free(SID_FARG score_lookup_table);
    SID_free(SID_FARG group_index_1);
    if(!flag_match_substructure)
       SID_free(SID_FARG particle_rank_1);
    for(i_group=0;i_group<n_groups_1_local;i_group++){
       if(hist_size_array[i_group]>0){
          SID_free(SID_FARG hist_list_array[i_group]);
          SID_free(SID_FARG hist_score_array[i_group]);
       }
    }
    SID_free(SID_FARG n_hist_array);
    SID_free(SID_FARG hist_size_array);
    SID_free(SID_FARG hist_list_array);
    SID_free(SID_FARG hist_score_array);
    if(SID.n_proc>1){
       SID_free(SID_FARG n_particles_group_2);
       SID_free(SID_FARG group_index_2);
       SID_free(SID_FARG file_index_2);
       SID_free(SID_FARG id_2);
       SID_free(SID_FARG index_2);
    }
    SID_free(SID_FARG group_index_2_local);
    if(mark_list_index_1!=NULL)
      SID_free(SID_FARG mark_list_index_1);
    if(mark_list_index_2!=NULL)
      SID_free(SID_FARG mark_list_index_2);
    if(!flag_store_score && match_score!=NULL)
      SID_free(SID_FARG match_score);

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

    SID_log("%d of %d matched to %d.",SID_LOG_COMMENT,
            n_match,
            n_match_1_all,
            n_groups_2_all);
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
    SID_log("NO GROUPS TO MATCH!",SID_LOG_COMMENT);
    SID_log("Done.",SID_LOG_CLOSE);
  }
}

