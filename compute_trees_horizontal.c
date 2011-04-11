#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void compute_trees_horizontal(char *filename_halo_root_in,
                              char *filename_cat_root_in,
                              char *filename_root_out,
                              int   i_read_start,
                              int   i_read_stop,
                              int   i_read_step,
                              int   n_search,
                              int  *flag_clean){
  plist_info  plist1;
  plist_info  plist2;
  plist_info  plist3;
  char        filename_cat1_1[256];
  char        filename_cat1_2[256];
  char        filename_cat1_3[256];
  char        filename_cat2_1[256];
  char        filename_cat2_2[256];
  char        filename_matches_out[256];
  char        filename_groups[256];
  char        filename_subgroups[256];
  char        filename_subgroup_properties_in[256];
  char        filename_group_properties_in[256];
  char        filename_subgroup_properties_out[256];
  char        filename_group_properties_out[256];
  char        text_temp[256];
  char        group_text_prefix[5];
  FILE       *fp;
  FILE       *fp_matches_out;
  SID_fp      fp_group_properties_out;
  SID_fp      fp_subgroup_properties_out;
  FILE       *fp_group_properties_in;
  FILE       *fp_subgroup_properties_in;
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
  int         k_write;
  int         i_file_1;
  int         i_file_2;
  int         i_file_3;
  int         j_file_1;
  int         j_file_2;
  int         i_read;
  int         j_read;
  int         i_read_1;
  int         i_read_2;
  int         i_read_3;
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
  int         my_descendant;
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
  int        *drop_list_1=NULL;
  int        *drop_list_2=NULL;
  int        *bridge_list_1=NULL;
  int        *bridge_list_2=NULL;
  int        *drop_index=NULL;
  int        *match_id_bridge=NULL;
  int         i_find,n_find;
  int         flag_continue;
  int         flag_drop;
  int        *match_id=NULL;
  int        *back_match_id=NULL;
  int         n_progenitors_max;
  int         i_search;
  int         flag_dropped;
  int         flag_first;
  int         n_particles;
  int         n_particles_max;
  int         trunk_index;
  int         biggest_stray;
  int         biggest_stray_drop;
  int         biggest_stray_bridge;
  int        *n_particles_1=NULL;
  int        *n_particles_2=NULL;
  int        *n_groups=NULL;
  int        *n_subgroups=NULL;
  int       **progenitor_id=NULL;
  int       **progenitor_id_group=NULL;
  int       **progenitor_id_subgroup=NULL;
  int       **descendant_id=NULL;
  int       **descendant_id_group=NULL;
  int       **descendant_id_subgroup=NULL;
  int       **file_offset=NULL;
  int       **file_offset_group=NULL;
  int       **file_offset_subgroup=NULL;
  short int **n_progenitors_subgroup=NULL;
  short int **n_progenitors_group=NULL;
  short int **n_progenitors=NULL;
  int       **tree_id=NULL;
  int       **tree_id_group=NULL;
  int       **tree_id_subgroup=NULL;
  int         max_tree_id_group;
  int         max_tree_id_subgroup;
  int         max_tree_id;
  int       **n_subgroups_group=NULL;
  int        *n_subgroups_group_1=NULL;
  size_t    **sort_id=NULL;
  size_t    **sort_group_id=NULL;
  size_t    **sort_subgroup_id=NULL;
  size_t     *match_index=NULL;
  size_t     *back_match_index=NULL;
  float      *match_score=NULL;
  int         flag_match_subgroups;
  char        cat_label_1[20];
  char        cat_label_2[20];
  char        cat_name_1[20];
  char        cat_name_2[20];
  char        filename_log[256];
  int         flag_keep_strays=FALSE;
  halo_info   properties;
  int         n_k_match=2;

  sprintf(filename_log,"%s.trees_horizontal_log",filename_root_out);

  SID_log("Constructing vertical merger trees for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_start,i_read_stop,i_read_step);

  if(n_search<1)
    SID_trap_error("n_search=%d but must be at least 1",ERROR_LOGIC);
    
  // Count the maximum number of substructures in one snapshot; needed to set array sizes
  SID_log("Counting indices...",SID_LOG_OPEN);
  if(SID.I_am_Master){
  for(i_read=i_read_stop,n_progenitors_max=0;i_read>=i_read_start;i_read-=i_read_step){

    sprintf(filename_groups,"%s_%03d.catalog_groups",filename_halo_root_in,i_read);
    fp=fopen(filename_groups,"r");
    fread(&n_groups_1,sizeof(int),1,fp);
    fclose(fp);

    sprintf(filename_subgroups,"%s_%03d.catalog_subgroups",filename_halo_root_in,i_read);
    fp=fopen(filename_subgroups,"r");
    fread(&n_subgroups_1,sizeof(int),1,fp);
    fclose(fp);

    n_progenitors_max=MAX(n_progenitors_max,n_subgroups_1);
    n_progenitors_max=MAX(n_progenitors_max,n_groups_1);
  }
  }
  SID_Bcast(&n_progenitors_max,sizeof(int),1,SID.COMM_WORLD);
  n_search+=2; // Need indices for current and last i_file as well
  SID_log("Done. (max_index=%d)",SID_LOG_CLOSE,n_progenitors_max);

  // Initialize arrays
  SID_log("Initializing arrays...",SID_LOG_OPEN);
  tree_id_subgroup       =(int       **)SID_malloc(sizeof(int       *)*n_search);
  tree_id_group          =(int       **)SID_malloc(sizeof(int       *)*n_search);
  n_progenitors_subgroup =(short int **)SID_malloc(sizeof(short int *)*n_search);
  n_progenitors_group    =(short int **)SID_malloc(sizeof(short int *)*n_search);
  progenitor_id_subgroup =(int       **)SID_malloc(sizeof(int       *)*n_search);
  progenitor_id_group    =(int       **)SID_malloc(sizeof(int       *)*n_search);
  descendant_id_subgroup =(int       **)SID_malloc(sizeof(int       *)*n_search);
  descendant_id_group    =(int       **)SID_malloc(sizeof(int       *)*n_search);
  file_offset_subgroup   =(int       **)SID_malloc(sizeof(int       *)*n_search);
  file_offset_group      =(int       **)SID_malloc(sizeof(int       *)*n_search);
  n_subgroups_group      =(int       **)SID_malloc(sizeof(int       *)*n_search);
  sort_subgroup_id       =(size_t    **)SID_malloc(sizeof(size_t    *)*n_search);
  sort_group_id          =(size_t    **)SID_malloc(sizeof(size_t    *)*n_search);
  drop_list_1            =(int        *)SID_malloc(sizeof(int)*n_progenitors_max);
  drop_list_2            =(int        *)SID_malloc(sizeof(int)*n_progenitors_max);
  n_groups               =(int        *)SID_malloc(sizeof(int)*n_search);
  n_subgroups            =(int        *)SID_malloc(sizeof(int)*n_search);
  for(i_search=0;i_search<n_search;i_search++){
    tree_id_subgroup[i_search]       =(int       *)SID_malloc(sizeof(int)      *n_progenitors_max);
    tree_id_group[i_search]          =(int       *)SID_malloc(sizeof(int)      *n_progenitors_max);                
    n_progenitors_subgroup[i_search] =(short int *)SID_malloc(sizeof(short int)*n_progenitors_max);                
    n_progenitors_group[i_search]    =(short int *)SID_malloc(sizeof(short int)*n_progenitors_max);                
    progenitor_id_subgroup[i_search] =(int       *)SID_malloc(sizeof(int)      *n_progenitors_max);                
    progenitor_id_group[i_search]    =(int       *)SID_malloc(sizeof(int)      *n_progenitors_max);                
    descendant_id_subgroup[i_search] =(int       *)SID_malloc(sizeof(int)      *n_progenitors_max);                
    descendant_id_group[i_search]    =(int       *)SID_malloc(sizeof(int)      *n_progenitors_max);                
    file_offset_subgroup[i_search]   =(int       *)SID_malloc(sizeof(int)      *n_progenitors_max);                
    file_offset_group[i_search]      =(int       *)SID_malloc(sizeof(int)      *n_progenitors_max);                
    n_subgroups_group[i_search]      =(int       *)SID_malloc(sizeof(int)      *n_progenitors_max);                
    sort_subgroup_id[i_search]       =(size_t    *)SID_malloc(sizeof(size_t)   *n_progenitors_max);
    sort_group_id[i_search]          =(size_t    *)SID_malloc(sizeof(size_t)   *n_progenitors_max);
    // Initialize to a unique negative number (leave -1 unused)
    for(i_group=0;i_group<n_progenitors_max;i_group++){
      progenitor_id_subgroup[i_search][i_group]=-i_group-2;
      progenitor_id_group[i_search][i_group]   =-i_group-2;
    }
  }
  init_plist(&plist1,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  init_plist(&plist2,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  init_plist(&plist3,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  SID_log("Done.",SID_LOG_CLOSE);

  // Process the first file separately
  //   (just give everything ids from a running index) ...
  SID_log("Processing snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_stop);
  i_file_1=i_read_stop;
  i_file_2=i_read_stop+1;
  i_file_3=i_read_stop+2;
  i_read_1=i_read_stop;
  i_read_2=i_read_stop+1*i_read_step;
  i_read_3=i_read_stop+2*i_read_step;
  sprintf(filename_cat1_1,"%03d",i_read_1);
  sprintf(filename_cat1_2,"%03d",i_read_2);
  
  // ... read groups
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
  read_groups(filename_halo_root_in,i_read_1,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES,&plist1,filename_cat1_1);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

  // ... process subgroups ...
  n_subgroups[i_file_3%n_search]=((int *)ADaPS_fetch(plist1.data,"n_subgroups_%s",filename_cat1_1))[0];
  n_subgroups[i_file_2%n_search]=n_subgroups[i_file_3%n_search];
  n_subgroups[i_file_1%n_search]=n_subgroups[i_file_3%n_search];
  SID_log("Assigning %d subgroup IDs...",SID_LOG_OPEN,n_subgroups[i_file_1%n_search]);
  for(i_group=0,max_id_subgroup=0,max_tree_id_subgroup=0;i_group<n_subgroups[i_file_3%n_search];i_group++,max_id_subgroup++,max_tree_id_subgroup++){
    n_progenitors_subgroup[i_file_3%n_search][i_group]=1; // Default is 1:1 merger history for first snapshot
    n_progenitors_subgroup[i_file_2%n_search][i_group]=1; // Default is 1:1 merger history for first snapshot
    n_progenitors_subgroup[i_file_1%n_search][i_group]=0; // Computed later
    progenitor_id_subgroup[i_file_3%n_search][i_group]=max_id_subgroup;
    progenitor_id_subgroup[i_file_2%n_search][i_group]=max_id_subgroup;
    progenitor_id_subgroup[i_file_1%n_search][i_group]=max_id_subgroup;
    descendant_id_subgroup[i_file_3%n_search][i_group]=max_id_subgroup;
    descendant_id_subgroup[i_file_2%n_search][i_group]=max_id_subgroup;
    descendant_id_subgroup[i_file_1%n_search][i_group]=max_id_subgroup;
    tree_id_subgroup[i_file_3%n_search][i_group]      =max_tree_id_subgroup;
    tree_id_subgroup[i_file_2%n_search][i_group]      =max_tree_id_subgroup;
    tree_id_subgroup[i_file_1%n_search][i_group]      =max_tree_id_subgroup;
    file_offset_subgroup[i_file_3%n_search][i_group]  =0;  // These carry the number of files that
    file_offset_subgroup[i_file_2%n_search][i_group]  =0;  //   the progenitor is offset from
    file_offset_subgroup[i_file_1%n_search][i_group]  =0;  //   each descendant
  }
  SID_free((void **)&sort_subgroup_id[i_file_3%n_search]);

  // Create sorting indices for the progenitor IDs
  merge_sort(progenitor_id_subgroup[i_file_3%n_search],
             n_subgroups[i_file_3%n_search],
             &(sort_subgroup_id[i_file_3%n_search]),
             SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  SID_free((void **)&sort_subgroup_id[i_file_2%n_search]);
  merge_sort(progenitor_id_subgroup[i_file_2%n_search],
             n_subgroups[i_file_2%n_search],
             &(sort_subgroup_id[i_file_2%n_search]),
             SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  SID_free((void **)&sort_subgroup_id[i_file_1%n_search]);
  merge_sort(progenitor_id_subgroup[i_file_1%n_search],
             n_subgroups[i_file_1%n_search],
             &(sort_subgroup_id[i_file_1%n_search]),
             SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  SID_log("Done.",SID_LOG_CLOSE);

  // ... process groups ...
  n_groups[i_file_3%n_search]=((int *)ADaPS_fetch(plist1.data,"n_groups_%s",filename_cat1_1))[0];
  n_groups[i_file_2%n_search]=n_groups[i_file_3%n_search];
  n_groups[i_file_1%n_search]=n_groups[i_file_3%n_search];
  n_subgroups_group_1        = (int *)ADaPS_fetch(plist1.data,"n_subgroups_group_%s",filename_cat1_1);
  SID_log("Assigning %d group IDs...",SID_LOG_OPEN,n_groups[i_file_1%n_search]);
  for(i_group=0,max_id_group=0,max_tree_id_group=0;i_group<n_groups[i_file_1%n_search];i_group++,max_id_group++,max_tree_id_group++){
    n_progenitors_group[i_file_3%n_search][i_group] = 1; // Default is 1:1 merger history
    n_progenitors_group[i_file_2%n_search][i_group] = 1; // Default is 1:1 merger history
    n_progenitors_group[i_file_1%n_search][i_group] = 0; // Computed later
    progenitor_id_group[i_file_3%n_search][i_group] =max_id_group;
    progenitor_id_group[i_file_2%n_search][i_group] =max_id_group;
    progenitor_id_group[i_file_1%n_search][i_group] =max_id_group;
    descendant_id_group[i_file_3%n_search][i_group] =max_id_group;
    descendant_id_group[i_file_2%n_search][i_group] =max_id_group;
    descendant_id_group[i_file_1%n_search][i_group] =max_id_group;
    tree_id_group[i_file_3%n_search][i_group]       =max_tree_id_group;
    tree_id_group[i_file_2%n_search][i_group]       =max_tree_id_group;
    tree_id_group[i_file_1%n_search][i_group]       =max_tree_id_group;
    n_subgroups_group[i_file_3%n_search][i_group]   =n_subgroups_group_1[i_group];
    n_subgroups_group[i_file_2%n_search][i_group]   =n_subgroups_group_1[i_group];
    n_subgroups_group[i_file_1%n_search][i_group]   =n_subgroups_group_1[i_group];
    file_offset_group[i_file_3%n_search][i_group]   =0; 
    file_offset_group[i_file_2%n_search][i_group]   =0; 
    file_offset_group[i_file_1%n_search][i_group]   =0; 
  }
  SID_free((void **)&sort_group_id[i_file_3%n_search]);
  merge_sort(progenitor_id_group[i_file_3%n_search],
             n_groups[i_file_3%n_search],
             &(sort_group_id[i_file_3%n_search]),
             SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  SID_free((void **)&sort_group_id[i_file_2%n_search]);
  merge_sort(progenitor_id_group[i_file_2%n_search],
             n_groups[i_file_2%n_search],
             &(sort_group_id[i_file_2%n_search]),
             SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  SID_free((void **)&sort_group_id[i_file_1%n_search]);
  merge_sort(progenitor_id_group[i_file_1%n_search],
             n_groups[i_file_1%n_search],
             &(sort_group_id[i_file_1%n_search]),
             SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_log("Done.",SID_LOG_CLOSE);

  // The first snapshot is done now (set to defaults) ... now loop over all other snapshots
  for(i_read=i_read_1-i_read_step, // Check this index
        i_file=i_file_1-1,         // Check this index
        j_file=0,                  // Check this index
        i_write=i_read_stop,       // Check this index
        j_write=0,                 // Check this index
        k_write=i_read_stop;       // Check this index
      i_read>=i_read_start;        // Check this index
      i_read-=i_read_step,         // Check this index
        i_file--,                  // Check this index
        j_file++){                 // Check this index
    SID_log("Processing snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read);

    // Shift the snapshot info up the list
    i_file_3=i_file_2;
    i_file_2=i_file_1;
    i_file_1=i_file;
    i_read_3=i_read_2;
    i_read_2=i_read_1;
    i_read_1=i_read;
    sprintf(filename_cat1_1,"%03d",i_read_1);
    sprintf(filename_cat1_2,"%03d",i_read_2);
    sprintf(filename_cat1_3,"%03d",i_read_3);
    ADaPS_free(&(plist2.data));
    plist2.data=plist1.data;
    plist1.data=NULL;

    // Read the group info
    SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
    read_groups(filename_halo_root_in,i_read_1,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES,&plist1,filename_cat1_1);
    SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

    // Retrieve the number of halos (copy to another array because we need to hold onto it for a while)
    n_subgroups_group_1           = (int *)ADaPS_fetch(plist1.data,"n_subgroups_group_%s",filename_cat1_1);
    n_subgroups[i_file_1%n_search]=((int *)ADaPS_fetch(plist1.data,"n_subgroups_%s",filename_cat1_1))[0];
    n_groups[i_file_1%n_search]   =((int *)ADaPS_fetch(plist1.data,"n_groups_%s",filename_cat1_1))[0];
    for(i_group=0;i_group<n_groups[i_file_1%n_search];i_group++)
      n_subgroups_group[i_file_1%n_search][i_group]=n_subgroups_group_1[i_group];

    // Loop twice (1st to process subgroups, 2nd to process groups)
    for(k_match=0;k_match<n_k_match;k_match++){

      // Initialize counters
      n_match_halos       =0;
      n_drop              =0;
      n_mergers_drop      =0;
      n_drop_found        =0;
      n_strays_drop       =0;
      n_multimatch        =0;
      n_bridge_candidates =0;
      n_bridge_systems    =0;
      n_found_bridge      =0;
      n_bridges           =0;
      n_mergers_bridge    =0;
      n_strays_bridge     =0;
      n_mergers           =0;
      n_found             =0;
      n_strays            =0;
      n_sputter           =0;
      biggest_stray       =0;
      biggest_stray_drop  =0;
      biggest_stray_bridge=0;

      // Set working array pointers to point to group or subgroup information (alternately, depending on k_match)
      switch(k_match){
      case 0:
        flag_match_subgroups=MATCH_SUBGROUPS;
        n_groups_1    =n_subgroups[i_file_1%n_search];
        n_groups_2    =n_subgroups[i_file_2%n_search];
        n_groups_3    =n_subgroups[i_file_3%n_search];
        n_particles_1 =(int *)ADaPS_fetch(plist1.data,"n_particles_subgroup_%s",filename_cat1_1);
        n_particles_2 =(int *)ADaPS_fetch(plist2.data,"n_particles_subgroup_%s",filename_cat1_2);
        n_progenitors =n_progenitors_subgroup;
        progenitor_id =progenitor_id_subgroup;
        descendant_id =descendant_id_subgroup;
        file_offset   =file_offset_subgroup;
        sort_id       =sort_subgroup_id;
        tree_id       =tree_id_subgroup;
        max_id        =max_id_subgroup;
        max_tree_id   =max_tree_id_subgroup;
        sprintf(group_text_prefix,"sub");
        break;
      case 1:
        flag_match_subgroups=MATCH_GROUPS;
        n_groups_1    =n_groups[i_file_1%n_search];
        n_groups_2    =n_groups[i_file_2%n_search];
        n_groups_3    =n_groups[i_file_3%n_search];
        n_particles_1 =(int *)ADaPS_fetch(plist1.data,"n_particles_group_%s",filename_cat1_1);
        n_particles_2 =(int *)ADaPS_fetch(plist2.data,"n_particles_group_%s",filename_cat1_2);
        n_progenitors =n_progenitors_group;
        progenitor_id =progenitor_id_group;
        descendant_id =descendant_id_group;
        file_offset   =file_offset_group;
        sort_id       =sort_group_id;
        tree_id       =tree_id_group;
        max_id        =max_id_group;
        max_tree_id   =max_tree_id_group;
        sprintf(group_text_prefix,"");
        break;
      }

      SID_log("Processing %d %sgroups of snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,n_groups_1,group_text_prefix,i_read_1);
      n_unprocessed=n_groups_1;
      if(n_groups_1>0 && n_groups_2>0 && n_groups_3>0){
        sprintf(cat_label_1,"1to2_%sgroup",group_text_prefix);
        sprintf(cat_label_2,"3to2_%sgroup",group_text_prefix);
        sprintf(cat_name_1, "match_1to2_%sgroup",group_text_prefix);
        sprintf(cat_name_2, "back_match_3to2_%sgroup",group_text_prefix);

        // Perform forward-matching
        SID_log("Performing matching...",SID_LOG_OPEN|SID_LOG_TIMER);
        match_halos(&plist1,i_read_1,NULL,0,&plist2,i_read_2,NULL,0,cat_label_1,flag_match_subgroups|MATCH_STORE_SCORE);
        match_id=(int  *)ADaPS_fetch(plist1.data,cat_name_1); // These are the ids in file_2 that match to each object in file_1 ...
        merge_sort((void *)match_id,(size_t)n_groups_1,&match_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE); // ... and these are their sort indices

        // Perform back-matching; needed when separating bridges from mergers in file_2
        match_halos(&plist1,i_read_1,NULL,0,&plist2,i_read_2,NULL,0,cat_label_2,flag_match_subgroups|MATCH_BACK|MATCH_STORE_2);
        if(j_file>0){
          back_match_id=(int  *)ADaPS_fetch(plist2.data,cat_name_2); // These are the ids in file_3 that back-match to each object in file_2 ...
          merge_sort((void *)back_match_id,(size_t)n_groups_3,&back_match_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE); // ... and these are their sort indices
        }

        // Initialize tree pointer-arrays with unique and negative dummy values
        for(i_group=0;i_group<n_progenitors_max;i_group++){
          n_progenitors[i_file_1%n_search][i_group] =         0;
          tree_id[i_file_1%n_search][i_group]       =-i_group-2; // Leave -1 free for dropped halos
          progenitor_id[i_file_1%n_search][i_group] =-i_group-2; // Leave -1 free for dropped halos
          descendant_id[i_file_1%n_search][i_group] =-i_group-2; // Leave -1 free for dropped halos
          file_offset[i_file_1%n_search][i_group]   =-i_group-2; // Leave -1 free for dropped halos
        }
        SID_log("Done.",SID_LOG_CLOSE);

// To Do: 1) Generate match_rank and match_rank_index arrays in match_halos; also, n_match_rank array saying how many matches to each rank
//        2) During loops over n_groups, loop over ranks (in a ring) instead, processing the matches made to each rank in turn
//        3) What info is needed from each rank? projenitor_id[][]  etc

        // Find i_file_2 halos that have been matched to by i_file_1 but which do not have
        //   assigned progenitors.  They were not matched to anything when i_file_3 was processed
        //   and we need to search for them in previous outputs (done after this).
        SID_log("Searching %d %sgroups for dropped halos...",SID_LOG_OPEN,n_unprocessed,group_text_prefix);
        for(i_group=0,n_drop=0;i_group<n_groups_1;i_group++){
           // This is the index of the group in file_2 that the i_group'th group in file_1 has been matched to (-1 if no match)
           my_descendant_index=match_id[i_group];
           // If my_descendant_index>=0 then this group (and maybe more) has
           //   been matched to a descendant in the previous snapshot.
           if(my_descendant_index>=0){
              n_match_halos++;
              my_descendant=progenitor_id[i_file_2%n_search][my_descendant_index];
              // This is the (previously assigned) progenitor id of the i_group'th's descendant halo -> this halo's descendant
              // If the descendant doesn't itself have a descendant (ie. it's progenitor id is not set),
              //   then it has been dropped between file_2 and file_3.  Store information so
              //   we can deal with this afterwards.
              if(my_descendant<-1){  
                 // Check to see if this unmatched halo has already been listed ...
                 for(i_drop=0,flag_drop=TRUE;i_drop<n_drop;i_drop++){
                    if(drop_list_2[i_drop]==my_descendant_index)
                       flag_drop=FALSE;
                 }
                 // ... if not, add it to the list
                 if(flag_drop){
                    drop_list_1[n_drop]  =i_group;
                    drop_list_2[n_drop++]=my_descendant_index;
                 }
              }
           }
        }
        if(n_drop>0)
          SID_log("Done. (%d found)",SID_LOG_CLOSE,n_drop);
        else
          SID_log("Done. (none found)",SID_LOG_CLOSE);

        // Try to fix any dropped halos found above.  Any not fixed will be labeled
        //   strays and will be discarded.  This needs to be done before we attempt
        //   to fix bridges so that we have progenitor_ids for all possible mergers/bridges.
        if(n_drop>0){
          SID_log("Attempting to resolve %d dropped %sgroup candidates...",SID_LOG_OPEN|SID_LOG_TIMER,n_drop,group_text_prefix);
          //SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
          n_strays_drop=n_drop;
          for(j_file_1=i_file_2,
                j_file_2=i_file_2+2, // Manifestly: there was no match at i_file_2+1
                j_read_1=i_read_2,
                j_read_2=i_read_2+2*i_read_step,
                i_search=0,
                n_drop_found=0;
              j_file_2<=i_read_stop && i_search<n_search-2 && n_strays_drop>0;
              j_file_2++,
                j_read_2+=i_read_step,
                i_search++){
            sprintf(filename_cat2_1,"%03d",j_read_1);
            sprintf(filename_cat2_2,"%03d",j_read_2);
            
            // Read catalog to test against
            read_groups(filename_halo_root_in,j_read_2,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES,&plist3,filename_cat2_2);
            
            // Perform matching
            match_halos(&plist2,j_read_1,drop_list_2,n_strays_drop,&plist3,j_read_2,NULL,0,"drop",flag_match_subgroups);
            n_match=((int *)ADaPS_fetch(plist2.data,"n_match_drop"))[0];

            // See if we have found any of the dropped halos
            if(n_match>0){
              drop_index=(int *)ADaPS_fetch(plist2.data,"match_drop"); // These are indices of j_file_2'th groups matched to drop list
              for(i_drop=0,j_drop=0,n_sputter=0;i_drop<n_strays_drop;i_drop++){
                my_descendant_index=drop_index[i_drop];
                // If the dropped halo was found in the j_file_2'th snapshot ...
                if(my_descendant_index>=0){
                  my_descendant=progenitor_id[j_file_2%n_search][my_descendant_index];
fprintf(stderr,"%d\n",my_descendant);
                  // ... and the matched halo has a valid id, then A DROPPED HALO HAS BEEN FOUND.
                  if(my_descendant>=0){
                    // Fix the dropped halo's pointers etc.
                    n_progenitors[j_file_2%n_search][my_descendant_index]++;
                    // If it is a merger, give it a new id
                    if(n_progenitors[j_file_2%n_search][my_descendant_index]>1)
                       progenitor_id[i_file_2%n_search][drop_list_2[i_drop]]=max_id++;
                    else
                       progenitor_id[i_file_2%n_search][drop_list_2[i_drop]]=my_descendant;
                    file_offset[i_file_2%n_search][drop_list_2[i_drop]]    =j_file_2-j_file_1;
                    descendant_id[i_file_2%n_search][drop_list_2[i_drop]]  =my_descendant;
                    tree_id[i_file_2%n_search][drop_list_2[i_drop]]        =tree_id[j_file_2%n_search][my_descendant_index];
                    n_drop_found++;
                  }
                  // ... else this is a sputtering halo ... it has been dropped, found 
                  //     and dropped again multiple times now
                  else if(my_descendant==-1){
                     progenitor_id[i_file_1%n_search][drop_list_1[i_drop]] =-1;
                     progenitor_id[i_file_2%n_search][drop_list_2[i_drop]] =-1;
                     descendant_id[i_file_1%n_search][drop_list_1[i_drop]] =-1;
                     descendant_id[i_file_2%n_search][drop_list_2[i_drop]] =-1;
                     file_offset[i_file_1%n_search][drop_list_1[i_drop]]   =-1;
                     file_offset[i_file_2%n_search][drop_list_2[i_drop]]   =-1;
                     tree_id[i_file_1%n_search][drop_list_1[i_drop]]       =-1;
                     tree_id[i_file_2%n_search][drop_list_2[i_drop]]       =-1;
                     match_id[drop_list_1[i_drop]]                         =-1;
                     n_unprocessed--;
                     n_sputter++; // This is a sputtering halo ... it's been dropped multiple times
                  }
                }
                // ... else it is missing in consecutive files and we need to keep looking for it
                else{
                   drop_list_1[j_drop]=drop_list_1[i_drop];
                   drop_list_2[j_drop]=drop_list_2[i_drop];
                   j_drop++;                   
                }
              }
              n_strays_drop=j_drop;
              ADaPS_free(&(plist3.data));
            }
          }

          // If any drop-candidates were not found, label them as dropped and assign them
          //   dummy values; ignored later during output.
          for(i_stray=0;i_stray<n_strays_drop;i_stray++){
            biggest_stray_drop=MAX(biggest_stray_drop,n_particles_1[drop_list_1[i_stray]]);
            progenitor_id[i_file_1%n_search][drop_list_1[i_stray]] =-1;
            progenitor_id[i_file_2%n_search][drop_list_2[i_stray]] =-1;
            descendant_id[i_file_1%n_search][drop_list_1[i_stray]] =-1;
            descendant_id[i_file_2%n_search][drop_list_2[i_stray]] =-1;
            file_offset[i_file_1%n_search][drop_list_1[i_stray]]   =-1;
            file_offset[i_file_2%n_search][drop_list_2[i_stray]]   =-1;
            tree_id[i_file_1%n_search][drop_list_1[i_stray]]       =-1;
            tree_id[i_file_2%n_search][drop_list_2[i_stray]]       =-1;
            match_id[drop_list_1[i_stray]]                         =-1;
            n_unprocessed--;
//SID_log("%07d",SID_LOG_COMMENT,drop_list_1[i_stray]);
          }
          SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
          if(n_strays_drop>0)
            SID_log("Done. (%d strays discarded and %d are sputtering; biggest new stray=%d)",SID_LOG_CLOSE,n_strays_drop,n_sputter,biggest_stray_drop);
          else
            SID_log("Done. (no strays)",SID_LOG_CLOSE);
        }

        // We may have changed the match_ids at this point, so we need to resort
        SID_free(SID_FARG match_index);
        merge_sort((void *)match_id,(size_t)n_groups_1,&match_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE); // ... and these are their sort indices

        // Search for and process multi-matched systems (mergers or bridges) ...
/*
        if(j_file>0){

          // ... Identify systems with multiple matches and do a first pass
          //     to find bridge candidates.  We will check afterwards if there 
          //     are any mergers in the candidate bridged mergers.
          SID_log("Searching %d unprocessed %sgroups for multi-match systems...",SID_LOG_OPEN,n_unprocessed,group_text_prefix);
          //SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
          bridge_list_1=drop_list_1; // Done to make the code more readable
          bridge_list_2=drop_list_2; // Done to make the code more readable
          for(i_group=0,n_bridge_candidates=0,n_bridge_systems=0;i_group<n_groups_2;i_group++){
            // Count the number of matches between the i_group'th halo in i_file_2 and halos in i_file_1.
            j_group=find_index_int(match_id,i_group,n_groups_1,match_index);
            n_match=0;
            while(match_id[match_index[j_group+n_match]]==i_group && j_group+n_match<n_groups_1-1){
              bridge_list_1[n_bridge_candidates+n_match]=match_index[j_group+n_match];
              bridge_list_2[n_bridge_candidates+n_match]=i_group;
              n_match++;
            }
            if(match_id[match_index[j_group+n_match]]==i_group && j_group+n_match==n_groups_1-1){
              bridge_list_1[n_bridge_candidates+n_match]=match_index[j_group+n_match];
              bridge_list_2[n_bridge_candidates+n_match]=i_group;
              n_match++;
            }
            // If n_match>1, then we've got a multi-matched system...
            if(n_match>1){
              // Count the number of back matches between halos in i_file_3 and the i_group'th halo in i_file_2.
              n_back_match=0;
              k_group     =find_index_int(back_match_id,i_group,n_groups_3,back_match_index);
              while(back_match_id[back_match_index[k_group+n_back_match]]==i_group && k_group+n_back_match<n_groups_3-1)
                n_back_match++;
              if(back_match_id[back_match_index[k_group+n_back_match]]==i_group && k_group+n_back_match==n_groups_3-1)
                n_back_match++;
              // If n_back_match>1, then we've got a bridge (there might still be mergers involved as well though; deal with this after) ...
              if(n_back_match>1){
                n_bridge_candidates+=n_match;
                n_bridge_systems++;
              }
              n_multimatch+=n_match;
            }
          }
//for(i_match=0;i_match<n_multimatch;i_match++) SID_log("%07d",SID_LOG_COMMENT,bridge_list_1[i_match]);
          SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
          if(n_bridge_candidates>0)
            SID_log("Done. (%d found; %d are in %d bridged systems)",SID_LOG_CLOSE,n_multimatch,n_bridge_candidates,n_bridge_systems);
          else
            SID_log("Done. (none found)",SID_LOG_CLOSE);

          // ... Compare multi-matched halos to the results of previous
          //   back-matching to decide if they're mergers or bridges.  All scans
          //   are done over the n_groups_2 halos in file_2 but we collect the
          //   systems in file_1 which are involved.
          if(n_bridge_candidates>0){
            SID_log("Searching %d %sgroups in %d bridged systems for mergers...",SID_LOG_OPEN|SID_LOG_TIMER,n_bridge_candidates,group_text_prefix,n_bridge_systems);
            //SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
            n_strays_bridge=n_bridge_candidates;
            for(j_file_1=i_file_1,
                  j_file_2=i_file_2+1,
                  j_read_1=i_read_1,
                  j_read_2=i_read_2+i_read_step,
                  i_search=0,
                  n_found_bridge=0,
                  n_mergers_bridge=0,
                  n_bridges=0;
                j_read_2<=i_read_stop && i_search<n_search-2 && n_strays_bridge>0;
                j_file_2++,j_read_2+=i_read_step,i_search++){
              sprintf(filename_cat2_1,"%03d",j_read_1);
              sprintf(filename_cat2_2,"%03d",j_read_2);
            
              // Read catalog to test against
              read_groups(filename_halo_root_in,j_read_2,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES,&plist3,filename_cat2_2);
            
              // Perform matching
              match_halos(&plist1,j_read_1,bridge_list_1,n_strays_bridge,&plist3,j_read_2,NULL,0,"check4bridges",flag_match_subgroups);
              n_match=((int *)ADaPS_fetch(plist1.data,"n_match_check4bridges"))[0];

              // See if we have found any of the bridge candidates
              if(n_match>0){
                match_id_bridge=(int *)ADaPS_fetch(plist1.data,"match_check4bridges"); // These are indices of j_file_2'th groups matched to bridge list
                // Process matches
                for(i_bridge=0,j_bridge=0;i_bridge<n_strays_bridge;i_bridge++){
                  // If the bridge candidate was found in the j_file_2'th snapshot and the
                  //   matched halo has a descendant, then WE HAVE IDENTIFIED THE DESCENDANT
                  //   OF THIS BRIDGE CANDIDATE IN j_file_2...
                  my_descendant_index=match_id_bridge[i_bridge];
                  if(my_descendant_index>=0){
                    my_descendant=progenitor_id[j_file_2%n_search][my_descendant_index];
                    if(my_descendant>=0){
                      // If this is a match to the trunk, then call it a merger
                      //   to the trunk between i_file_1 and i_file_2
                      if(progenitor_id[i_file_2%n_search][bridge_list_2[i_bridge]]==my_descendant){
                        n_progenitors[i_file_2%n_search][bridge_list_2[i_bridge]]++;
                        if(n_progenitors[i_file_2%n_search][bridge_list_2[i_bridge]]>1){
                          progenitor_id[i_file_1%n_search][bridge_list_1[i_bridge]]=max_id++;
                          n_mergers_bridge++;
                        }
                        else
                          progenitor_id[i_file_1%n_search][bridge_list_1[i_bridge]]=my_descendant;
                        file_offset[i_file_1%n_search][bridge_list_1[i_bridge]]    =i_file_2-i_file_1;
                        descendant_id[i_file_1%n_search][bridge_list_1[i_bridge]]  =my_descendant;
                        tree_id[i_file_1%n_search][bridge_list_1[i_bridge]]        =tree_id[i_file_2%n_search][bridge_list_2[i_bridge]];
                      }
                      // ... else, this is a bridged halo and treat it as such.
                      else{
                        n_progenitors[j_file_2%n_search][my_descendant_index]++;
                        if(n_progenitors[j_file_2%n_search][my_descendant_index]>1){
                          progenitor_id[i_file_1%n_search][bridge_list_1[i_bridge]]=max_id++;
                          n_mergers_bridge++;
                        }
                        else
                          progenitor_id[i_file_1%n_search][bridge_list_1[i_bridge]]=my_descendant;
                        file_offset[i_file_1%n_search][bridge_list_1[i_bridge]]    =j_file_2-i_file_1;
                        descendant_id[i_file_1%n_search][bridge_list_1[i_bridge]]  =my_descendant;
                        tree_id[i_file_1%n_search][bridge_list_1[i_bridge]]        =tree_id[j_file_2%n_search][my_descendant_index];
                        n_bridges++;                    
                      }
                      n_unprocessed--;
                      n_found_bridge++;
                    }
                    // ... else it is a dropped halo and we need to continue to
                    //     look in the next file (or give-up if we have scanned
                    //     n_match catalogs, and call it an stray; handle after)
                    else{
                      bridge_list_1[j_bridge]=bridge_list_1[i_bridge];
                      bridge_list_2[j_bridge]=bridge_list_2[i_bridge];
                      j_bridge++;
                    }
                  }
                  else{
                    bridge_list_1[j_bridge]=bridge_list_1[i_bridge];
                    bridge_list_2[j_bridge]=bridge_list_2[i_bridge];
                    j_bridge++;
                  }
                }
                n_strays_bridge=j_bridge;
                ADaPS_free(&(plist3.data));
              }
            }

            // If any bridge-candidate halos were not found in previous snaps, call them mergers
            //   to the trunk (ie. their original matching)
            for(i_bridge=0;i_bridge<n_strays_bridge;i_bridge++){
               n_progenitors[i_file_2%n_search][bridge_list_2[i_bridge]]++;
               if(n_progenitors[i_file_2%n_search][bridge_list_2[i_bridge]]>1){
                 progenitor_id[i_file_1%n_search][bridge_list_1[i_bridge]]=max_id++;
                 n_mergers_bridge++;
               }
               else
                 progenitor_id[i_file_1%n_search][bridge_list_1[i_bridge]]=my_descendant;
              n_unprocessed--;
            }
            n_strays_bridge=0;

            SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
            SID_log("Done. (%d were called mergers, %d were called bridges)",SID_LOG_CLOSE,n_mergers_bridge,n_bridges);
          }
        }
        */

        // Propagate progenitor and tree ids for the simple matches
        SID_log("Assigning IDs to %d %sgroups...",SID_LOG_OPEN,n_unprocessed,group_text_prefix);
        for(i_group=0,j_group=0,n_mergers=0;i_group<n_groups_1;i_group++){
          // This is the index of the group in file_2 that the i_group'th group in file_1 has been matched to (-1 if no match)
          my_descendant_index=match_id[i_group];

          // If my_descendant_index>=0 then this group (and maybe more) has
          //   been matched to a descendant in the previous snapshot.  If
          //   progenitor_id==-1, it was dropped and if it is less than 
          //   that, it has not been matched to anything yet
          if(my_descendant_index>=0){
            my_descendant=progenitor_id[i_file_2%n_search][my_descendant_index];
            if(my_descendant>=0 && progenitor_id[i_file_1%n_search][i_group]<(-1)){
              n_progenitors[i_file_2%n_search][my_descendant_index]++;
              if(n_progenitors[i_file_2%n_search][my_descendant_index]>1){
                progenitor_id[i_file_1%n_search][i_group]=max_id++;
                n_mergers++;
              }
              else
                progenitor_id[i_file_1%n_search][i_group]=my_descendant;                  
              file_offset[i_file_1%n_search][i_group]    =1;
              descendant_id[i_file_1%n_search][i_group]  =my_descendant;
              tree_id[i_file_1%n_search][i_group]        =tree_id[i_file_2%n_search][my_descendant_index];
              n_unprocessed--;
            }
            else if(progenitor_id[i_file_1%n_search][i_group]<(-1)){ // Don't set it's progenitor_ID because we need the <-1 value for looking for dropped halos
              file_offset[i_file_1%n_search][i_group]  =-1;
              //descendant_id[i_file_1%n_search][i_group]=-1;
              tree_id[i_file_1%n_search][i_group]      =-1;            
              n_unprocessed--;
            }
          }
          // This halo was not matched to anything
          else if(progenitor_id[i_file_1%n_search][i_group]!=(-1)){ // Don't re-process strays
             progenitor_id[i_file_1%n_search][i_group]=-1;
             file_offset[i_file_1%n_search][i_group]  =-1;
             //descendant_id[i_file_1%n_search][i_group]=-1;
             tree_id[i_file_1%n_search][i_group]      =-1;            
             n_unprocessed--;
          }
        }
        SID_log("Done. (%d mergers)",SID_LOG_CLOSE,n_mergers);

        // Generate sort indices for progenitor ids
        SID_free((void **)&sort_id[i_file_1%n_search]);
        merge_sort(progenitor_id[i_file_1%n_search],
                   n_groups_1,
                   &(sort_id[i_file_1%n_search]),
                   SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

        // Free match sort indices
        SID_free((void **)&match_index);
        if(j_file>0)
          SID_free((void **)&back_match_index);

        if(n_unprocessed!=0)
          SID_trap_error("n_unprocessed=%d when it should be zero!",ERROR_LOGIC,n_unprocessed);
      }
     
      // Write a log file
      n_mergers   +=n_mergers_bridge;
      n_found      =n_drop_found+n_found_bridge;
      n_strays     =n_strays_drop+n_strays_bridge;
      biggest_stray=MAX(biggest_stray_drop,biggest_stray_bridge);
      if(k_match==0){
        if(j_file==0){
          fp=fopen(filename_log,"w");
          fprintf(fp,"# (1):  filenumber\n");
          fprintf(fp,"# (2):  Maximum subgroup ID\n");
          fprintf(fp,"# (3):  # of subgroups\n");
          fprintf(fp,"# (4):  # of matched subgroups\n");
          fprintf(fp,"# (5):  # of dropped subgroups\n");
          fprintf(fp,"# (6):  # of dropped subgroups found\n");
          fprintf(fp,"# (7):  # of sputtering subgroups\n");
          fprintf(fp,"# (8):  # of subgroups in multiple-match systems (MMSs)\n");
          fprintf(fp,"# (9):  # of MMSs that are bridge candidates\n");
          fprintf(fp,"# (10): # of MMSs that ARE in bridge systems\n");
          fprintf(fp,"# (11): # of MMSs that are mergers in bridges\n");
          fprintf(fp,"# (12): # of subgroup mergers (total)\n");
          fprintf(fp,"# (13): # of stray subgroups\n");
          fprintf(fp,"# (14): biggest strayed subgroup\n");
          if(n_k_match>1){
             fprintf(fp,"# (15): Maximum group ID\n");
             fprintf(fp,"# (16): # of groups\n");
             fprintf(fp,"# (17): # of matched groups\n");
             fprintf(fp,"# (18): # of dropped groups\n");
             fprintf(fp,"# (19): # of dropped groups found\n");
             fprintf(fp,"# (20): # of sputtering subgroups\n");
             fprintf(fp,"# (21): # of groups in multiple-match systems (MMSs)\n");
             fprintf(fp,"# (22): # of MMSs that are bridge candidates\n");
             fprintf(fp,"# (23): # of MMSs that ARE in bridge systems\n");
             fprintf(fp,"# (24): # of MMSs that are mergers in bridges\n");
             fprintf(fp,"# (25): # of group mergers (total)\n");
             fprintf(fp,"# (26): # of stray groups\n");
             fprintf(fp,"# (27): biggest strayed group\n");
          }
        }
        else
          fp=fopen(filename_log,"a");
        fprintf(fp,"%4d",i_read);
      }
      else
        fp=fopen(filename_log,"a");
      fprintf(fp," %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d",
              max_id-1,
              n_groups_1,
              n_match_halos,
              n_drop,
              n_drop_found,
              n_sputter,
              n_multimatch,
              n_bridge_candidates,
              n_bridge_systems,
              n_mergers_bridge,
              n_bridges,
              n_mergers,
              n_strays,
              biggest_stray);
      if(k_match==n_k_match-1)
        fprintf(fp,"\n");
      fclose(fp);

      switch(k_match){
      case 0:
        max_id_subgroup=max_id;
        break;
      case 1:
        max_id_group=max_id;
        break;
      }      
      SID_log("Done.",SID_LOG_CLOSE);
    }

    // Write tree info once a few files have been processed
    //   and no more dropped groups need to be given ids
    if(j_file>=n_search-2){
      SID_log("Writing results for snapshot #%d...",SID_LOG_OPEN,k_write);
      sprintf(filename_matches_out,            "%s.trees_horizontal_%d",                     filename_root_out,k_write);
      sprintf(filename_group_properties_out,   "%s.trees_horizontal_groups_properties_%d",   filename_root_out,k_write);
      sprintf(filename_subgroup_properties_out,"%s.trees_horizontal_subgroups_properties_%d",filename_root_out,k_write);
      sprintf(filename_group_properties_in,    "%s_%03d.catalog_groups_properties",          filename_cat_root_in, k_write);
      sprintf(filename_subgroup_properties_in, "%s_%03d.catalog_subgroups_properties",       filename_cat_root_in, k_write);
      fp_matches_out=fopen(filename_matches_out,"w");
      SID_fopen(filename_group_properties_out,   "w",&fp_group_properties_out);
      SID_fopen(filename_subgroup_properties_out,"w",&fp_subgroup_properties_out);
      fprintf(fp_matches_out,"%d %d %d %d %d\n",n_groups[i_write%n_search],n_subgroups[i_write%n_search],n_progenitors_max,max_tree_id_subgroup,max_tree_id_group); // max_tree_ids' are actually n_tree_ids ... so no -1
      if(n_groups[i_write%n_search]>0){
        fp_group_properties_in   =fopen(filename_group_properties_in,"r");
        fp_subgroup_properties_in=fopen(filename_subgroup_properties_in,"r");
        fseeko(fp_group_properties_in,   4*sizeof(int),SEEK_CUR);
        fseeko(fp_subgroup_properties_in,4*sizeof(int),SEEK_CUR);
        for(i_group=0,i_subgroup=0;i_group<n_groups[i_write%n_search];i_group++){
          fprintf(fp_matches_out,"%d %d %d %d %d",progenitor_id_group[i_write%n_search][i_group],descendant_id_group[i_write%n_search][i_group],tree_id_group[i_write%n_search][i_group],file_offset_group[i_write%n_search][i_group],n_subgroups_group[i_write%n_search][i_group]);
          read_group_properties(fp_group_properties_in,&properties,i_group,k_write);
          if(progenitor_id_group[i_write%n_search][i_group]>=0)
            SID_fwrite(&properties,sizeof(halo_info),1,&fp_group_properties_out);
          for(j_subgroup=0;j_subgroup<n_subgroups_group[i_write%n_search][i_group];j_subgroup++,i_subgroup++){
            fprintf(fp_matches_out,"  %d %d %d %d",progenitor_id_subgroup[i_write%n_search][i_subgroup],descendant_id_subgroup[i_write%n_search][i_subgroup],tree_id_subgroup[i_write%n_search][i_subgroup],file_offset_subgroup[i_write%n_search][i_subgroup]);          
            read_group_properties(fp_subgroup_properties_in,&properties,i_subgroup,k_write);
            if(progenitor_id_subgroup[i_write%n_search][i_subgroup]>=0)
              SID_fwrite(&properties,sizeof(halo_info),1,&fp_subgroup_properties_out);
          }
          fprintf(fp_matches_out,"\n");
        }
        fclose(fp_group_properties_in);
        fclose(fp_subgroup_properties_in);
      }
      fclose(fp_matches_out);
      SID_fclose(&fp_group_properties_out);
      SID_fclose(&fp_subgroup_properties_out);
      i_write--;
      j_write++;
      k_write-=i_read_step;
      SID_log("Done.",SID_LOG_CLOSE);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
  free_plist(&plist1);
  free_plist(&plist2);
  free_plist(&plist3);

  // Write the remaining tree info
  for(;k_write>=i_read_start;i_write--,j_write++,k_write-=i_read_step){
    SID_log("Writing results for snapshot #%d...",SID_LOG_OPEN,k_write);
    sprintf(filename_matches_out,            "%s.trees_horizontal_%d",                     filename_root_out,   k_write);
    sprintf(filename_group_properties_out,   "%s.trees_horizontal_groups_properties_%d",   filename_root_out,   k_write);
    sprintf(filename_subgroup_properties_out,"%s.trees_horizontal_subgroups_properties_%d",filename_root_out,   k_write);
    sprintf(filename_group_properties_in,    "%s_%03d.catalog_groups_properties",          filename_cat_root_in,k_write);
    sprintf(filename_subgroup_properties_in, "%s_%03d.catalog_subgroups_properties",       filename_cat_root_in,k_write);
    fp_matches_out=fopen(filename_matches_out,"w");
    SID_fopen(filename_group_properties_out,   "w",&fp_group_properties_out);
    SID_fopen(filename_subgroup_properties_out,"w",&fp_subgroup_properties_out);
    fprintf(fp_matches_out,"%d %d %d %d %d\n",n_groups[i_write%n_search],n_subgroups[i_write%n_search],n_progenitors_max,max_tree_id_subgroup,max_tree_id_group); // max_tree_ids' are actually n_tree_ids ... so no -1
    if(n_groups[i_write%n_search]>0){
      fp_group_properties_in   =fopen(filename_group_properties_in,"r");
      fp_subgroup_properties_in=fopen(filename_subgroup_properties_in,"r");
      fseeko(fp_group_properties_in,   4*sizeof(int),SEEK_CUR);
      fseeko(fp_subgroup_properties_in,4*sizeof(int),SEEK_CUR);
      for(i_group=0,i_subgroup=0;i_group<n_groups[i_write%n_search];i_group++){
        fprintf(fp_matches_out,"%d %d %d %d %d",progenitor_id_group[i_write%n_search][i_group],descendant_id_group[i_write%n_search][i_group],tree_id_group[i_write%n_search][i_group],file_offset_group[i_write%n_search][i_group],n_subgroups_group[i_write%n_search][i_group]);
        read_group_properties(fp_group_properties_in,&properties,i_group,k_write);        
        if(progenitor_id_group[i_write%n_search][i_group]>=0)
          SID_fwrite(&properties,sizeof(halo_info),1,&fp_group_properties_out);
        for(j_subgroup=0;j_subgroup<n_subgroups_group[i_write%n_search][i_group];j_subgroup++,i_subgroup++){
          fprintf(fp_matches_out,"  %d %d %d %d",progenitor_id_subgroup[i_write%n_search][i_subgroup],descendant_id_subgroup[i_write%n_search][i_subgroup],tree_id_subgroup[i_write%n_search][i_subgroup],file_offset_subgroup[i_write%n_search][i_subgroup]);          
          read_group_properties(fp_subgroup_properties_in,&properties,i_subgroup,k_write);        
          if(progenitor_id_subgroup[i_write%n_search][i_subgroup]>=0)
            SID_fwrite(&properties,sizeof(halo_info),1,&fp_subgroup_properties_out);
        }
        fprintf(fp_matches_out,"\n");
      }
      fclose(fp_group_properties_in);
      fclose(fp_subgroup_properties_in);
    }
    fclose(fp_matches_out);
    SID_fclose(&fp_group_properties_out);
    SID_fclose(&fp_subgroup_properties_out);
    SID_log("Done.",SID_LOG_CLOSE);
  }

  // Clean-up
  SID_log("Freeing arrays...",SID_LOG_OPEN);
  for(i_search=0;i_search<n_search;i_search++){
    SID_free((void **)&tree_id_subgroup[i_search]);
    SID_free((void **)&tree_id_group[i_search]);
    SID_free((void **)&n_progenitors_subgroup[i_search]);
    SID_free((void **)&n_progenitors_group[i_search]);
    SID_free((void **)&progenitor_id_subgroup[i_search]);
    SID_free((void **)&progenitor_id_group[i_search]);
    SID_free((void **)&descendant_id_subgroup[i_search]);
    SID_free((void **)&descendant_id_group[i_search]);
    SID_free((void **)&file_offset_subgroup[i_search]);
    SID_free((void **)&file_offset_group[i_search]);
    SID_free((void **)&n_subgroups_group[i_search]);
    SID_free((void **)&sort_subgroup_id[i_search]);
    SID_free((void **)&sort_group_id[i_search]);
  }
  SID_free((void **)&tree_id_subgroup);
  SID_free((void **)&tree_id_group);
  SID_free((void **)&n_progenitors_subgroup);
  SID_free((void **)&n_progenitors_group);
  SID_free((void **)&progenitor_id_subgroup);
  SID_free((void **)&progenitor_id_group);
  SID_free((void **)&descendant_id_subgroup);
  SID_free((void **)&descendant_id_group);
  SID_free((void **)&file_offset_subgroup);
  SID_free((void **)&file_offset_group);
  SID_free((void **)&n_subgroups_group);
  SID_free((void **)&sort_subgroup_id);
  SID_free((void **)&sort_group_id);
  SID_free((void **)&drop_list_1);
  SID_free((void **)&drop_list_2);
  SID_free((void **)&n_groups);
  SID_free((void **)&n_subgroups);
  SID_log("Done.",SID_LOG_CLOSE);

  // Set flag_clean=TRUE so that the output files generated here
  //  will be treated as temporary in the next step (if called)
  *flag_clean=TRUE;

  SID_log("Done.",SID_LOG_CLOSE);

}
