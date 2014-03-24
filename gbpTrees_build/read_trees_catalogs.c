#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void read_trees_catalogs(tree_info *trees,
                         char      *filename_SSimPL_dir,
                         char      *filename_catalog_name,
                         int        mode){
  int i_read;
  int j_read;
  int i_snap;
  int j_snap;
  int k_snap;
  int i_neighbour;

  SID_log("Reading tree properties...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Create the data array(s) where stuff will be stored
  halo_properties_SAGE_info **SAGE_properties_groups_local;
  halo_properties_info      **properties_groups_local;
  halo_profile_info         **profiles_groups_local;
  halo_properties_SAGE_info **SAGE_properties_subgroups_local;
  halo_properties_info      **properties_subgroups_local;
  halo_profile_info         **profiles_subgroups_local;
  int                    read_catalog_mode=READ_CATALOG_PROPERTIES;
  if(check_mode_for_flag(mode,READ_TREES_CATALOGS_GROUPS)){
     if(check_mode_for_flag(mode,READ_TREES_CATALOGS_SAGE)){
        init_trees_data(trees,(void ***)&SAGE_properties_groups_local,sizeof(halo_properties_SAGE_info),INIT_TREE_DATA_GROUPS,"properties_groups_SAGE");
        properties_groups_local=NULL;
        trees->group_properties_SAGE=SAGE_properties_groups_local;
     }
     else{
        init_trees_data(trees,(void ***)&properties_groups_local,sizeof(halo_properties_info),INIT_TREE_DATA_GROUPS,"properties_groups");
        SAGE_properties_groups_local=NULL;
        trees->group_properties     =properties_groups_local;
     }
     if(check_mode_for_flag(mode,READ_TREES_CATALOGS_PROFILES)){
        init_trees_data(trees,(void ***)&profiles_groups_local,sizeof(halo_profile_info),INIT_TREE_DATA_GROUPS,"profiles_groups");
        read_catalog_mode=read_catalog_mode|READ_CATALOG_PROFILES;
     }
     else
        profiles_groups_local=NULL;
  }
  if(check_mode_for_flag(mode,READ_TREES_CATALOGS_SUBGROUPS)){
     if(check_mode_for_flag(mode,READ_TREES_CATALOGS_SAGE)){
        init_trees_data(trees,(void ***)&SAGE_properties_subgroups_local,sizeof(halo_properties_SAGE_info),INIT_TREE_DATA_SUBGROUPS,"properties_subgroups_SAGE");
        properties_subgroups_local     =NULL;
        trees->subgroup_properties_SAGE=SAGE_properties_subgroups_local;
     }
     else{
        init_trees_data(trees,(void ***)&properties_subgroups_local,sizeof(halo_properties_info),INIT_TREE_DATA_SUBGROUPS,"properties_subgroups");
        SAGE_properties_subgroups_local=NULL;
        trees->subgroup_properties     =properties_subgroups_local;
     }
     if(check_mode_for_flag(mode,READ_TREES_CATALOGS_PROFILES)){
        init_trees_data(trees,(void ***)&profiles_subgroups_local,sizeof(halo_profile_info),INIT_TREE_DATA_SUBGROUPS,"profiles_subgroups");
        read_catalog_mode=read_catalog_mode|READ_CATALOG_PROFILES;
     }
     else
        profiles_subgroups_local=NULL;
  }

  // Process each snapshot in turn
  int    *nebr_idx_list_local;
  int    *file_idx_list_local;
  int    *list_init_local;
  size_t *file_idx_list_local_index;
  int     n_list_local;
  int     n_alloc=MAX(trees->max_n_subgroups_snap_local,trees->max_n_groups_snap_local);
  file_idx_list_local=(int *)SID_malloc(sizeof(int)*n_alloc);
  nebr_idx_list_local=(int *)SID_malloc(sizeof(int)*n_alloc);
  list_init_local    =(int *)SID_malloc(sizeof(int)*n_alloc);
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
     i_read=trees->snap_list[i_snap];
     SID_log("Processing snapshot %03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read);

     // Loop twice; once for subgroups and then once for groups
     for(int i_type=0;i_type<2;i_type++){
        // Set some group/subgroup specific things
        int                   n_halos;
        tree_node_info       *first_neighbour;
        char                  group_text_prefix[5];
        int                   open_catalog_mode;
        halo_properties_SAGE_info            *SAGE_properties_in;
        halo_properties_info *properties_in;
        halo_profile_info    *profiles_in;
        int                   flag_proceed=TRUE;
        switch(i_type){
           case 0:
              flag_proceed=check_mode_for_flag(mode,READ_TREES_CATALOGS_SUBGROUPS);
              if(flag_proceed){
                 sprintf(group_text_prefix,"sub");
                 open_catalog_mode  =READ_CATALOG_SUBGROUPS|read_catalog_mode;
                 n_halos            =trees->n_subgroups_snap_local[i_snap];
                 first_neighbour    =trees->first_neighbour_subgroups[i_snap];
                 if(SAGE_properties_subgroups_local!=NULL)
                    SAGE_properties_in=SAGE_properties_subgroups_local[i_snap];
                 else
                    SAGE_properties_in=NULL;
                 if(properties_subgroups_local!=NULL)
                    properties_in=properties_subgroups_local[i_snap];
                 else
                    properties_in=NULL;
                 if(profiles_subgroups_local!=NULL)
                    profiles_in=profiles_subgroups_local[i_snap];
                 else
                    profiles_in=NULL;
              }
              break;
           case 1:
              flag_proceed=check_mode_for_flag(mode,READ_TREES_CATALOGS_GROUPS);
              if(flag_proceed){
                 sprintf(group_text_prefix,"");
                 open_catalog_mode  =READ_CATALOG_GROUPS|read_catalog_mode;
                 n_halos            =trees->n_groups_snap_local[i_snap];
                 first_neighbour    =trees->first_neighbour_groups[i_snap];
                 if(SAGE_properties_groups_local!=NULL)
                    SAGE_properties_in=SAGE_properties_groups_local[i_snap];
                 else
                    SAGE_properties_in=NULL;
                 if(properties_groups_local!=NULL)
                    properties_in=properties_groups_local[i_snap];
                 else
                    properties_in=NULL;
                 if(profiles_groups_local!=NULL)
                    profiles_in=profiles_groups_local[i_snap];
                 else
                    profiles_in=NULL;
              }
              break;
        }
        if(flag_proceed){
           // Initialize the validation array
           for(i_neighbour=0;i_neighbour<n_halos;i_neighbour++) 
              list_init_local[i_neighbour]=0;

           // Create a sorted list of locally stored halos
           tree_node_info *current;
           n_list_local=0;
           current     =first_neighbour;
           while(current!=NULL){
              file_idx_list_local[n_list_local]=current->file_index;
              nebr_idx_list_local[n_list_local]=current->neighbour_index;
              n_list_local++;
              current=current->next_neighbour;
           }
           if(n_list_local!=n_halos)
              SID_trap_error("Count mismatch for %sgroups in read_trees_properties() (ie. %d!=%d)",ERROR_LOGIC,
                             group_text_prefix,n_list_local,n_halos);
           merge_sort(file_idx_list_local,(size_t)n_list_local,&file_idx_list_local_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);

           // Set filenames and open files
           fp_catalog_info fp_properties;
           char   filename_cat_root_in[256];
           sprintf(filename_cat_root_in,"%s/catalogs/%s",
                                        filename_SSimPL_dir,filename_catalog_name);
           fopen_catalog(filename_cat_root_in,
                         i_read,
                         open_catalog_mode,
                         &fp_properties);

           // Create some read buffer structures
           halo_properties_SAGE_info            *SAGE_properties_buffer=NULL;
           halo_properties_info *properties_buffer      =NULL;
           halo_profile_info    *profiles_buffer        =NULL;
           if(SAGE_properties_in!=NULL)
              SAGE_properties_buffer=(halo_properties_SAGE_info *)SID_malloc(sizeof(halo_properties_SAGE_info));
           else
              SAGE_properties_buffer=NULL;
           if(properties_in!=NULL)
              properties_buffer=(halo_properties_info *)SID_malloc(sizeof(halo_properties_info));
           else
              properties_buffer=NULL;
           if(profiles_in!=NULL)
              profiles_buffer=(halo_profile_info *)SID_malloc(sizeof(halo_profile_info));
           else
              profiles_buffer=NULL;

           // Perform read
           int k_read;
           int l_read;
           for(k_read=0,l_read=0;k_read<fp_properties.n_halos_total;k_read++){
              fread_catalog_file(&fp_properties,SAGE_properties_buffer,properties_buffer,profiles_buffer,k_read);
              if(l_read<n_list_local){
                 if(k_read==file_idx_list_local[file_idx_list_local_index[l_read]]){
                    if(SAGE_properties_buffer!=NULL)
                       memcpy(&(SAGE_properties_in[nebr_idx_list_local[file_idx_list_local_index[l_read]]]),SAGE_properties_buffer,sizeof(halo_properties_SAGE_info));
                    if(properties_buffer!=NULL)
                       memcpy(&(properties_in[nebr_idx_list_local[file_idx_list_local_index[l_read]]]),properties_buffer,sizeof(halo_properties_info));
                    if(profiles_buffer!=NULL)
                       memcpy(&(profiles_in[nebr_idx_list_local[file_idx_list_local_index[l_read]]]),profiles_buffer,sizeof(halo_profile_info));
                    list_init_local[nebr_idx_list_local[file_idx_list_local_index[l_read]]]++;
                    l_read++;
                 }
              }
           }
           SID_free(SID_FARG SAGE_properties_buffer);
           SID_free(SID_FARG properties_buffer);
           SID_free(SID_FARG profiles_buffer);

           // Check that everything was read properly
           if(l_read!=n_list_local)
              SID_trap_error("An incorrect number of matches were read (ie. %d!=%d) for snaphot %d",ERROR_LOGIC,
                             l_read,n_list_local,i_read);
           for(i_neighbour=0;i_neighbour<n_halos;i_neighbour++)
              if(list_init_local[i_neighbour]!=1 && list_init_local[i_neighbour]>=0)
                 SID_trap_error("Neighbour %d of %d was not processed correctly (init=%d) for snapshot %d.",ERROR_LOGIC,
                                i_neighbour,n_halos,list_init_local[i_neighbour],i_read);

           // Clean-up
           fclose_catalog(&fp_properties);
           SID_free(SID_FARG file_idx_list_local_index);
        } // if(flag_proceed)
     } // i_type

     SID_log("Done.",SID_LOG_CLOSE);
  }

  // Clean-up
  SID_free(SID_FARG nebr_idx_list_local);
  SID_free(SID_FARG file_idx_list_local);
  SID_free(SID_FARG list_init_local);

  SID_log("Done.",SID_LOG_CLOSE);
}

