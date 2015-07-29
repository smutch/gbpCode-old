#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void read_trees_data(tree_info    *trees,
                     char         *filename_root,
                     int           mode,
                     SID_Datatype  data_type,
                     const char   *name){
  int i_read;
  int j_read;
  int i_snap;
  int j_snap;
  int k_snap;
  int i_neighbour;

  SID_log("Reading %s...",SID_LOG_OPEN|SID_LOG_TIMER,name);

  // Create the data array(s) where stuff will be stored
  void **data_groups_local   =NULL;
  void **data_subgroups_local=NULL;
  int    data_type_size;
  SID_Type_size(data_type,&data_type_size);
  if(check_mode_for_flag(mode,READ_TREES_DATA_GROUPS))
     init_trees_data(trees,&data_groups_local,data_type_size,INIT_TREE_DATA_GROUPS,"%s_groups",name);
  if(check_mode_for_flag(mode,READ_TREES_DATA_SUBGROUPS))
     init_trees_data(trees,&data_subgroups_local,data_type_size,INIT_TREE_DATA_SUBGROUPS,"%s_subgroups",name);

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
     void *data_read=SID_malloc(data_type_size);
     for(int i_type=0;i_type<2;i_type++){
        // Set some group/subgroup specific things
        int             n_halos;
        int             n_halos_all;
        tree_node_info *first_neighbour;
        char            group_text_prefix[5];
        int             open_catalog_mode;
        char           *data_in;
        int             flag_proceed=TRUE;
        switch(i_type){
           case 0:
              flag_proceed=check_mode_for_flag(mode,READ_TREES_DATA_SUBGROUPS);
              if(flag_proceed){
                 sprintf(group_text_prefix,"sub");
                 n_halos            =trees->n_subgroups_snap_local[i_snap];
                 n_halos_all        =trees->n_subgroups_catalog[i_snap];
                 first_neighbour    =trees->first_neighbour_subgroups[i_snap];
                 if(data_subgroups_local!=NULL)
                    data_in=(char *)data_subgroups_local[i_snap];
                 else
                    data_in=NULL;
              }
              break;
           case 1:
              flag_proceed=check_mode_for_flag(mode,READ_TREES_DATA_GROUPS);
              if(flag_proceed){
                 sprintf(group_text_prefix,"");
                 n_halos            =trees->n_groups_snap_local[i_snap];
                 n_halos_all        =trees->n_groups_catalog[i_snap];
                 first_neighbour    =trees->first_neighbour_groups[i_snap];
                 if(data_groups_local!=NULL)
                    data_in=(char *)data_groups_local[i_snap];
                 else
                    data_in=NULL;
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
              SID_trap_error("Count mismatch for %sgroups in read_trees_data() (ie. %d!=%d)",ERROR_LOGIC,
                             group_text_prefix,n_list_local,n_halos);
           merge_sort(file_idx_list_local,(size_t)n_list_local,&file_idx_list_local_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);

           // Set filename and open file
           char  filename_cat_in[256];
           sprintf(filename_cat_in,"%s_%03d.%sgroups",filename_root,trees->snap_list[i_snap],group_text_prefix);
           FILE *fp_data=fopen(filename_cat_in,"r");

           if(fp_data!=NULL){
              // Read header
              int i_file;
              int n_files;
              int n_halos_file;
              int n_halos_total;
              fread(&i_file,       sizeof(int),1,fp_data);
              fread(&n_files,      sizeof(int),1,fp_data);
              fread(&n_halos_file, sizeof(int),1,fp_data);
              fread(&n_halos_total,sizeof(int),1,fp_data);
//              if(n_halos_total!=n_halos_all)
//                 SID_trap_error("Count of halos in data file does not match trees (ie. %d!=%d).",ERROR_LOGIC,n_halos_total,n_halos_all);

              // Perform read
              int k_read;
              int l_read;
              for(k_read=0,l_read=0;k_read<n_halos_total;k_read++){
                 fread(data_read,data_type_size,1,fp_data);
                 if(l_read<n_list_local){
                    if(k_read==file_idx_list_local[file_idx_list_local_index[l_read]]){
                       if(data_in!=NULL)
                          memcpy(&(data_in[data_type_size*nebr_idx_list_local[file_idx_list_local_index[l_read]]]),data_read,data_type_size);
                       list_init_local[nebr_idx_list_local[file_idx_list_local_index[l_read]]]++;
                       l_read++;
                    }
                 }
              }
              fclose(fp_data);

              // Check that everything was read properly
if(l_read!=n_list_local) SID_log_warning("Invalid read count: %d!=%d",SID_WARNING_DEFAULT,l_read,n_list_local);
//              if(l_read!=n_list_local)
//                 SID_trap_error("An incorrect number of matches were read (ie. %d!=%d) for snaphot %d",ERROR_LOGIC,
//                                l_read,n_list_local,i_read);
//              for(i_neighbour=0;i_neighbour<n_halos;i_neighbour++)
//                 if(list_init_local[i_neighbour]!=1 && list_init_local[i_neighbour]>=0)
//                    SID_trap_error("Neighbour %d of %d was not processed correctly (init=%d) for snapshot %d.",ERROR_LOGIC,
//                                   i_neighbour,n_halos,list_init_local[i_neighbour],i_read);
           }
           else{
              if(data_in!=NULL){
                 for(int k_read=0;k_read<n_halos*data_type_size;k_read++)
                    data_in[k_read]=0;
              }
           }

           // Clean-up
           SID_free(SID_FARG file_idx_list_local_index);
        } // if(flag_proceed)
     } // i_type
     SID_free(SID_FARG data_read);

     SID_log("Done.",SID_LOG_CLOSE);
  }

  // Clean-up
  SID_free(SID_FARG nebr_idx_list_local);
  SID_free(SID_FARG file_idx_list_local);
  SID_free(SID_FARG list_init_local);

  SID_log("Done.",SID_LOG_CLOSE);
}

