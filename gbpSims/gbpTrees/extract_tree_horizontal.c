#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees.h>

int main(int argc,char *argv[]){
  SID_init(&argc,&argv,NULL,NULL);

  char filename_root[MAX_FILENAME_LENGTH];
  char filename_catalog_root[MAX_FILENAME_LENGTH];
  char filename_out_root[MAX_FILENAME_LENGTH];
  char filename_out_groups[MAX_FILENAME_LENGTH];
  char filename_out_subgroups[MAX_FILENAME_LENGTH];
  char filename_out_list[MAX_FILENAME_LENGTH];
  int  snap_start,snap_stop,snap_step;
  int  tree_type;

  if(argc!=6)
     SID_trap_error("Invalid syntax.",ERROR_SYNTAX);
  int i_tree_extract_group   =-1;
  int i_tree_extract_subgroup=-1;
  strcpy(filename_root,        argv[1]);
  strcpy(filename_catalog_root,argv[4]);
  strcpy(filename_out_root,    argv[5]);
  FILE *fp_out_groups   =NULL;
  FILE *fp_out_subgroups=NULL;
  if(!strcmp(argv[2],"group") ||
     !strcmp(argv[2],"groups") ||
     !strcmp(argv[2],"GROUP") ||
     !strcmp(argv[2],"GROUPS")){
     i_tree_extract_group=atoi(argv[3]);
     sprintf(filename_out_groups,"%s_group_tree_%09d.ascii",filename_out_root,i_tree_extract_group);
     sprintf(filename_out_list,  "%s_group_forest.ascii",   filename_out_root);
     fp_out_groups=fopen(filename_out_groups,"w");
  }
  else if(!strcmp(argv[2],"subgroup") ||
          !strcmp(argv[2],"subgroups") ||
          !strcmp(argv[2],"SUBGROUP") ||
          !strcmp(argv[2],"SUBGROUPS")){
     i_tree_extract_subgroup=atoi(argv[3]);
     sprintf(filename_out_subgroups,"%s_subgroup_tree_%09d.ascii",filename_out_root,i_tree_extract_subgroup);
     sprintf(filename_out_list,  "%s_subgroup_forest.ascii",      filename_out_root);
     fp_out_subgroups=fopen(filename_out_subgroups,"w");
  }
  else
     SID_trap_error("Invalid tree type specified {%s}.",ERROR_SYNTAX,argv[2]);
 
  SID_log("Extracting horizontal tree...",SID_LOG_OPEN|SID_LOG_TIMER);
  char   filename_file_root[MAX_FILENAME_LENGTH];
  char   filename_dir_horizontal[MAX_FILENAME_LENGTH];
  char   filename_dir_horizontal_trees[MAX_FILENAME_LENGTH];
  char   filename_dir_horizontal_groups[MAX_FILENAME_LENGTH];
  char   filename_dir_horizontal_subgroups[MAX_FILENAME_LENGTH];
  char   filename_in[MAX_FILENAME_LENGTH];
  SID_fp fp_in;
  strcpy(filename_file_root,filename_root);
  strip_path(filename_file_root);

  int i_read;
  int i_read_start;
  int i_read_stop;
  int i_read_step;
  int n_search;
  int flag_fix_bridges;
  int flag_compute_fragmented;
  int flag_compute_ghosts;
  read_trees_run_parameters(filename_root,
                            &i_read_start,
                            &i_read_stop,
                            &i_read_step,
                            &n_search,
                            &flag_fix_bridges,
                            &flag_compute_fragmented,
                            &flag_compute_ghosts);

  // Fetch the final number of trees
  int i_read_last=i_read_stop;
  while(i_read_last>=n_search) i_read_last-=n_search;
  sprintf(filename_dir_horizontal,      "%s/horizontal",filename_root);
  sprintf(filename_dir_horizontal_trees,"%s/trees",     filename_dir_horizontal);
  sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_dir_horizontal_trees,i_read_last);
  SID_fopen(filename_in,"r",&fp_in);
  int n_step_in;
  int n_search_in;
  int n_groups_max_in;
  int n_subgroups_max_in;
  int n_groups;
  int n_subgroups;
  int n_trees_subgroup_total;
  int n_trees_group_total;
  SID_fread_all(&n_step_in,         sizeof(int),1,&fp_in);
  SID_fread_all(&n_search_in,       sizeof(int),1,&fp_in);
  SID_fread_all(&n_groups,          sizeof(int),1,&fp_in);
  SID_fread_all(&n_subgroups,       sizeof(int),1,&fp_in);
  SID_fread_all(&n_groups_max_in,   sizeof(int),1,&fp_in);
  SID_fread_all(&n_subgroups_max_in,sizeof(int),1,&fp_in);
  SID_fread_all(&n_trees_subgroup_total,sizeof(int),1,&fp_in);
  SID_fread_all(&n_trees_group_total,   sizeof(int),1,&fp_in);
  SID_fclose(&fp_in);
  SID_log("n_trees_groups    = %d",SID_LOG_COMMENT,n_trees_group_total);
  SID_log("n_trees_subgroups = %d",SID_LOG_COMMENT,n_trees_subgroup_total);

  int  i_tree;
  int *i_forest_subgroups;
  int *i_forest_subgroups_start;
  int *i_forest_groups;
  int *i_forest_list;
  int *i_forest_list_snap;
  int *i_forest_list_size;
  int  n_forest_list=0;
  int  max_link_length=n_search*i_read_step;
  i_forest_subgroups      =(int *)SID_calloc(sizeof(int)*n_trees_subgroup_total);
  i_forest_subgroups_start=(int *)SID_calloc(sizeof(int)*n_trees_subgroup_total);
  i_forest_groups         =(int *)SID_calloc(sizeof(int)*n_trees_group_total);
  i_forest_list           =(int *)SID_calloc(sizeof(int)*n_trees_group_total);
  i_forest_list_snap      =(int *)SID_calloc(sizeof(int)*n_trees_group_total);
  i_forest_list_size      =(int *)SID_calloc(sizeof(int)*n_trees_group_total);
  for(i_tree=0;i_tree<n_trees_subgroup_total;i_tree++)
     i_forest_subgroups[i_tree]=n_trees_subgroup_total+1;
  for(i_tree=0;i_tree<n_trees_group_total;i_tree++)
     i_forest_groups[i_tree]=i_tree;
  for(i_read=i_read_stop;i_read>=i_read_start;i_read-=i_read_step){
     // Open catalog files and read headers
     fp_catalog_info  fp_group_properties;
     fp_catalog_info  fp_subgroup_properties;
     halo_properties_SAGE_info        group_properties;
     halo_properties_SAGE_info        subgroup_properties;
     fopen_catalog(filename_catalog_root,
                   i_read,
                   READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES,
                   &fp_group_properties);
     fopen_catalog(filename_catalog_root,
                   i_read,
                   READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES,
                   &fp_subgroup_properties);

     int n_trees_subgroup_i;
     int n_trees_group_i;
     sprintf(filename_dir_horizontal,      "%s/horizontal",filename_root);
     sprintf(filename_dir_horizontal_trees,"%s/trees",     filename_dir_horizontal);
     sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_dir_horizontal_trees,i_read);
     SID_log("Reading {%s}...",SID_LOG_OPEN,filename_in);
     SID_fopen(filename_in,"r",&fp_in);
     SID_fread_all(&n_step_in,         sizeof(int),1,&fp_in);
     SID_fread_all(&n_search_in,       sizeof(int),1,&fp_in);
     SID_fread_all(&n_groups,          sizeof(int),1,&fp_in);
     SID_fread_all(&n_subgroups,       sizeof(int),1,&fp_in);
     SID_fread_all(&n_groups_max_in,   sizeof(int),1,&fp_in);
     SID_fread_all(&n_subgroups_max_in,sizeof(int),1,&fp_in);
     SID_fread_all(&n_trees_subgroup_i,sizeof(int),1,&fp_in);
     SID_fread_all(&n_trees_group_i,   sizeof(int),1,&fp_in);
     SID_log("(%d groups, %d subgroups)...",SID_LOG_CONTINUE,n_groups,n_subgroups);
     int i_group;
     int i_subgroup;
     int j_subgroup;
     for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
       // Read group properties
       fread_catalog_file(&fp_group_properties,&group_properties,NULL,NULL,i_group);

       int group_id;
       int group_type;
       int group_descendant_id;
       int group_tree_id;
       int group_file_offset;
       int group_file_index;
       int n_subgroups_group;
       SID_fread(&(group_id),           sizeof(int),1,&fp_in);
       SID_fread(&(group_type),         sizeof(int),1,&fp_in);
       SID_fread(&(group_descendant_id),sizeof(int),1,&fp_in);
       SID_fread(&(group_tree_id),      sizeof(int),1,&fp_in);
       SID_fread(&(group_file_offset),  sizeof(int),1,&fp_in);
       SID_fread(&(group_file_index),   sizeof(int),1,&fp_in);
       SID_fread(&(n_subgroups_group),  sizeof(int),1,&fp_in);
       int subgroup_tree_id_min=n_trees_subgroup_total+1;
       int subgroup_tree_id_max=0;
       for(j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
         // Read group properties
         fread_catalog_file(&fp_subgroup_properties,&subgroup_properties,NULL,NULL,i_subgroup);

         int subgroup_id;
         int subgroup_type;
         int subgroup_descendant_id;
         int subgroup_tree_id;
         int subgroup_file_offset;
         int subgroup_file_index;
         SID_fread(&(subgroup_id),           sizeof(int),1,&fp_in);
         SID_fread(&(subgroup_type),         sizeof(int),1,&fp_in);
         SID_fread(&(subgroup_descendant_id),sizeof(int),1,&fp_in);
         SID_fread(&(subgroup_tree_id),      sizeof(int),1,&fp_in);
         SID_fread(&(subgroup_file_offset),  sizeof(int),1,&fp_in);
         SID_fread(&(subgroup_file_index),   sizeof(int),1,&fp_in);
         if(i_forest_subgroups[subgroup_tree_id]>n_trees_subgroup_total){
            i_forest_subgroups_start[subgroup_tree_id]=i_read;
            i_forest_subgroups[subgroup_tree_id]      =group_tree_id;
         }
         subgroup_tree_id_min=MIN(subgroup_tree_id_min,i_forest_subgroups[subgroup_tree_id]);
         subgroup_tree_id_max=MAX(subgroup_tree_id_min,i_forest_subgroups[subgroup_tree_id]);
         if(i_tree_extract_group==group_tree_id){
            int flag_found=FALSE;
            for(i_tree=0;i_tree<n_forest_list && !flag_found;i_tree++){
               if(i_forest_list[i_tree]==i_forest_subgroups[subgroup_tree_id])
                  flag_found=TRUE;
            }
            //if(!flag_found && i_forest_subgroups_start[subgroup_tree_id]<(i_read+max_link_length)){
            if(!flag_found){
               i_forest_list[n_forest_list]     =i_forest_subgroups[subgroup_tree_id];
               i_forest_list_snap[n_forest_list]=i_read;
               i_forest_list_size[n_forest_list]=n_subgroups_group;
               n_forest_list++;
            }
         }
         if(i_tree_extract_subgroup==subgroup_tree_id)
           fprintf(fp_out_subgroups,"%6d %6d %6d %6d %6d %6d\n",i_read,group_tree_id,i_group,group_id,i_subgroup,subgroup_id);
       }
       if(subgroup_tree_id_min==(n_trees_subgroup_total+1)){
          subgroup_tree_id_min=-1;
          subgroup_tree_id_max=-1;
       }
       if(i_tree_extract_group==group_tree_id)
          fprintf(fp_out_groups,"%4d %8d %8d %8d %8d %8d %d %8d %8d %le %le %le\n",i_read,i_group,group_file_index,group_id,n_subgroups_group,group_type,check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_STRAYED)||check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_RETURNED)||check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_EXCHANGED),subgroup_tree_id_min,subgroup_tree_id_max,group_properties.pos[0],group_properties.pos[1],group_properties.pos[2]);
     }
     SID_fclose(&fp_in);
     fclose_catalog(&fp_group_properties);
     fclose_catalog(&fp_subgroup_properties);
     SID_log("Done.",SID_LOG_CLOSE);
  }

  // Write forest
  FILE   *fp_out_list;
  size_t *i_forest_list_index=NULL;
  SID_log("Writing %d group trees IDs in this one's forest...",SID_LOG_OPEN,n_forest_list); 
  fp_out_list=fopen(filename_out_list,"w");
  merge_sort(i_forest_list,(size_t)n_forest_list,&i_forest_list_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  for(i_tree=0;i_tree<n_forest_list;i_tree++){
     size_t i_tree_i=i_forest_list_index[i_tree];
     fprintf(fp_out_list,"%9d %9d %9d\n",i_forest_list[i_tree_i],i_forest_list_snap[i_tree_i],i_forest_list_size[i_tree_i]);
  }
  fclose(fp_out_list);
  SID_log("Done.",SID_LOG_CLOSE); 

  // Clean-up
  if(fp_out_groups!=NULL)    fclose(fp_out_groups);
  if(fp_out_subgroups!=NULL) fclose(fp_out_subgroups);
  SID_log("Done.",SID_LOG_CLOSE); 
  SID_exit(ERROR_NONE);
}

