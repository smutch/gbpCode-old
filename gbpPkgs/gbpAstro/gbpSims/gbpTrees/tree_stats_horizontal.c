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
  int  snap_start,snap_stop,snap_step;
  int  tree_type;

  strcpy(filename_root,argv[1]);
 
  SID_log("Computing horizontal tree stats...",SID_LOG_OPEN|SID_LOG_TIMER);
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

  // Allocate arrays
  size_t *n_halos_tree_groups;
  size_t *n_halos_tree_subgroups;
  n_halos_tree_groups   =(size_t *)SID_calloc(sizeof(size_t)*n_trees_group_total);
  n_halos_tree_subgroups=(size_t *)SID_calloc(sizeof(size_t)*n_trees_subgroup_total);
  for(i_read=i_read_stop;i_read>=i_read_start;i_read-=i_read_step){
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
     int i_group;
     int i_subgroup;
     int j_subgroup;
     for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
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
       if(group_tree_id<0 || group_tree_id>=n_trees_group_total)
          SID_trap_error("group_tree_id (%d) exceeds range (0->%d)",group_tree_id,n_trees_group_total);
       n_halos_tree_groups[group_tree_id]++;
//       fprintf(stdout,"%7d %7d %7d %7d %7d\n",i_read,i_group,n_subgroups_group,group_id,group_tree_id);
       for(j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
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
         if(subgroup_tree_id<0 || subgroup_tree_id>=n_trees_subgroup_total)
            SID_trap_error("subgroup_tree_id (%d) exceeds range (0->%d)",subgroup_tree_id,n_trees_subgroup_total);
         n_halos_tree_subgroups[subgroup_tree_id]++;
//         fprintf(stdout,"%7d %7d %7d %7d %7d %7d %7d %7d\n",
//                        i_read,i_subgroup,i_group,j_subgroup,
//                        group_id,group_tree_id,subgroup_id,subgroup_tree_id);
       }
     }
     SID_fclose(&fp_in);
     SID_log("Done.",SID_LOG_CLOSE);
  }

  // Write results
  int i_tree;
  char  filename_out[256];
  FILE *fp_out;
  sprintf(filename_out,"groups.out");
  SID_log("Writing {%s}...",SID_LOG_OPEN,filename_out);
  fp_out=fopen(filename_out,"w");
  for(i_tree=0;i_tree<n_trees_group_total;i_tree++){
     fprintf(fp_out,"%7d %7lld\n",i_tree,n_halos_tree_groups[i_tree]);
  }
  fclose(fp_out);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_log("Writing {%s}...",SID_LOG_OPEN,filename_out);
  sprintf(filename_out,"subgroups.out");
  fp_out=fopen(filename_out,"w");
  for(i_tree=0;i_tree<n_trees_subgroup_total;i_tree++){
     fprintf(fp_out,"%7d %7lld\n",i_tree,n_halos_tree_subgroups[i_tree]);
  }
  fclose(fp_out);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_log("Done.",SID_LOG_CLOSE); 

  // Clean-up
  SID_free(SID_FARG n_halos_tree_groups);
  SID_free(SID_FARG n_halos_tree_subgroups);
  SID_exit(ERROR_NONE);
}

