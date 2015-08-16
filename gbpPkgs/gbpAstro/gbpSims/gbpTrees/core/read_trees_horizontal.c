#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void read_trees_horizontal(void **groups_in,   int *n_groups_in,
                           void **subgroups_in,int *n_subgroups_in,
                           int   *n_subgroups_group,
                           int   *n_trees_subgroup_in,
                           int   *n_trees_group_in,
                           int    i_file, // tree snapshot index
                           int    i_read, // actual snapshot index
                           int    j_file,
                           int    n_wrap,
                           char  *filename_output_dir,
                           int    mode){
   char   filename_output_dir_horizontal[MAX_FILENAME_LENGTH];
   char   filename_output_dir_horizontal_trees[MAX_FILENAME_LENGTH];
   char   filename_in[MAX_FILENAME_LENGTH];
   char   filename_output_file_root[MAX_FILENAME_LENGTH];
   SID_fp fp_in;
   SID_log("Reading horizontal trees for snapshot #%03d...",SID_LOG_OPEN,i_read);

   // Parse the mode and set things up accordingly
   tree_horizontal_extended_info            *groups_extended;
   tree_horizontal_extended_info            *subgroups_extended;
   tree_horizontal_extended_info           **subgroups_extended_all;
   if(check_mode_for_flag(mode,TREE_HORIZONTAL_STORE_EXTENDED) &&
      check_mode_for_flag(mode,TREE_HORIZONTAL_STORE_GHOSTS))
      SID_trap_error("Too many storage modes chosen in read_trees_horizontal.",ERROR_LOGIC);
   int flag_read_extended =check_mode_for_flag(mode,TREE_HORIZONTAL_READ_EXTENDED);
   int flag_store_extended=check_mode_for_flag(mode,TREE_HORIZONTAL_STORE_EXTENDED);
   int flag_store_ghosts  =check_mode_for_flag(mode,TREE_HORIZONTAL_STORE_GHOSTS);
   if(flag_store_ghosts)
      SID_trap_error("Ghost processing not supported in read_trees_horizontal().",ERROR_LOGIC);

   // Set filename and open file
   strcpy(filename_output_file_root,filename_output_dir);
   strip_path(filename_output_file_root);
   sprintf(filename_output_dir_horizontal,      "%s/horizontal",filename_output_dir);
   sprintf(filename_output_dir_horizontal_trees,"%s/trees",     filename_output_dir_horizontal);
   if(flag_read_extended)
      sprintf(filename_in,"%s/horizontal_trees_tmp_%03d.dat",filename_output_dir_horizontal_trees,i_read);
   else
      sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_output_dir_horizontal_trees,i_read);
   SID_fopen(filename_in,"r",&fp_in);

   // Read header
   int n_step_in;
   int n_search_in;
   int n_groups_max_in;
   int n_subgroups_max_in;
   SID_fread_all(&n_step_in,         sizeof(int),1,&fp_in);
   SID_fread_all(&n_search_in,       sizeof(int),1,&fp_in);
   SID_fread_all(n_groups_in,        sizeof(int),1,&fp_in);
   SID_fread_all(n_subgroups_in,     sizeof(int),1,&fp_in);
   SID_fread_all(&n_groups_max_in,   sizeof(int),1,&fp_in);
   SID_fread_all(&n_subgroups_max_in,sizeof(int),1,&fp_in);
   SID_fread_all(n_trees_subgroup_in,sizeof(int),1,&fp_in);
   SID_fread_all(n_trees_group_in,   sizeof(int),1,&fp_in);
   if(flag_store_extended){
      groups_extended       =(tree_horizontal_extended_info  *)groups_in[i_file%n_wrap];
      subgroups_extended    =(tree_horizontal_extended_info  *)subgroups_in[i_file%n_wrap];
      subgroups_extended_all=(tree_horizontal_extended_info **)subgroups_in;
   }

   // Perform read
   int i_group;
   int i_subgroup;
   for(i_group=0,i_subgroup=0;i_group<(*n_groups_in);i_group++){

      // Read groups
      int group_id;
      int group_type;
      int group_descendant_id;
      int group_tree_id;
      int group_file_offset;
      int group_index;
      SID_fread_all(&group_id,            sizeof(int),1,&fp_in);
      SID_fread_all(&group_type,          sizeof(int),1,&fp_in);
      SID_fread_all(&group_descendant_id, sizeof(int),1,&fp_in);
      SID_fread_all(&group_tree_id,       sizeof(int),1,&fp_in);
      SID_fread_all(&group_file_offset,   sizeof(int),1,&fp_in);
      SID_fread_all(&group_index,         sizeof(int),1,&fp_in);
      if(!check_validity_of_tree_case_flag(group_type)) 
         SID_trap_error("Invalid match type (%d) for i_group=%d",ERROR_LOGIC,
                        group_type,i_group);
      if(flag_store_extended){
         groups_extended[i_group].id           =group_id;
         groups_extended[i_group].type         =group_type;
         groups_extended[i_group].descendant_id=group_descendant_id;
         groups_extended[i_group].tree_id      =group_tree_id;
         groups_extended[i_group].file_offset  =group_file_offset;
         groups_extended[i_group].index        =group_index;
      }
      if(flag_read_extended){
         int   group_n_particles;
         int   group_n_particles_parent;
         int   group_n_particles_desc;
         int   group_n_particles_proj;
         float group_score_desc;
         float group_score_prog;
         int   group_snap_bridge;
         int   group_file_bridge;
         int   group_index_bridge;
         int   group_id_bridge;
         SID_fread(&group_n_particles,       sizeof(int),  1,&fp_in);
         SID_fread(&group_n_particles_parent,sizeof(int),  1,&fp_in);
         SID_fread(&group_n_particles_desc,  sizeof(int),  1,&fp_in);
         SID_fread(&group_n_particles_proj,  sizeof(int),  1,&fp_in);
         SID_fread(&group_score_desc,        sizeof(float),1,&fp_in);
         SID_fread(&group_score_prog,        sizeof(float),1,&fp_in);
         SID_fread(&group_snap_bridge,       sizeof(int),  1,&fp_in);
         SID_fread(&group_file_bridge,       sizeof(int),  1,&fp_in);
         SID_fread(&group_index_bridge,      sizeof(int),  1,&fp_in);
         SID_fread(&group_id_bridge,         sizeof(int),  1,&fp_in);
         if(flag_store_extended){
            groups_extended[i_group].n_particles       =group_n_particles;
            groups_extended[i_group].n_particles_parent=group_n_particles_parent;
            groups_extended[i_group].n_particles_desc  =group_n_particles_desc;
            groups_extended[i_group].n_particles_proj  =group_n_particles_proj;
            groups_extended[i_group].score_desc        =group_score_desc;
            groups_extended[i_group].score_prog        =group_score_prog;
            groups_extended[i_group].snap_bridge       =group_snap_bridge;
            groups_extended[i_group].file_bridge       =group_file_bridge;
            groups_extended[i_group].index_bridge      =group_index_bridge;
            groups_extended[i_group].id_bridge         =group_id_bridge;
         }
      }
      SID_fread_all(&(n_subgroups_group[i_group]),sizeof(int),1,&fp_in);

      int j_subgroup;
      for(j_subgroup=0;j_subgroup<n_subgroups_group[i_group];i_subgroup++,j_subgroup++){

         // Read subgroups
         int subgroup_id;
         int subgroup_type;
         int subgroup_descendant_id;
         int subgroup_tree_id;
         int subgroup_file_offset;
         int subgroup_index;
         SID_fread_all(&subgroup_id,           sizeof(int),1,&fp_in);
         SID_fread_all(&subgroup_type,         sizeof(int),1,&fp_in);
         SID_fread_all(&subgroup_descendant_id,sizeof(int),1,&fp_in);
         SID_fread_all(&subgroup_tree_id,      sizeof(int),1,&fp_in);
         SID_fread_all(&subgroup_file_offset,  sizeof(int),1,&fp_in);
         SID_fread_all(&subgroup_index,        sizeof(int),1,&fp_in);
         if(!check_validity_of_tree_case_flag(subgroup_type)) 
            SID_trap_error("Invalid match type (%d) for i_group,j_subgroup,i_subgroup=%d,%d,%d",ERROR_LOGIC,
                           subgroup_type,i_group,j_subgroup,i_subgroup);
         if(flag_store_extended){
            subgroups_extended[i_subgroup].id           =subgroup_id;
            subgroups_extended[i_subgroup].type         =subgroup_type;
            subgroups_extended[i_subgroup].descendant_id=subgroup_descendant_id;
            subgroups_extended[i_subgroup].tree_id      =subgroup_tree_id;
            subgroups_extended[i_subgroup].file_offset  =subgroup_file_offset;
            subgroups_extended[i_subgroup].index        =subgroup_index;
         }
         if(flag_read_extended){
            int   subgroup_n_particles;
            int   subgroup_n_particles_parent;
            int   subgroup_n_particles_desc;
            int   subgroup_n_particles_proj;
            float subgroup_score_desc;
            float subgroup_score_prog;
            int   subgroup_snap_bridge;
            int   subgroup_file_bridge;
            int   subgroup_index_bridge;
            int   subgroup_id_bridge;
            SID_fread(&subgroup_n_particles,       sizeof(int),  1,&fp_in);
            SID_fread(&subgroup_n_particles_parent,sizeof(int),  1,&fp_in);
            SID_fread(&subgroup_n_particles_desc,  sizeof(int),  1,&fp_in);
            SID_fread(&subgroup_n_particles_proj,  sizeof(int),  1,&fp_in);
            SID_fread(&subgroup_score_desc,        sizeof(float),1,&fp_in);
            SID_fread(&subgroup_score_prog,        sizeof(float),1,&fp_in);
            SID_fread(&subgroup_snap_bridge,       sizeof(int),  1,&fp_in);
            SID_fread(&subgroup_file_bridge,       sizeof(int),  1,&fp_in);
            SID_fread(&subgroup_index_bridge,      sizeof(int),  1,&fp_in);
            SID_fread(&subgroup_id_bridge,         sizeof(int),  1,&fp_in);
            if(flag_store_extended){
               subgroups_extended[i_subgroup].n_particles       =subgroup_n_particles;
               subgroups_extended[i_subgroup].n_particles_parent=subgroup_n_particles_parent;
               subgroups_extended[i_subgroup].n_particles_desc  =subgroup_n_particles_desc;
               subgroups_extended[i_subgroup].n_particles_proj  =subgroup_n_particles_proj;
               subgroups_extended[i_subgroup].score_desc        =subgroup_score_desc;
               subgroups_extended[i_subgroup].score_prog        =subgroup_score_prog;
               subgroups_extended[i_subgroup].snap_bridge       =subgroup_snap_bridge;
               subgroups_extended[i_subgroup].file_bridge       =subgroup_file_bridge;
               subgroups_extended[i_subgroup].index_bridge      =subgroup_index_bridge;
               subgroups_extended[i_subgroup].id_bridge         =subgroup_id_bridge;
            }
         }
      }
   }
   if(flag_store_extended){
      for(;i_group<n_groups_max_in;i_group++)
         groups_extended[i_group].type=TREE_CASE_INVALID;
      for(;i_subgroup<n_subgroups_max_in;i_subgroup++)
         subgroups_extended[i_subgroup].type=TREE_CASE_INVALID;
   }
   SID_fclose(&fp_in);
   SID_log("Done.",SID_LOG_CLOSE);
}

