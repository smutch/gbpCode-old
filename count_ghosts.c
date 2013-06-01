#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees.h>

void count_ghosts(int  *n_groups_in,    int *n_group_ghosts,
                  int  *n_subgroups_in, int *n_subgroup_ghosts,
                  char *filename_output_dir,
                  int   i_file,int i_read){
   SID_fp fp_in;
   char   filename_output_dir_horizontal[MAX_FILENAME_LENGTH];
   char   filename_output_dir_horizontal_trees[MAX_FILENAME_LENGTH];
   char   filename_in[MAX_FILENAME_LENGTH];
   char   filename_output_file_root[MAX_FILENAME_LENGTH];

   // Set filename and open file
   strcpy(filename_output_file_root,filename_output_dir);
   strip_path(filename_output_file_root);
   sprintf(filename_output_dir_horizontal,      "%s/horizontal",filename_output_dir);
   sprintf(filename_output_dir_horizontal_trees,"%s/trees",     filename_output_dir_horizontal);
   sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_output_dir_horizontal_trees,i_read);
   SID_fopen(filename_in,"r",&fp_in);

   // Read header
   int n_progenitors_max_in;
   int n_trees_subgroup_in;
   int n_trees_group_in;
   int n_step_in;
   int n_search_in;
   int n_groups_max_in;
   int n_subgroups_max_in;
   SID_fread_all(&n_step_in,          sizeof(int),1,&fp_in);
   SID_fread_all(&n_search_in,        sizeof(int),1,&fp_in);
   SID_fread_all(n_groups_in,         sizeof(int),1,&fp_in);
   SID_fread_all(n_subgroups_in,      sizeof(int),1,&fp_in);
   SID_fread_all(&n_groups_max_in,    sizeof(int),1,&fp_in);
   SID_fread_all(&n_subgroups_max_in, sizeof(int),1,&fp_in);
   SID_fread_all(&n_trees_subgroup_in,sizeof(int),1,&fp_in);
   SID_fread_all(&n_trees_group_in,   sizeof(int),1,&fp_in);

   // Perform read
   int i_group;
   int i_subgroup;
   int n_group_ghosts_count   =0;
   int n_subgroup_ghosts_count=0;
   int n_subgroup_strays_count=0;
   for(i_group=0,i_subgroup=0;i_group<(*n_groups_in);i_group++){

      // Read groups
      int group_id;
      int group_type;
      int group_descendant_id;
      int group_tree_id;
      int group_file_offset;
      int group_index;
      int n_subgroups_group;
      SID_fread_all(&group_id,           sizeof(int),1,&fp_in);
      SID_fread_all(&group_type,         sizeof(int),1,&fp_in);
      SID_fread_all(&group_descendant_id,sizeof(int),1,&fp_in);
      SID_fread_all(&group_tree_id,      sizeof(int),1,&fp_in);
      SID_fread_all(&group_file_offset,  sizeof(int),1,&fp_in);
      SID_fread_all(&group_index,        sizeof(int),1,&fp_in);
      SID_fread_all(&n_subgroups_group,  sizeof(int),1,&fp_in);

      // Count the number of real groups and ghost groups
      if(group_file_offset<0)
         n_group_ghosts[i_file]++;
      else{
         int i_offset;
         for(i_offset=0;i_offset<group_file_offset;i_offset++) 
            n_group_ghosts[i_file+i_offset]++;
         n_group_ghosts_count+=(group_file_offset-1);
      }

      int j_subgroup;
      for(j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){

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

         // Count the number of real subgroups and ghost subgroups
         if(subgroup_file_offset<0){
            n_subgroup_ghosts[i_file]++;
            n_subgroup_strays_count++;
         }
         else if(subgroup_file_offset==0)
            SID_log_warning("subgroup_file_offset==0!",ERROR_NONE);
         else{
            int i_offset;
            for(i_offset=0;i_offset<subgroup_file_offset;i_offset++) 
               n_subgroup_ghosts[i_file+i_offset]++;
            n_subgroup_ghosts_count+=(subgroup_file_offset-1);
         }
      }
   }
   SID_fclose(&fp_in);
fprintf(stderr,"count_ghosts: i_file=%d n_g=%d\n",i_file,n_subgroup_ghosts_count);
   if(i_subgroup!=(*n_subgroups_in))
      SID_trap_error("There's a problem with the substructure counts while reading trees for ghost-counting (ie. %d!=%d)",ERROR_LOGIC,
                     i_subgroup,(*n_subgroups_in));
}

