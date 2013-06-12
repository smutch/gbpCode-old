#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees.h>

void write_trees_horizontal(void  **groups_in, 
                            void  **subgroups_in,
                            int     n_groups,    int n_groups_max,   
                            int     n_subgroups, int n_subgroups_max,
                            int   **n_subgroups_group,
                            int     max_tree_id_subgroup,
                            int     max_tree_id_group,
                            int     i_write,
                            int     j_write,
                            int     l_write,
                            int     n_step,
                            int     n_search,
                            int     n_wrap,
                            int     i_file_start,
                            char   *filename_cat_root_in,
                            char   *filename_output_dir,
                            double *a_list,
                            cosmo_info **cosmo,
                            int     n_k_match,
                            int     flag_init_write,
                            int     mode){
   char        filename_output_dir_horizontal[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_trees[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_cases[MAX_FILENAME_LENGTH];
   char        filename_output_file_root[MAX_FILENAME_LENGTH];
   char        filename_matches_out[MAX_FILENAME_LENGTH];
   char        filename_groups[MAX_FILENAME_LENGTH];
   char        filename_subgroups[MAX_FILENAME_LENGTH];
   char        filename_log[MAX_FILENAME_LENGTH];
   char        filename_matching_out[MAX_FILENAME_LENGTH];
   char        filename_mergers_out[MAX_FILENAME_LENGTH];
   char        filename_strayed_out[MAX_FILENAME_LENGTH];
   char        filename_sputtered_out[MAX_FILENAME_LENGTH];
   char        filename_dropped_out[MAX_FILENAME_LENGTH];
   char        filename_bridged_out[MAX_FILENAME_LENGTH];
   char        filename_emerged_out[MAX_FILENAME_LENGTH];
   char        filename_fragmented_out[MAX_FILENAME_LENGTH];
   FILE       *fp_matching_out;
   FILE       *fp_mergers_out;
   FILE       *fp_strayed_out;
   FILE       *fp_sputtered_out;
   FILE       *fp_dropped_out;
   FILE       *fp_bridged_out;
   FILE       *fp_emerged_out;
   FILE       *fp_fragmented_out;
   int         i_halo;
   char        group_text_prefix[5];
   SID_fp      fp_matches_out;
   FILE       *fp;
   int         n_halos;
   int         n_emerged;
   int         n_fragmented_lost;
   int         n_fragmented_lost_main;
   int         n_fragmented_returned;
   int         n_fragmented_returned_main;
   int         n_fragmented_exchanged;
   int         n_fragmented_exchanged_main;
   int         j_halo;
   int         k_halo;
   int         file_offset;

   int         file_index;
   int         i_subgroup;
   int         j_subgroup;
   int         i_group;
   int         j_group;
   int         i_k_match;
   int         j_k_match;
   int         i_column;
   int         n_particles_emerged;
   int         n_particles_emerged_main;
   int         n_particles_fragmented_lost;
   int         n_particles_fragmented_lost_main;
   int         n_particles_fragmented_returned;
   int         n_particles_fragmented_returned_main;
   int         n_particles_fragmented_exchanged;
   int         n_particles_fragmented_exchanged_main;
   int         n_p_largest;
   int         n_p_largest_main;
   int         n_p_largest_i;
   int         n_p_largest_index;
   double      dt_descendant;
   double      dt_progenitor;
   char       *line=NULL;
   int         line_length=0;
   tree_horizontal_stats_info  stats;
   int                         desc_id;

   SID_log("Writing results for snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,j_write);

   // Interpret the mode
   int flag_check_fragmented;
   int flag_write_allcases;
   int flag_write_nocases;
   int flag_write_extended;
   int flag_write_ghosts;
   flag_check_fragmented=check_mode_for_flag(mode,TREE_HORIZONTAL_WRITE_CHECK_FRAGMENTED);
   flag_write_allcases  =check_mode_for_flag(mode,TREE_HORIZONTAL_WRITE_ALLCASES);
   flag_write_nocases   =check_mode_for_flag(mode,TREE_HORIZONTAL_WRITE_NOCASES);
   flag_write_extended  =check_mode_for_flag(mode,TREE_HORIZONTAL_WRITE_EXTENDED);
   flag_write_ghosts    =check_mode_for_flag(mode,TREE_HORIZONTAL_WRITE_GHOSTS);
   if(flag_write_extended && flag_write_ghosts) SID_trap_error("Incompatible mode flags set in write_trees_horizontal.",ERROR_LOGIC);
   if(flag_write_ghosts){
      flag_write_nocases =TRUE;
      flag_write_allcases=FALSE;
   }

   // Check for fragemented halos
   if(flag_check_fragmented){
      check_for_fragmented_halos(0,(tree_horizontal_info **)subgroups_in,n_subgroups,i_write,j_write,l_write,n_wrap);
      check_for_fragmented_halos(1,(tree_horizontal_info **)groups_in,   n_groups,   i_write,j_write,l_write,n_wrap);
   }

   // Set filenames and open file for writing
   sprintf(filename_output_dir_horizontal,      "%s/horizontal",   filename_output_dir);
   sprintf(filename_output_dir_horizontal_trees,"%s/trees",        filename_output_dir_horizontal);
   sprintf(filename_output_dir_horizontal_cases,"%s/special_cases",filename_output_dir_horizontal);
   mkdir(filename_output_dir,                 02755);
   mkdir(filename_output_dir_horizontal,      02755);
   mkdir(filename_output_dir_horizontal_trees,02755);
   mkdir(filename_output_dir_horizontal_cases,02755);
   strcpy(filename_output_file_root,filename_output_dir);
   strip_path(filename_output_file_root);
   if(flag_write_extended)
      sprintf(filename_matches_out,"%s/horizontal_trees_tmp_%03d.dat",filename_output_dir_horizontal_trees,j_write);
   else if(flag_write_ghosts)
      sprintf(filename_matches_out,"%s/horizontal_trees_ghosts_%03d.dat",filename_output_dir_horizontal_trees,j_write);
   else
      sprintf(filename_matches_out,"%s/horizontal_trees_%03d.dat",filename_output_dir_horizontal_trees,j_write);
   SID_fopen(filename_matches_out,"w",&fp_matches_out);

   // Write the header information for the horizontal trees files
   SID_log("Writing tree files...",SID_LOG_OPEN);
   int n_halos_max_header;
   n_halos_max_header=MAX(n_subgroups_max,n_groups_max);
   SID_fwrite(&n_step,              sizeof(int),1,&fp_matches_out);
   SID_fwrite(&n_search,            sizeof(int),1,&fp_matches_out);
   SID_fwrite(&n_groups,            sizeof(int),1,&fp_matches_out);
   SID_fwrite(&n_subgroups,         sizeof(int),1,&fp_matches_out);
   SID_fwrite(&n_groups_max,        sizeof(int),1,&fp_matches_out);
   SID_fwrite(&n_subgroups_max,     sizeof(int),1,&fp_matches_out);
   SID_fwrite(&max_tree_id_subgroup,sizeof(int),1,&fp_matches_out);
   SID_fwrite(&max_tree_id_group,   sizeof(int),1,&fp_matches_out);

   // Write the horizontal tree files
   if(flag_write_ghosts){

      // Create indices
      int j_group;
      int n_subgroups_indexed;
      tree_horizontal_ghost_group_info *groups;
      groups=(tree_horizontal_ghost_group_info *)groups_in[i_write%n_wrap];
      for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
         tree_horizontal_ghost_subgroup_info *current;
         current=groups[i_group].first_substructure;
         while(current!=NULL){
            current->halo_index=i_subgroup++;
            if(current->descendant!=NULL)
               current->file_index=current->descendant->halo_index;
            else
               current->file_index=-1;
            current=current->next_substructure;
         }
      }
      n_subgroups_indexed=i_subgroup;
      if(n_subgroups_indexed!=n_subgroups)
         SID_trap_error("The number of substructures indexed do not match the given substructure counts (ie. %d!=%d).",
                        ERROR_LOGIC,
                        n_subgroups_indexed,n_subgroups);

      // Loop over the groups to perform the write
      int n_subgroups_written=0;
      int n_ghosts_written   =0;
      for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
         int   group_id;
         int   group_type;
         int   group_descendant_id;
         int   group_tree_id;
         int   group_file_offset;
         int   group_file_index;
         group_id           =groups[i_group].id;
         group_type         =groups[i_group].type;
         group_descendant_id=groups[i_group].descendant_id;
         group_tree_id      =groups[i_group].tree_id;
         group_file_index   =groups[i_group].halo_index;

         // Write group information to the horizontal tree files
         SID_fwrite(&group_id,           sizeof(int),1,&fp_matches_out);
         SID_fwrite(&group_type,         sizeof(int),1,&fp_matches_out);
         SID_fwrite(&group_descendant_id,sizeof(int),1,&fp_matches_out);
         SID_fwrite(&group_tree_id,      sizeof(int),1,&fp_matches_out);
         SID_fwrite(&group_file_index,   sizeof(int),1,&fp_matches_out);
         SID_fwrite(&(groups[i_group].n_subgroups),sizeof(int),1,&fp_matches_out);
 
         // Write the substructure information
         tree_horizontal_ghost_subgroup_info *current;
         int n_subgroups_i=0;
         int n_ghosts_i   =0;
         current=groups[i_group].first_substructure;
         while(current!=NULL){
            int subgroup_id;
            int subgroup_type;
            int subgroup_descendant_id;
            int subgroup_tree_id;
            int subgroup_file_offset;
            int subgroup_file_index;
            subgroup_id           =current->id;
            subgroup_type         =current->type;
            subgroup_descendant_id=current->descendant_id;
            subgroup_tree_id      =current->tree_id;
            subgroup_file_index   =current->halo_index;

            // Write subgroup information to the horizontal trees files
            SID_fwrite(&subgroup_id,           sizeof(int),1,&fp_matches_out);
            SID_fwrite(&subgroup_type,         sizeof(int),1,&fp_matches_out);
            SID_fwrite(&subgroup_descendant_id,sizeof(int),1,&fp_matches_out);
            SID_fwrite(&subgroup_tree_id,      sizeof(int),1,&fp_matches_out);
            SID_fwrite(&subgroup_file_index,   sizeof(int),1,&fp_matches_out);            

            // Increment counters and move to the next substructure
            n_subgroups_i++;
            n_ghosts_i+=check_mode_for_flag(current->type,TREE_CASE_GHOST);
            current=current->next_substructure;
         }
         n_subgroups_written+=n_subgroups_i;
         if(n_subgroups_i!=groups[i_group].n_subgroups)
            SID_trap_error("The number of substructures written for i_group=%d does not match the given count (ie. %d!=%d).",
                           ERROR_LOGIC,
                           i_group,n_subgroups_i,groups[i_group].n_subgroups);
         n_ghosts_written+=n_ghosts_i;
      }
      if(n_subgroups_written!=n_subgroups)
         SID_trap_error("The number of substructures written do not match the given substructure counts (ie. %d!=%d; n_ghosts_written=%d n_groups=%d).",
                        ERROR_LOGIC,
                        n_subgroups_written,n_subgroups,n_ghosts_written,n_groups);
   }
   else{
      if(n_groups>0){
        // Loop over the groups
        for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){

           // Gather the information we want to write
           int   group_id;
           int   group_type;
           int   group_descendant_id;
           int   group_tree_id;
           int   group_file_offset;
           int   group_file_index;
           int   group_n_particles;
           int   group_n_particles_parent;
           int   group_n_particles_desc;
           int   group_n_particles_proj;
           int   group_snap_bridge;
           int   group_index_bridge;
           int   group_id_bridge;
           float group_score_desc;
           float group_score_prog;
           if(flag_write_extended){
              tree_horizontal_info **groups;
              groups             =(tree_horizontal_info **)groups_in;
              group_id           =groups[i_write%n_wrap][i_group].id;
              group_type         =groups[i_write%n_wrap][i_group].type;
              group_descendant_id=set_match_id(&(groups[i_write%n_wrap][i_group].descendant));
              group_tree_id      =groups[i_write%n_wrap][i_group].tree_id;
              group_file_index   =set_match_index(&(groups[i_write%n_wrap][i_group].descendant));
              // compute_trees_verticle wants the file offsets to be -ve for the roots, +ve everywhere else
              if(i_write<i_file_start){
                 group_file_offset=set_match_file(&(groups[i_write%n_wrap][i_group].descendant));
                 if(group_file_offset>i_write)
                    group_file_offset-=i_write;
              }
              else
                 group_file_offset=-1;
              group_n_particles       =groups[i_write%n_wrap][i_group].n_particles;
              group_n_particles_parent=groups[i_write%n_wrap][i_group].n_particles_parent;
              group_n_particles_desc  =set_match_n_particles(&(groups[i_write%n_wrap][i_group].descendant));
              if((groups[i_write%n_wrap][i_group].descendant.halo)!=NULL)
                 group_n_particles_proj=set_match_n_particles(&(groups[i_write%n_wrap][i_group].descendant.halo->first_progenitor));
              else
                 group_n_particles_proj=-1;
              group_score_desc        =set_match_score(   &(groups[i_write%n_wrap][i_group].descendant));
              group_score_prog        =set_match_score(   &(groups[i_write%n_wrap][i_group].first_progenitor));
              group_snap_bridge       =set_match_snapshot(&(groups[i_write%n_wrap][i_group].bridge_backmatch));
              group_index_bridge      =set_match_index(   &(groups[i_write%n_wrap][i_group].bridge_backmatch));
              group_id_bridge         =set_match_id(      &(groups[i_write%n_wrap][i_group].bridge_backmatch));
           }
           else{
              tree_horizontal_extended_info **groups;
              groups             =(tree_horizontal_extended_info **)groups_in;
              group_id           =groups[i_write%n_wrap][i_group].id;
              group_type         =groups[i_write%n_wrap][i_group].type;
              group_descendant_id=groups[i_write%n_wrap][i_group].descendant_id;
              group_tree_id      =groups[i_write%n_wrap][i_group].tree_id;
              group_file_offset  =groups[i_write%n_wrap][i_group].file_offset;
              group_file_index   =groups[i_write%n_wrap][i_group].index;
           }

           // Write group information to the horizontal tree files
           SID_fwrite(&group_id,           sizeof(int),1,&fp_matches_out);
           SID_fwrite(&group_type,         sizeof(int),1,&fp_matches_out);
           SID_fwrite(&group_descendant_id,sizeof(int),1,&fp_matches_out);
           SID_fwrite(&group_tree_id,      sizeof(int),1,&fp_matches_out);
           SID_fwrite(&group_file_offset,  sizeof(int),1,&fp_matches_out);
           SID_fwrite(&group_file_index,   sizeof(int),1,&fp_matches_out);
           if(flag_write_extended){
              SID_fwrite(&group_n_particles,       sizeof(int),  1,&fp_matches_out);
              SID_fwrite(&group_n_particles_parent,sizeof(int),  1,&fp_matches_out);
              SID_fwrite(&group_n_particles_desc,  sizeof(int),  1,&fp_matches_out);
              SID_fwrite(&group_n_particles_proj,  sizeof(int),  1,&fp_matches_out);
              SID_fwrite(&group_score_desc,        sizeof(float),1,&fp_matches_out);
              SID_fwrite(&group_score_prog,        sizeof(float),1,&fp_matches_out);
              SID_fwrite(&group_snap_bridge,       sizeof(int),  1,&fp_matches_out);
              SID_fwrite(&group_index_bridge,      sizeof(int),  1,&fp_matches_out);
              SID_fwrite(&group_id_bridge,         sizeof(int),  1,&fp_matches_out);
           }
           SID_fwrite(&(n_subgroups_group[i_write%n_wrap][i_group]),sizeof(int),1,&fp_matches_out);

           // Write the subgroup information to the horizontal tree files
           for(j_subgroup=0;j_subgroup<n_subgroups_group[i_write%n_wrap][i_group];j_subgroup++,i_subgroup++){

              // Gather the information we want to write
              int   subgroup_id;
              int   subgroup_type;
              int   subgroup_descendant_id;
              int   subgroup_tree_id;
              int   subgroup_file_offset;
              int   subgroup_file_index;
              int   subgroup_n_particles;
              int   subgroup_n_particles_parent;
              int   subgroup_n_particles_desc;
              int   subgroup_n_particles_proj;
              int   subgroup_snap_bridge;
              int   subgroup_index_bridge;
              int   subgroup_id_bridge;
              float subgroup_score_desc;
              float subgroup_score_prog;
              if(flag_write_extended){
                 tree_horizontal_info **subgroups;
                 subgroups                  =(tree_horizontal_info **)subgroups_in;
                 subgroup_id                =subgroups[i_write%n_wrap][i_subgroup].id;
                 subgroup_type              =subgroups[i_write%n_wrap][i_subgroup].type;
                 subgroup_descendant_id     =set_match_id(&(subgroups[i_write%n_wrap][i_subgroup].descendant));
                 subgroup_tree_id           =subgroups[i_write%n_wrap][i_subgroup].tree_id;
                 subgroup_file_index        =set_match_index(&(subgroups[i_write%n_wrap][i_subgroup].descendant));
                 // compute_trees_verticle wants the file offsets to be -ve for the roots, +ve everywhere else
                 if(i_write<i_file_start){
                    subgroup_file_offset=set_match_file(&(subgroups[i_write%n_wrap][i_subgroup].descendant));
                    if(subgroup_file_offset>i_write)
                       subgroup_file_offset-=i_write;
                 }
                 else
                    subgroup_file_offset=-1;
                 subgroup_n_particles       =subgroups[i_write%n_wrap][i_subgroup].n_particles;
                 subgroup_n_particles_parent=subgroups[i_write%n_wrap][i_subgroup].n_particles_parent;
                 subgroup_n_particles_desc  =set_match_n_particles(&(subgroups[i_write%n_wrap][i_subgroup].descendant));
                 if((subgroups[i_write%n_wrap][i_subgroup].descendant.halo)!=NULL)
                    subgroup_n_particles_proj=set_match_n_particles(&(subgroups[i_write%n_wrap][i_subgroup].descendant.halo->first_progenitor));
                 else
                    subgroup_n_particles_proj=-1;
                 subgroup_score_desc        =set_match_score(   &(subgroups[i_write%n_wrap][i_subgroup].descendant));
                 subgroup_score_prog        =set_match_score(   &(subgroups[i_write%n_wrap][i_subgroup].first_progenitor));
                 subgroup_snap_bridge       =set_match_snapshot(&(subgroups[i_write%n_wrap][i_subgroup].bridge_backmatch));
                 subgroup_index_bridge      =set_match_index(   &(subgroups[i_write%n_wrap][i_subgroup].bridge_backmatch));
                 subgroup_id_bridge         =set_match_id(      &(subgroups[i_write%n_wrap][i_subgroup].bridge_backmatch));
              }
              else{
                 tree_horizontal_extended_info **subgroups;
                 subgroups             =(tree_horizontal_extended_info **)subgroups_in;
                 subgroup_id           =subgroups[i_write%n_wrap][i_subgroup].id;
                 subgroup_type         =subgroups[i_write%n_wrap][i_subgroup].type;
                 subgroup_descendant_id=subgroups[i_write%n_wrap][i_subgroup].descendant_id;
                 subgroup_tree_id      =subgroups[i_write%n_wrap][i_subgroup].tree_id;
                 subgroup_file_offset  =subgroups[i_write%n_wrap][i_subgroup].file_offset;
                 subgroup_file_index   =subgroups[i_write%n_wrap][i_subgroup].index;
              }

              // Write subgroup information to the horizontal trees files
              SID_fwrite(&subgroup_id,           sizeof(int),1,&fp_matches_out);
              SID_fwrite(&subgroup_type,         sizeof(int),1,&fp_matches_out);
              SID_fwrite(&subgroup_descendant_id,sizeof(int),1,&fp_matches_out);
              SID_fwrite(&subgroup_tree_id,      sizeof(int),1,&fp_matches_out);
              SID_fwrite(&subgroup_file_offset,  sizeof(int),1,&fp_matches_out);
              SID_fwrite(&subgroup_file_index,   sizeof(int),1,&fp_matches_out);
              if(flag_write_extended){
                 SID_fwrite(&subgroup_n_particles,       sizeof(int),  1,&fp_matches_out);
                 SID_fwrite(&subgroup_n_particles_parent,sizeof(int),  1,&fp_matches_out);
                 SID_fwrite(&subgroup_n_particles_desc,  sizeof(int),  1,&fp_matches_out);
                 SID_fwrite(&subgroup_n_particles_proj,  sizeof(int),  1,&fp_matches_out);
                 SID_fwrite(&subgroup_score_desc,        sizeof(float),1,&fp_matches_out);
                 SID_fwrite(&subgroup_score_prog,        sizeof(float),1,&fp_matches_out);
                 SID_fwrite(&subgroup_snap_bridge,       sizeof(int),  1,&fp_matches_out);
                 SID_fwrite(&subgroup_index_bridge,      sizeof(int),  1,&fp_matches_out);
                 SID_fwrite(&subgroup_id_bridge,         sizeof(int),  1,&fp_matches_out);
              }
           }
        }
      }
   }
   SID_fclose(&fp_matches_out);
   SID_log("Done.",SID_LOG_CLOSE);

   // Write statistics to ascii files
   SID_log("Writing statistics etc...",SID_LOG_OPEN|SID_LOG_TIMER);
   if(!flag_write_ghosts){
      for(i_k_match=0;i_k_match<n_k_match;i_k_match++){
         // Initialize a bunch of stuff depending on whether
         //   we are processing groups or subgroups
         switch(i_k_match){
            case 0:
               compute_trees_horizontal_stats((void *)(subgroups_in[i_write%n_wrap]),n_subgroups,n_subgroups_max,&stats,flag_write_allcases);
               break;
            case 1:
               compute_trees_horizontal_stats((void *)(groups_in[i_write%n_wrap]),   n_groups,   n_groups_max,   &stats,flag_write_allcases);
               break;
         }

         // Write snapshot summary statistics
         if(flag_write_extended)
            sprintf(filename_log,"%s/log_preprop.txt",filename_output_dir_horizontal);
         else
            sprintf(filename_log,"%s/log.txt",filename_output_dir_horizontal);
         if(flag_init_write && i_k_match==0)
            write_trees_horizontal_log_file(filename_log,l_write,j_write,i_k_match,n_k_match,&stats,a_list,cosmo,TRUE);
         else
            write_trees_horizontal_log_file(filename_log,l_write,j_write,i_k_match,n_k_match,&stats,a_list,cosmo,FALSE);
      }
   }

   if(flag_write_allcases){
      tree_horizontal_info  *halos;
      tree_horizontal_info **halos_all;
      for(i_k_match=0;i_k_match<n_k_match;i_k_match++){
         // Initialize a bunch of stuff depending on whether
         //   we are processing groups or subgroups
         switch(i_k_match){
            case 0:
               halos    =(tree_horizontal_info  *)subgroups_in[i_write%n_wrap];
               halos_all=(tree_horizontal_info **)subgroups_in;
               n_halos  =n_subgroups;
               sprintf(group_text_prefix,"sub");
               break;
            case 1:
               halos    =(tree_horizontal_info  *)groups_in[i_write%n_wrap];
               halos_all=(tree_horizontal_info **)groups_in;
               n_halos  =n_groups;
               sprintf(group_text_prefix,"");
               break;
         }

         // Write matching and special case information
         sprintf(filename_matching_out,  "%s/%sgroup_progenitors.txt",   filename_output_dir_horizontal_cases,group_text_prefix);
         sprintf(filename_strayed_out,   "%s/%sgroup_strays.txt",        filename_output_dir_horizontal_cases,group_text_prefix);
         sprintf(filename_sputtered_out, "%s/%sgroup_sputters.txt",      filename_output_dir_horizontal_cases,group_text_prefix);
         sprintf(filename_dropped_out,   "%s/%sgroup_drops.txt",         filename_output_dir_horizontal_cases,group_text_prefix);
         sprintf(filename_bridged_out,   "%s/%sgroup_bridges.txt",       filename_output_dir_horizontal_cases,group_text_prefix);
         sprintf(filename_emerged_out,   "%s/%sgroup_emerged.txt",       filename_output_dir_horizontal_cases,group_text_prefix);
         sprintf(filename_mergers_out,   "%s/%sgroup_mergers.txt",       filename_output_dir_horizontal_cases,group_text_prefix);
         sprintf(filename_fragmented_out,"%s/%sgroup_fragmented_new.txt",filename_output_dir_horizontal_cases,group_text_prefix);
         if(flag_init_write){
            fp_matching_out  =fopen(filename_matching_out,  "w");
            fp_strayed_out   =fopen(filename_strayed_out,   "w");
            fp_sputtered_out =fopen(filename_sputtered_out, "w");
            fp_dropped_out   =fopen(filename_dropped_out,   "w");
            fp_bridged_out   =fopen(filename_bridged_out,   "w");
            fp_emerged_out   =fopen(filename_emerged_out,   "w");
         }
         else{
            fp_matching_out  =fopen(filename_matching_out,  "a");
            fp_strayed_out   =fopen(filename_strayed_out,   "a");
            fp_sputtered_out =fopen(filename_sputtered_out, "a");
            fp_dropped_out   =fopen(filename_dropped_out,   "a");
            fp_bridged_out   =fopen(filename_bridged_out,   "a");
            fp_emerged_out   =fopen(filename_emerged_out,   "a");
         }
         if(flag_init_write){
            fp=fp_matching_out;
            i_column=1;
            fprintf(fp,"# (%02d): Tree ID\n",                          i_column++);
            fprintf(fp,"# (%02d): Halo file\n",                        i_column++);
            fprintf(fp,"# (%02d): Halo index\n",                       i_column++);
            fprintf(fp,"# (%02d): Halo ID\n",                          i_column++);
            fprintf(fp,"# (%02d): Progenitor number\n",                i_column++);
            fprintf(fp,"# (%02d): Progenitor file\n",                  i_column++);
            fprintf(fp,"# (%02d): Progenitor index\n",                 i_column++);
            fprintf(fp,"# (%02d): Progenitor ID\n",                    i_column++);
            fprintf(fp,"# (%02d): Progenitor match score\n",           i_column++);
            fprintf(fp,"# (%02d): Number of particles in halo\n",      i_column++);
            fprintf(fp,"# (%02d): Number of particles in progenitor\n",i_column++);
            fp=fp_strayed_out;
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                       i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                        i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                         i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n", i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",             i_column++,group_text_prefix);
            fp=fp_sputtered_out;
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                       i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                        i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                         i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                 i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n", i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",i_column++,group_text_prefix);
            fp=fp_dropped_out;
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup->descendant interval [yrs]\n",                  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                              i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup id\n",                                          i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                      i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n",      i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",                  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",     i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's main progenitor\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's descendant\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's main progenitor\n",       i_column++,group_text_prefix);
            fp=fp_bridged_out;
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                                                                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                                                                 i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                                                                  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup id\n",                                                                              i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of back-matched %sgroups identified with this bridge (excludes descendant)\n",       i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of back-matched %sgroups identified as emerging halos\n",                            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of back-matched %sgroups identified as lost      fragmented halos\n",                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of back-matched %sgroups identified as returning fragmented halos\n",                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of back-matched %sgroups identified as exchanged fragmented halos\n",                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n",                                          i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",                                                      i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's main progenitor\n",                                    i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the largest              back-matched %sgroup\n",                    i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the largest              back-matched %sgroup's main progenitor\n",  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the emerged              back-matched %sgroups\n",                   i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the emerged              back-matched %sgroups's main progenitors\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the returned  fragmented back-matched %sgroups\n",                   i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the returned  fragmented back-matched %sgroups's main progenitors\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the exchanged fragmented back-matched %sgroups\n",                   i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the exchanged fragmented back-matched %sgroups's main progenitors\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the lost      fragmented back-matched %sgroups\n",                   i_column++,group_text_prefix);
            fp=fp_emerged_out;
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup->progenitor interval [yrs]\n",                      i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                          i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                                 i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                                  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup id\n",                                              i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup bridge snapshot number\n",                          i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup bridge snapshot index\n",                           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup bridge id\n",                                       i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's bridge\n",                    i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's descendant\n",                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's main progenitor\n",           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n",          i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's bridge\n",             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",                      i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",         i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's main progenitor\n",    i_column++,group_text_prefix);
         }

         // Loop over each halo
         for(i_halo=0;i_halo<n_halos;i_halo++){

            // Compute the time between the halo and it's descendant
            if(halos[i_halo].descendant.halo!=NULL){
               if(halos[i_halo].descendant.halo->file>i_write){ // This isn't true for first snapshot for instance
                 if(l_write<(halos[i_halo].descendant.halo->file-i_write))
                    SID_trap_error("Unreasonable descendant file offset: %d %d %d",ERROR_LOGIC,l_write,halos[i_halo].descendant.halo->file,i_write);
                 dt_descendant=deltat_a(cosmo,a_list[l_write],a_list[l_write-(halos[i_halo].descendant.halo->file-i_write)])/S_PER_YEAR;
               }
               else
                 dt_descendant=0.;
            }
            else
              dt_descendant=0.;

            // Compute the time between the halo and it's progenitor
            if(halos[i_halo].first_progenitor.halo!=NULL){
               if(halos[i_halo].first_progenitor.halo->file<i_write) // This isn't true for first snapshot for instance
                 dt_progenitor=deltat_a(cosmo,a_list[l_write+(i_write-halos[i_halo].first_progenitor.halo->file)],a_list[l_write])/S_PER_YEAR;
               else
                 dt_progenitor=0.;
            }
            else
              dt_progenitor=0.;

            // Print progenitor information
            match_info *current_match;
            current_match=&(halos[i_halo].first_progenitor);
            j_halo=0;
            while(current_match->halo!=NULL){
               fprintf(fp_matching_out,"%3d %7d   %4d %7d %7d   %4d %7d %7d   %10.4f   %7d %7d\n",
                       halos[i_halo].tree_id,
                       j_write,
                       i_halo,
                       halos[i_halo].id,
                       j_halo,
                       set_match_snapshot(current_match),
                       set_match_index(   current_match),
                       set_match_id(      current_match),
                       set_match_score(   current_match),
                       halos[i_halo].n_particles,
                       set_match_n_particles(current_match));
               j_halo++;
               current_match=&(current_match->halo->next_progenitor);
            }

            // Write strays
            fp=fp_strayed_out;
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_STRAYED))
               fprintf(fp,"%10.3le %8d %8d %8d %8d\n",
                  a_list[l_write],
                  j_write,
                  i_halo,
                  halos[i_halo].n_particles_parent,
                  halos[i_halo].n_particles);

            // Write sputtered halos
            fp=fp_sputtered_out;
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_SPUTTERED))
               fprintf(fp,"%10.3le %8d %8d %8d %8d %8d %8d\n",
                  a_list[l_write],
                  j_write,
                  i_halo,
                  set_match_file(&(halos[i_halo].descendant))-i_write,
                  halos[i_halo].n_particles_parent,
                  halos[i_halo].n_particles,
                  set_match_n_particles(&(halos[i_halo].descendant)));

            // Write dropped halos
            fp=fp_dropped_out;
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_DROPPED))
               fprintf(fp,"%10.3le %10.3le %8d %8d %8d %8d %8d %8d %8d %8d %10.3f %10.3f\n",
                       a_list[l_write],
                       dt_descendant,
                       j_write,
                       i_halo,
                       halos[i_halo].id,
                       set_match_file(&(halos[i_halo].descendant))-i_write,
                       halos[i_halo].n_particles_parent,
                       halos[i_halo].n_particles,
                       set_match_n_particles(&(halos[i_halo].descendant)),
                       set_match_n_particles(&(halos[i_halo].first_progenitor)),
                       set_match_score(&(halos[i_halo].descendant)),
                       set_match_score(&(halos[i_halo].first_progenitor)));

            // Write bridged halos
            fp=fp_bridged_out;
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_BRIDGED)){
               // Find the largest back-matched halo
               j_halo=0;
               n_p_largest     =set_match_n_particles(&(halos[i_halo].bridges[j_halo]));
               n_p_largest_main=set_match_n_particles(&(halos_all[set_match_file(&(halos[i_halo].bridges[j_halo]))%n_wrap]
                                                                   [set_match_index(&(halos[i_halo].bridges[j_halo]))].first_progenitor));
               n_p_largest_index=j_halo;
               for(j_halo=1;j_halo<halos[i_halo].n_bridges;j_halo++){
                  n_p_largest_i=set_match_n_particles(&(halos[i_halo].bridges[j_halo]));
                  if(n_p_largest_i>n_p_largest){
                     n_p_largest      =n_p_largest_i;
                     n_p_largest_main =set_match_n_particles(&(halos_all[set_match_file(&(halos[i_halo].bridges[j_halo]))%n_wrap]
                                                                          [set_match_index(&(halos[i_halo].bridges[j_halo]))].first_progenitor));
                     n_p_largest_index=j_halo;
                  }
               }

               // Count the number of emerged/dropped halos and the number of particles
               //   in each set, as well as the number of particles in their main progenitors
               n_emerged                            =0;
               n_fragmented_lost                    =0;
               n_fragmented_lost_main               =0;
               n_fragmented_returned                =0;
               n_fragmented_returned_main           =0;
               n_fragmented_exchanged               =0;
               n_fragmented_exchanged_main          =0;
               n_particles_emerged                  =0;
               n_particles_emerged_main             =0;
               n_particles_fragmented_lost          =0;
               n_particles_fragmented_lost_main     =0;
               n_particles_fragmented_returned      =0;
               n_particles_fragmented_returned_main =0;
               n_particles_fragmented_exchanged     =0;
               n_particles_fragmented_exchanged_main=0;
               for(j_halo=0;j_halo<halos[i_halo].n_bridges;j_halo++){
                  if(check_mode_for_flag(set_match_type(&(halos[i_halo].bridges[j_halo])),TREE_CASE_EMERGED)){
                     n_particles_emerged     +=set_match_n_particles(&(halos[i_halo].bridges[j_halo]));
                     n_particles_emerged_main+=set_match_n_particles(&(halos_all[set_match_file(&(halos[i_halo].bridges[j_halo]))%n_wrap]
                                                                                  [set_match_index(&(halos[i_halo].bridges[j_halo]))].first_progenitor));
                     n_emerged++;
                  }
                  if(check_mode_for_flag(set_match_type(&(halos[i_halo].bridges[j_halo])),TREE_CASE_FRAGMENTED_LOST)){
                     n_particles_fragmented_lost     +=set_match_n_particles(&(halos[i_halo].bridges[j_halo]));
                     n_particles_fragmented_lost_main+=set_match_n_particles(&(halos_all[set_match_file(&(halos[i_halo].bridges[j_halo]))%n_wrap]
                                                                                          [set_match_index(&(halos[i_halo].bridges[j_halo]))].first_progenitor));
                     n_fragmented_lost++;
                  }
                  if(check_mode_for_flag(set_match_type(&(halos[i_halo].bridges[j_halo])),TREE_CASE_FRAGMENTED_RETURNED)){
                     n_particles_fragmented_returned     +=set_match_n_particles(&(halos[i_halo].bridges[j_halo]));
                     n_particles_fragmented_returned_main+=set_match_n_particles(&(halos_all[set_match_file(&(halos[i_halo].bridges[j_halo]))%n_wrap]
                                                                                              [set_match_index(&(halos[i_halo].bridges[j_halo]))].first_progenitor));
                     n_fragmented_returned++;
                  }
                  if(check_mode_for_flag(set_match_type(&(halos[i_halo].bridges[j_halo])),TREE_CASE_FRAGMENTED_EXCHANGED)){
                     n_particles_fragmented_exchanged     +=set_match_n_particles(&(halos[i_halo].bridges[j_halo]));
                     n_particles_fragmented_exchanged_main+=set_match_n_particles(&(halos_all[set_match_file(&(halos[i_halo].bridges[j_halo]))%n_wrap]
                                                                                               [set_match_index(&(halos[i_halo].bridges[j_halo]))].first_progenitor));
                     n_fragmented_exchanged++;
                  }
               }
               fprintf(fp,"%10.3le %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d\n",
                       a_list[l_write],
                       j_write,
                       i_halo,
                       halos[i_halo].id,
                       halos[i_halo].n_bridges,
                       n_emerged,
                       n_fragmented_lost,
                       n_fragmented_returned,
                       n_fragmented_exchanged,
                       halos[i_halo].n_particles_parent,
                       halos[i_halo].n_particles,
                       set_match_n_particles(&(halos[i_halo].first_progenitor)),
                       set_match_n_particles(&(halos[i_halo].bridges[n_p_largest_index])),
                       n_p_largest_main,
                       n_particles_emerged,
                       n_particles_emerged_main,
                       n_particles_fragmented_returned,
                       n_particles_fragmented_returned_main,
                       n_particles_fragmented_exchanged,
                       n_particles_fragmented_exchanged_main,
                       n_particles_fragmented_lost);
            }

            // Write emerged halos
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_EMERGED)){
               int offset;
               if(set_match_file(&(halos[i_halo].descendant))>=0)
                  offset=set_match_file(&(halos[i_halo].descendant))-i_write;
               else
                  offset=-1;
               fp=fp_emerged_out;
               fprintf(fp,"%10.3le %10.3le %3d %8d %8d %8d %8d %8d %8d %10.3f %10.3f %10.3f %8d %8d %8d %8d %8d\n",
                       a_list[l_write],
                       dt_progenitor,
                       offset,
                       j_write,
                       i_halo,
                       halos[i_halo].id,
                       set_match_snapshot(&(halos[i_halo].bridge_backmatch)),
                       set_match_index(   &(halos[i_halo].bridge_backmatch)),
                       set_match_id(      &(halos[i_halo].bridge_backmatch)),
                       set_match_score(   &(halos[i_halo].bridge_backmatch)),
                       set_match_score(   &(halos[i_halo].descendant)),
                       set_match_score(   &(halos[i_halo].first_progenitor)),
                       halos[i_halo].n_particles_parent,
                       set_match_n_particles(&(halos[i_halo].bridge_backmatch)),
                       halos[i_halo].n_particles,
                       set_match_n_particles(&(halos[i_halo].descendant)),
                       set_match_n_particles(&(halos[i_halo].first_progenitor)));
            }
         } // Loop over i_halo
         fclose(fp_matching_out);
         fclose(fp_strayed_out);
         fclose(fp_sputtered_out);
         fclose(fp_dropped_out);
         fclose(fp_bridged_out);
         fclose(fp_emerged_out);
      } // Loop of i_k_match
   } // If flag_write_allcases

   // The following two cases are written for both truth values of flag_write_allcases
   //    because the propagation of fragmented halos can change these things.  Note,
   //    the output final will be in order of expansion factor in these two cases, opposite
   //    from those above.
   if(!flag_write_nocases){
      for(i_k_match=0;i_k_match<n_k_match;i_k_match++){
         switch(i_k_match){
            case 0:
               n_halos=n_subgroups;
               sprintf(group_text_prefix,"sub");
               break;
            case 1:
               n_halos=n_groups;
               sprintf(group_text_prefix,"");
               break;
         }

         // Set filenames and open files
         sprintf(filename_mergers_out,   "%s/%sgroup_mergers.txt",       filename_output_dir_horizontal_cases,group_text_prefix);
         sprintf(filename_fragmented_out,"%s/%sgroup_fragmented_new.txt",filename_output_dir_horizontal_cases,group_text_prefix);
         if(flag_init_write){
            fp_mergers_out   =fopen(filename_mergers_out,   "w");
            fp_fragmented_out=fopen(filename_fragmented_out,"w");
         }
         else{
            fp_mergers_out   =fopen(filename_mergers_out,   "a");
            fp_fragmented_out=fopen(filename_fragmented_out,"a");
         }

         // Write header information
         if(flag_init_write){
            fp=fp_mergers_out;
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                                   i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                                    i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                                     i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup id\n",                                                 i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n",             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",                         i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup merger's main progenitor\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup's matching score to it's descendant\n",                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup's matching score to it's main progenitor\n",           i_column++,group_text_prefix);
            fp=fp_fragmented_out;
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                               i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                                 i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup id\n",                                             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup fragmented type (0=lost,1=returned,2=exchanged)\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                         i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup bridge snapshot\n",                                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup bridge index\n",                                   i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup bridge id\n",                                      i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's descendant\n",               i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n",         i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",                     i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",        i_column++,group_text_prefix);
         }

         // Loop over each halo
         for(i_halo=0;i_halo<n_halos;i_halo++){
            // Organise the data we will write for these last two cases
            int   halo_id;
            int   halo_type;
            int   halo_file_offset;
            int   halo_n_particles;
            int   halo_n_particles_parent;
            int   halo_n_particles_desc;
            int   halo_n_particles_proj;
            float halo_score_desc;
            float halo_score_prog;
            int   halo_snap_bridge;
            int   halo_index_bridge;
            int   halo_id_bridge;
            int   halo_main_progenitor_id;
            int   halo_backmatch_id;
            if(flag_write_extended){
               tree_horizontal_info  *halos;
               tree_horizontal_info **halos_all;
               switch(i_k_match){
                  case 0:
                     halos    =(tree_horizontal_info  *)subgroups_in[i_write%n_wrap];
                     halos_all=(tree_horizontal_info **)subgroups_in;
                     break;
                  case 1:
                     halos    =(tree_horizontal_info  *)groups_in[i_write%n_wrap];
                     halos_all=(tree_horizontal_info **)groups_in;
                     break;
               }
               halo_id                =halos[i_halo].id;
               halo_type              =halos[i_halo].type;
               halo_file_offset       =set_match_file(&(halos[i_halo].descendant))-i_write;
               halo_n_particles       =halos[i_halo].n_particles,
               halo_n_particles_parent=halos[i_halo].n_particles_parent,
               halo_n_particles_desc  =set_match_n_particles(&(halos[i_halo].descendant));
               if((halos[i_halo].descendant.halo)!=NULL)
                  halo_n_particles_proj=set_match_n_particles(&(halos[i_halo].descendant.halo->first_progenitor));
               else
                  halo_n_particles_proj=-1;
               halo_score_desc        =set_match_score(&(halos[i_halo].descendant));
               halo_score_prog        =set_match_score(&(halos[i_halo].first_progenitor));
               halo_snap_bridge       =set_match_file(&(halos[i_halo].bridge_backmatch));
               halo_index_bridge      =set_match_index(&(halos[i_halo].bridge_backmatch));
               halo_id_bridge         =set_match_id(&(halos[i_halo].bridge_backmatch));
               halo_main_progenitor_id=halos[i_halo].main_progenitor_id;
               halo_backmatch_id      =set_match_id(&(halos[i_halo].bridge_backmatch));
            }
            else{
               tree_horizontal_extended_info  *halos;
               tree_horizontal_extended_info **halos_all;
               // Initialize a bunch of stuff depending on whether
               //   we are processing groups or subgroups
               switch(i_k_match){
                  case 0:
                     halos    =(tree_horizontal_extended_info  *)subgroups_in[i_write%n_wrap];
                     halos_all=(tree_horizontal_extended_info **)subgroups_in;
                     n_halos  =n_subgroups;
                     sprintf(group_text_prefix,"sub");
                     break;
                  case 1:
                     halos    =(tree_horizontal_extended_info  *)groups_in[i_write%n_wrap];
                     halos_all=(tree_horizontal_extended_info **)groups_in;
                     n_halos  =n_groups;
                     sprintf(group_text_prefix,"");
                     break;
               }
               halo_id                =halos[i_halo].id;
               halo_type              =halos[i_halo].type;
               halo_file_offset       =halos[i_halo].file_offset;
               halo_n_particles       =halos[i_halo].n_particles;
               halo_n_particles_parent=halos[i_halo].n_particles_parent;
               halo_n_particles_desc  =halos[i_halo].n_particles_desc;
               halo_n_particles_proj  =halos[i_halo].n_particles_proj;
               halo_score_desc        =halos[i_halo].score_desc;
               halo_score_prog        =halos[i_halo].score_prog;
               halo_snap_bridge       =halos[i_halo].snap_bridge;
               halo_index_bridge      =halos[i_halo].index_bridge;
               halo_id_bridge         =halos[i_halo].id_bridge;
               halo_main_progenitor_id=99;
               halo_backmatch_id      =99;
            }

            // Write mergers (don't include fragmented halos in the list)
            if(check_mode_for_flag(halo_type,TREE_CASE_MERGER)){
               if(!(check_mode_for_flag(halo_type,TREE_CASE_FRAGMENTED_LOST)     ||
                    check_mode_for_flag(halo_type,TREE_CASE_FRAGMENTED_RETURNED) ||
                    check_mode_for_flag(halo_type,TREE_CASE_FRAGMENTED_EXCHANGED))){
                  fp=fp_mergers_out;
                  fprintf(fp,"%10.3le %8d %8d %8d %8d %8d %8d %8d %8d %10.3f %10.3f\n",
                          a_list[l_write],
                          j_write,
                          i_halo,
                          halo_id,
                          halo_file_offset,
                          halo_n_particles_parent,
                          halo_n_particles,
                          halo_n_particles_proj,
                          halo_n_particles_desc,
                          halo_score_desc,
                          halo_score_prog);
               }
            }

            // Write fragmented halos
            if(check_mode_for_flag(halo_type,TREE_CASE_FRAGMENTED_LOST)     ||
               check_mode_for_flag(halo_type,TREE_CASE_FRAGMENTED_RETURNED) ||
               check_mode_for_flag(halo_type,TREE_CASE_FRAGMENTED_EXCHANGED)){
               int type_fragmented;
               type_fragmented=-1;
               if(check_mode_for_flag(halo_type,TREE_CASE_FRAGMENTED_LOST)){
                  if(type_fragmented>0)
                     SID_trap_error("Multiple TREE_CASE_FRAGMENTED switches activated in mode (type=%d)",ERROR_LOGIC,halo_type);
                  type_fragmented=0;
               }
               if(check_mode_for_flag(halo_type,TREE_CASE_FRAGMENTED_RETURNED)){
                  if(type_fragmented>0)
                     SID_trap_error("Multiple TREE_CASE_FRAGMENTED switches activated in mode (type=%d)",ERROR_LOGIC,halo_type);
                  type_fragmented=1;
               }
               if(check_mode_for_flag(halo_type,TREE_CASE_FRAGMENTED_EXCHANGED)){
                  if(type_fragmented>0)
                     SID_trap_error("Multiple TREE_CASE_FRAGMENTED switches activated in mode (type=%d)",ERROR_LOGIC,halo_type);
                  type_fragmented=2;
               }
               fp=fp_fragmented_out;
               fprintf(fp,"%10.3le %8d %8d %8d %8d %8d %8d %8d %8d %10.3f %8d %8d %8d\n",
                       a_list[l_write],
                       j_write,
                       i_halo,
                       halo_id,
                       type_fragmented,
                       halo_file_offset,
                       halo_snap_bridge,
                       halo_index_bridge,
                       halo_id_bridge,
                       halo_score_desc,
                       halo_n_particles_parent,
                       halo_n_particles,
                       halo_n_particles_desc);
            }
         } // Loop over i_halo
         fclose(fp_mergers_out);
         fclose(fp_fragmented_out);
      } // Loop over i_k_match
   }
   SID_log("Done.",SID_LOG_CLOSE);
   SID_log("Done.",SID_LOG_CLOSE);
}

