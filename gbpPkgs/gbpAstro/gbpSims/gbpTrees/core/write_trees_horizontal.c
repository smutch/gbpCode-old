#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

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
   char        filename_pointers_out[MAX_FILENAME_LENGTH];
   char        filename_dom_prog_out[MAX_FILENAME_LENGTH];
   char        filename_groups[MAX_FILENAME_LENGTH];
   char        filename_subgroups[MAX_FILENAME_LENGTH];
   char        filename_log[MAX_FILENAME_LENGTH];
   char        filename_matching_out[MAX_FILENAME_LENGTH];
   char        filename_mergers_out[MAX_FILENAME_LENGTH];
   char        filename_strayed_out[MAX_FILENAME_LENGTH];
   char        filename_dropped_out[MAX_FILENAME_LENGTH];
   char        filename_bridged_out[MAX_FILENAME_LENGTH];
   char        filename_emerged_out[MAX_FILENAME_LENGTH];
   char        filename_fragmented_out[MAX_FILENAME_LENGTH];
   FILE       *fp_matching_out;
   FILE       *fp_mergers_out;
   FILE       *fp_strayed_out;
   FILE       *fp_dropped_out;
   FILE       *fp_bridged_out;
   FILE       *fp_emerged_out;
   FILE       *fp_fragmented_out;
   int         i_halo;
   char        group_text_prefix[5];
   SID_fp      fp_trees_out;
   SID_fp      fp_backmatch_ptrs_out;
   SID_fp      fp_bridge_ptrs_out;
   SID_fp      fp_dom_prog_out;
   FILE       *fp;
   int         n_halos;
   int         n_emerged;
   int         n_fragmented_strayed;
   int         n_fragmented_strayed_main;
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
   int         n_particles_fragmented_strayed;
   int         n_particles_fragmented_strayed_main;
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
   int flag_check_fragmented=check_mode_for_flag(mode,TREE_HORIZONTAL_WRITE_CHECK_FRAGMENTED);
   int flag_write_allcases  =check_mode_for_flag(mode,TREE_HORIZONTAL_WRITE_ALLCASES);
   int flag_write_nocases   =check_mode_for_flag(mode,TREE_HORIZONTAL_WRITE_NOCASES);
   int flag_write_extended  =check_mode_for_flag(mode,TREE_HORIZONTAL_WRITE_EXTENDED);
   int flag_write_ghosts    =check_mode_for_flag(mode,TREE_HORIZONTAL_WRITE_GHOSTS);
   if(flag_write_extended && flag_write_ghosts) SID_trap_error("Incompatible mode flags set in write_trees_horizontal.",ERROR_LOGIC);
   if(flag_write_ghosts){
      flag_write_nocases =TRUE;
      flag_write_allcases=FALSE;
      SID_trap_error("Ghost processing not supported in write_trees_horizontal().",ERROR_LOGIC);
   }

   // Check for fragemented halos
   if(flag_check_fragmented){
      check_for_fragmented_halos(0,(tree_horizontal_info **)subgroups_in,n_subgroups,i_write,j_write,l_write,n_wrap);
      check_for_fragmented_halos(1,(tree_horizontal_info **)groups_in,   n_groups,   i_write,j_write,l_write,n_wrap);
   }

   // Set filenames and open files for writing
   int flag_write_pointers=FALSE;
   int flag_write_dom_prog=FALSE;
   sprintf(filename_output_dir_horizontal,      "%s/horizontal",   filename_output_dir);
   sprintf(filename_output_dir_horizontal_trees,"%s/trees",        filename_output_dir_horizontal);
   mkdir(filename_output_dir,                 02755);
   mkdir(filename_output_dir_horizontal,      02755);
   mkdir(filename_output_dir_horizontal_trees,02755);
   strcpy(filename_output_file_root,filename_output_dir);
   strip_path(filename_output_file_root);
   if(flag_write_extended){
      flag_write_pointers=TRUE;
      flag_write_dom_prog=TRUE;
      sprintf(filename_matches_out,"%s/horizontal_trees_tmp_%03d.dat",filename_output_dir_horizontal_trees,j_write);
   }
   else if(flag_write_ghosts)
      sprintf(filename_matches_out,"%s/horizontal_trees_ghosts_%03d.dat",filename_output_dir_horizontal_trees,j_write);
   else{
      sprintf(filename_matches_out,"%s/horizontal_trees_%03d.dat",filename_output_dir_horizontal_trees,j_write);
   }
   SID_fopen(filename_matches_out, "w",&fp_trees_out);
   if(flag_write_pointers){
      sprintf(filename_pointers_out,"%s/horizontal_trees_bridge_backmatch_pointers_%03d.dat",filename_output_dir_horizontal_trees,j_write);
      SID_fopen(filename_pointers_out,"w",&fp_backmatch_ptrs_out);
      sprintf(filename_pointers_out,"%s/horizontal_trees_bridge_forematch_pointers_%03d.dat",filename_output_dir_horizontal_trees,j_write);
      SID_fopen(filename_pointers_out,"w",&fp_bridge_ptrs_out);
   }
   if(flag_write_dom_prog){
      sprintf(filename_dom_prog_out,"%s/horizontal_trees_dominant_progenitor_pointers_%03d.dat",filename_output_dir_horizontal_trees,j_write);
      SID_fopen(filename_dom_prog_out,"w",&fp_dom_prog_out);
   }

   // Write the header information for the horizontal trees files
   SID_log("Writing tree files...",SID_LOG_OPEN);
   int n_halos_max_header;
   n_halos_max_header=MAX(n_subgroups_max,n_groups_max);
   SID_fwrite(&n_step,              sizeof(int),1,&fp_trees_out);
   SID_fwrite(&n_search,            sizeof(int),1,&fp_trees_out);
   SID_fwrite(&n_groups,            sizeof(int),1,&fp_trees_out);
   SID_fwrite(&n_subgroups,         sizeof(int),1,&fp_trees_out);
   SID_fwrite(&n_groups_max,        sizeof(int),1,&fp_trees_out);
   SID_fwrite(&n_subgroups_max,     sizeof(int),1,&fp_trees_out);
   SID_fwrite(&max_tree_id_subgroup,sizeof(int),1,&fp_trees_out);
   SID_fwrite(&max_tree_id_group,   sizeof(int),1,&fp_trees_out);
   if(flag_write_pointers){
      SID_fwrite(&n_step,              sizeof(int),1,&fp_backmatch_ptrs_out);
      SID_fwrite(&n_search,            sizeof(int),1,&fp_backmatch_ptrs_out);
      SID_fwrite(&n_groups,            sizeof(int),1,&fp_backmatch_ptrs_out);
      SID_fwrite(&n_subgroups,         sizeof(int),1,&fp_backmatch_ptrs_out);
      SID_fwrite(&n_groups_max,        sizeof(int),1,&fp_backmatch_ptrs_out);
      SID_fwrite(&n_subgroups_max,     sizeof(int),1,&fp_backmatch_ptrs_out);
      SID_fwrite(&max_tree_id_subgroup,sizeof(int),1,&fp_backmatch_ptrs_out);
      SID_fwrite(&max_tree_id_group,   sizeof(int),1,&fp_backmatch_ptrs_out);
      SID_fwrite(&n_step,              sizeof(int),1,&fp_bridge_ptrs_out);
      SID_fwrite(&n_search,            sizeof(int),1,&fp_bridge_ptrs_out);
      SID_fwrite(&n_groups,            sizeof(int),1,&fp_bridge_ptrs_out);
      SID_fwrite(&n_subgroups,         sizeof(int),1,&fp_bridge_ptrs_out);
      SID_fwrite(&n_groups_max,        sizeof(int),1,&fp_bridge_ptrs_out);
      SID_fwrite(&n_subgroups_max,     sizeof(int),1,&fp_bridge_ptrs_out);
      SID_fwrite(&max_tree_id_subgroup,sizeof(int),1,&fp_bridge_ptrs_out);
      SID_fwrite(&max_tree_id_group,   sizeof(int),1,&fp_bridge_ptrs_out);
   }
   if(flag_write_dom_prog){
      SID_fwrite(&n_step,              sizeof(int),1,&fp_dom_prog_out);
      SID_fwrite(&n_search,            sizeof(int),1,&fp_dom_prog_out);
      SID_fwrite(&n_groups,            sizeof(int),1,&fp_dom_prog_out);
      SID_fwrite(&n_subgroups,         sizeof(int),1,&fp_dom_prog_out);
      SID_fwrite(&n_groups_max,        sizeof(int),1,&fp_dom_prog_out);
      SID_fwrite(&n_subgroups_max,     sizeof(int),1,&fp_dom_prog_out);
      SID_fwrite(&max_tree_id_subgroup,sizeof(int),1,&fp_dom_prog_out);
      SID_fwrite(&max_tree_id_group,   sizeof(int),1,&fp_dom_prog_out);
   }

   // Write the horizontal tree files
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
        int   group_snap_backmatch;
        int   group_file_backmatch;
        int   group_index_backmatch;
        int   group_snap_forematch;
        int   group_file_forematch;
        int   group_index_forematch;
        int   group_id_backmatch;
        float group_score_desc;
        float group_score_prog;
        int   group_first_progenitor_file;
        int   group_first_progenitor_index;
        int   group_next_progenitor_file;
        int   group_next_progenitor_index;
        int   group_dominant_progenitor_file =-1;
        int   group_dominant_progenitor_index=-1;

        // Parse the input.  If we are writing extended files,
        //    then the input is not extended, and vice versa
        if(flag_write_extended || flag_write_pointers){
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

           // If the file offset is negative, so should the index be
           if(group_file_offset<0)
              group_file_index=-1;

           // Compute other stats
           group_n_particles       =groups[i_write%n_wrap][i_group].n_particles;
           group_n_particles_parent=groups[i_write%n_wrap][i_group].n_particles_parent;
           group_n_particles_desc  =set_match_n_particles(&(groups[i_write%n_wrap][i_group].descendant));
           if((groups[i_write%n_wrap][i_group].descendant.halo)!=NULL)
              group_n_particles_proj=set_match_n_particles(&(groups[i_write%n_wrap][i_group].descendant.halo->first_progenitor));
           else
              group_n_particles_proj=-1;
           group_score_desc     =set_match_score(   &(groups[i_write%n_wrap][i_group].descendant));
           group_score_prog     =set_match_score(   &(groups[i_write%n_wrap][i_group].first_progenitor));
           group_snap_backmatch =set_match_snapshot(&(groups[i_write%n_wrap][i_group].bridge_backmatch));
           group_file_backmatch =set_match_file(    &(groups[i_write%n_wrap][i_group].bridge_backmatch));
           group_index_backmatch=set_match_index(   &(groups[i_write%n_wrap][i_group].bridge_backmatch));
           group_id_backmatch   =set_match_id(      &(groups[i_write%n_wrap][i_group].bridge_backmatch));
           group_snap_forematch =set_match_snapshot(&(groups[i_write%n_wrap][i_group].forematch_first));
           group_file_forematch =set_match_file(    &(groups[i_write%n_wrap][i_group].forematch_first));
           group_index_forematch=set_match_index(   &(groups[i_write%n_wrap][i_group].forematch_first));
           group_first_progenitor_file =set_match_file( &(groups[i_write%n_wrap][i_group].first_progenitor));
           group_first_progenitor_index=set_match_index(&(groups[i_write%n_wrap][i_group].first_progenitor));
           group_next_progenitor_file  =set_match_file( &(groups[i_write%n_wrap][i_group].next_progenitor));
           group_next_progenitor_index =set_match_index(&(groups[i_write%n_wrap][i_group].next_progenitor));
        }
        else{
           tree_horizontal_extended_info **groups;
           groups                         =(tree_horizontal_extended_info **)groups_in;
           group_id                       =groups[i_write%n_wrap][i_group].id;
           group_type                     =groups[i_write%n_wrap][i_group].type;
           group_descendant_id            =groups[i_write%n_wrap][i_group].descendant_id;
           group_tree_id                  =groups[i_write%n_wrap][i_group].tree_id;
           group_file_offset              =groups[i_write%n_wrap][i_group].descendant_file_offset;
           group_file_index               =groups[i_write%n_wrap][i_group].descendant_index;
           group_dominant_progenitor_file =groups[i_write%n_wrap][i_group].dominant_progenitor_file;
           group_dominant_progenitor_index=groups[i_write%n_wrap][i_group].dominant_progenitor_index;
        }

        // Write group information to the horizontal tree files
        SID_fwrite(&group_id,           sizeof(int),1,&fp_trees_out);
        SID_fwrite(&group_type,         sizeof(int),1,&fp_trees_out);
        SID_fwrite(&group_descendant_id,sizeof(int),1,&fp_trees_out);
        SID_fwrite(&group_tree_id,      sizeof(int),1,&fp_trees_out);
        SID_fwrite(&group_file_offset,  sizeof(int),1,&fp_trees_out);
        SID_fwrite(&group_file_index,   sizeof(int),1,&fp_trees_out);
        if(flag_write_extended){
           SID_fwrite(&group_n_particles,           sizeof(int),  1,&fp_trees_out);
           SID_fwrite(&group_n_particles_parent,    sizeof(int),  1,&fp_trees_out);
           SID_fwrite(&group_n_particles_desc,      sizeof(int),  1,&fp_trees_out);
           SID_fwrite(&group_n_particles_proj,      sizeof(int),  1,&fp_trees_out);
           SID_fwrite(&group_score_desc,            sizeof(float),1,&fp_trees_out);
           SID_fwrite(&group_score_prog,            sizeof(float),1,&fp_trees_out);
           SID_fwrite(&group_snap_backmatch,        sizeof(int),  1,&fp_trees_out);
           SID_fwrite(&group_file_backmatch,        sizeof(int),  1,&fp_trees_out);
           SID_fwrite(&group_index_backmatch,       sizeof(int),  1,&fp_trees_out);
           SID_fwrite(&group_id_backmatch,          sizeof(int),  1,&fp_trees_out);
           SID_fwrite(&group_first_progenitor_file, sizeof(int),  1,&fp_trees_out);  
           SID_fwrite(&group_first_progenitor_index,sizeof(int),  1,&fp_trees_out);
           SID_fwrite(&group_next_progenitor_file,  sizeof(int),  1,&fp_trees_out);
           SID_fwrite(&group_next_progenitor_index, sizeof(int),  1,&fp_trees_out);
        }
        SID_fwrite(&(n_subgroups_group[i_write%n_wrap][i_group]),sizeof(int),1,&fp_trees_out);

        // Write the forward and back-match pointers
        //   needed, for example, to investigate
        //   emerged and fragmented halos
        if(flag_write_pointers){
           SID_fwrite(&group_tree_id,        sizeof(int),1,&fp_backmatch_ptrs_out);
           SID_fwrite(&group_file_backmatch, sizeof(int),1,&fp_backmatch_ptrs_out);
           SID_fwrite(&group_index_backmatch,sizeof(int),1,&fp_backmatch_ptrs_out);
           SID_fwrite(&(n_subgroups_group[i_write%n_wrap][i_group]),sizeof(int),1,&fp_backmatch_ptrs_out);
           SID_fwrite(&group_tree_id,        sizeof(int),1,&fp_bridge_ptrs_out);
           SID_fwrite(&group_file_forematch, sizeof(int),1,&fp_bridge_ptrs_out);
           SID_fwrite(&group_index_forematch,sizeof(int),1,&fp_bridge_ptrs_out);
           SID_fwrite(&(n_subgroups_group[i_write%n_wrap][i_group]),sizeof(int),1,&fp_bridge_ptrs_out);
        }

        // Write dominant progenitor files
        if(flag_write_dom_prog){
           SID_fwrite(&group_dominant_progenitor_file, sizeof(int),1,&fp_dom_prog_out);
           SID_fwrite(&group_dominant_progenitor_index,sizeof(int),1,&fp_dom_prog_out);
        }

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
           int   subgroup_snap_backmatch;
           int   subgroup_file_backmatch;
           int   subgroup_index_backmatch;
           int   subgroup_id_backmatch;
           int   subgroup_snap_forematch;
           int   subgroup_file_forematch;
           int   subgroup_index_forematch;
           float subgroup_score_desc;
           float subgroup_score_prog;
           int   subgroup_first_progenitor_file;
           int   subgroup_first_progenitor_index;
           int   subgroup_next_progenitor_file;
           int   subgroup_next_progenitor_index;
           int   subgroup_dominant_progenitor_file =-1;
           int   subgroup_dominant_progenitor_index=-1;

           // Parse the input.  If we are writing extended files,
           //    then the input is not extended, and vice versa
           if(flag_write_extended || flag_write_pointers){
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

              // If the file offset is negative, so should the index be
              if(subgroup_file_offset<0)
                 subgroup_file_index=-1;

              // Compute other stats
              subgroup_n_particles       =subgroups[i_write%n_wrap][i_subgroup].n_particles;
              subgroup_n_particles_parent=subgroups[i_write%n_wrap][i_subgroup].n_particles_parent;
              subgroup_n_particles_desc  =set_match_n_particles(&(subgroups[i_write%n_wrap][i_subgroup].descendant));
              if((subgroups[i_write%n_wrap][i_subgroup].descendant.halo)!=NULL)
                 subgroup_n_particles_proj=set_match_n_particles(&(subgroups[i_write%n_wrap][i_subgroup].descendant.halo->first_progenitor));
              else
                 subgroup_n_particles_proj=-1;
              subgroup_score_desc     =set_match_score(   &(subgroups[i_write%n_wrap][i_subgroup].descendant));
              subgroup_score_prog     =set_match_score(   &(subgroups[i_write%n_wrap][i_subgroup].first_progenitor));
              subgroup_snap_backmatch =set_match_snapshot(&(subgroups[i_write%n_wrap][i_subgroup].bridge_backmatch));
              subgroup_file_backmatch =set_match_file(    &(subgroups[i_write%n_wrap][i_subgroup].bridge_backmatch));
              subgroup_index_backmatch=set_match_index(   &(subgroups[i_write%n_wrap][i_subgroup].bridge_backmatch));
              subgroup_id_backmatch   =set_match_id(      &(subgroups[i_write%n_wrap][i_subgroup].bridge_backmatch));
              subgroup_snap_forematch =set_match_snapshot(&(subgroups[i_write%n_wrap][i_subgroup].forematch_first));
              subgroup_file_forematch =set_match_file(    &(subgroups[i_write%n_wrap][i_subgroup].forematch_first));
              subgroup_index_forematch=set_match_index(   &(subgroups[i_write%n_wrap][i_subgroup].forematch_first));
              subgroup_first_progenitor_file =set_match_file( &(subgroups[i_write%n_wrap][i_subgroup].first_progenitor));
              subgroup_first_progenitor_index=set_match_index(&(subgroups[i_write%n_wrap][i_subgroup].first_progenitor));
              subgroup_next_progenitor_file  =set_match_file( &(subgroups[i_write%n_wrap][i_subgroup].next_progenitor));
              subgroup_next_progenitor_index =set_match_index(&(subgroups[i_write%n_wrap][i_subgroup].next_progenitor));
           }
           else{
              tree_horizontal_extended_info **subgroups;
              subgroups                         =(tree_horizontal_extended_info **)subgroups_in;
              subgroup_id                       =subgroups[i_write%n_wrap][i_subgroup].id;
              subgroup_type                     =subgroups[i_write%n_wrap][i_subgroup].type;
              subgroup_descendant_id            =subgroups[i_write%n_wrap][i_subgroup].descendant_id;
              subgroup_tree_id                  =subgroups[i_write%n_wrap][i_subgroup].tree_id;
              subgroup_file_offset              =subgroups[i_write%n_wrap][i_subgroup].descendant_file_offset;
              subgroup_file_index               =subgroups[i_write%n_wrap][i_subgroup].descendant_index;
              subgroup_dominant_progenitor_file =subgroups[i_write%n_wrap][i_subgroup].dominant_progenitor_file;
              subgroup_dominant_progenitor_index=subgroups[i_write%n_wrap][i_subgroup].dominant_progenitor_index;
           }

           // Write subgroup information to the horizontal trees files
           SID_fwrite(&subgroup_id,           sizeof(int),1,&fp_trees_out);
           SID_fwrite(&subgroup_type,         sizeof(int),1,&fp_trees_out);
           SID_fwrite(&subgroup_descendant_id,sizeof(int),1,&fp_trees_out);
           SID_fwrite(&subgroup_tree_id,      sizeof(int),1,&fp_trees_out);
           SID_fwrite(&subgroup_file_offset,  sizeof(int),1,&fp_trees_out);
           SID_fwrite(&subgroup_file_index,   sizeof(int),1,&fp_trees_out);
           if(flag_write_extended){
              SID_fwrite(&subgroup_n_particles,           sizeof(int),  1,&fp_trees_out);
              SID_fwrite(&subgroup_n_particles_parent,    sizeof(int),  1,&fp_trees_out);
              SID_fwrite(&subgroup_n_particles_desc,      sizeof(int),  1,&fp_trees_out);
              SID_fwrite(&subgroup_n_particles_proj,      sizeof(int),  1,&fp_trees_out);
              SID_fwrite(&subgroup_score_desc,            sizeof(float),1,&fp_trees_out);
              SID_fwrite(&subgroup_score_prog,            sizeof(float),1,&fp_trees_out);
              SID_fwrite(&subgroup_snap_backmatch,        sizeof(int),  1,&fp_trees_out);
              SID_fwrite(&subgroup_file_backmatch,        sizeof(int),  1,&fp_trees_out);
              SID_fwrite(&subgroup_index_backmatch,       sizeof(int),  1,&fp_trees_out);
              SID_fwrite(&subgroup_id_backmatch,          sizeof(int),  1,&fp_trees_out);
              SID_fwrite(&subgroup_first_progenitor_file, sizeof(int),  1,&fp_trees_out);  
              SID_fwrite(&subgroup_first_progenitor_index,sizeof(int),  1,&fp_trees_out);
              SID_fwrite(&subgroup_next_progenitor_file,  sizeof(int),  1,&fp_trees_out);
              SID_fwrite(&subgroup_next_progenitor_index, sizeof(int),  1,&fp_trees_out);
           }

           // Write the forward and back-match pointers
           //   needed, for example, to investigate
           //   emerged and fragmented halos
           if(flag_write_pointers){
              SID_fwrite(&subgroup_tree_id,        sizeof(int),1,&fp_backmatch_ptrs_out);
              SID_fwrite(&subgroup_file_backmatch, sizeof(int),1,&fp_backmatch_ptrs_out);
              SID_fwrite(&subgroup_index_backmatch,sizeof(int),1,&fp_backmatch_ptrs_out);
              SID_fwrite(&subgroup_tree_id,        sizeof(int),1,&fp_bridge_ptrs_out);
              SID_fwrite(&subgroup_file_forematch, sizeof(int),1,&fp_bridge_ptrs_out);
              SID_fwrite(&subgroup_index_forematch,sizeof(int),1,&fp_bridge_ptrs_out);
           }

           //Write the dominanant progenitor pointers
           if(flag_write_dom_prog){
              SID_fwrite(&subgroup_dominant_progenitor_file, sizeof(int),1,&fp_dom_prog_out);
              SID_fwrite(&subgroup_dominant_progenitor_index,sizeof(int),1,&fp_dom_prog_out);
           }
        }
     }
   }
   SID_fclose(&fp_trees_out);
   if(flag_write_pointers){
      SID_fclose(&fp_backmatch_ptrs_out);
      SID_fclose(&fp_bridge_ptrs_out);
   }
   if(flag_write_dom_prog){
      SID_fclose(&fp_dom_prog_out);
   }
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
            write_trees_horizontal_log_file(filename_log,i_write,j_write,i_k_match,n_k_match,&stats,a_list,cosmo,TRUE);
         else
            write_trees_horizontal_log_file(filename_log,i_write,j_write,i_k_match,n_k_match,&stats,a_list,cosmo,FALSE);
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
   SID_log("Done.",SID_LOG_CLOSE);
}

