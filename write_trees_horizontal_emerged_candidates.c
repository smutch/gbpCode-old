#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees.h>

void write_trees_horizontal_emerged_candidates(int                   i_read,
                                               int                   n_halos_i,
                                               tree_horizontal_info *halos_i,
                                               char                 *group_text_prefix,
                                               char                 *filename_output_dir,
                                               int                   flag_start_new_file){
   char  filename_output_dir_horizontal[MAX_FILENAME_LENGTH];
   char  filename_output_dir_horizontal_cases[MAX_FILENAME_LENGTH];
   char  filename_output_file_root[MAX_FILENAME_LENGTH];
   char  filename_matching_out[MAX_FILENAME_LENGTH];
   FILE *fp_matching_out;

   SID_log("Writing candidate emerged halo information...",SID_LOG_OPEN|SID_LOG_TIMER);
   sprintf(filename_output_dir_horizontal,      "%s/horizontal",   filename_output_dir);
   sprintf(filename_output_dir_horizontal_cases,"%s/special_cases",filename_output_dir_horizontal);
   mkdir(filename_output_dir,                   02755);
   mkdir(filename_output_dir_horizontal,        02755);
   mkdir(filename_output_dir_horizontal_cases,  02755);
   strcpy(filename_output_file_root,filename_output_dir);
   strip_path(filename_output_file_root);
   sprintf(filename_matching_out,"%s/%sgroup_emerged_candidates.txt",filename_output_dir_horizontal_cases,group_text_prefix);
   if(flag_start_new_file){
      fp_matching_out=fopen(filename_matching_out,"w");
      int i_column=1;
      fprintf(fp_matching_out,"# (%02d): Bridged halo snapshot\n",                           i_column++);
      fprintf(fp_matching_out,"# (%02d): Bridged halo index\n",                              i_column++);
      fprintf(fp_matching_out,"# (%02d): Bridged halo ID\n",                                 i_column++);
      fprintf(fp_matching_out,"# (%02d): Bridged halo descendant snapshot\n",                i_column++);
      fprintf(fp_matching_out,"# (%02d): Bridged halo descendant index\n",                   i_column++);
      fprintf(fp_matching_out,"# (%02d): Bridged halo descendant ID\n",                      i_column++);
      fprintf(fp_matching_out,"# (%02d): Emerged candidate number\n",                        i_column++);
      fprintf(fp_matching_out,"# (%02d): Emerged halo candidate snapshot\n",                 i_column++);
      fprintf(fp_matching_out,"# (%02d): Emerged halo candidate index\n",                    i_column++);
      fprintf(fp_matching_out,"# (%02d): Emerged halo candidate ID\n",                       i_column++);
      fprintf(fp_matching_out,"# (%02d): Bridge->descendant match score\n",                  i_column++);
      fprintf(fp_matching_out,"# (%02d): Emerged candidate's back-match score\n",            i_column++);
      fprintf(fp_matching_out,"# (%02d): Number of particles in bridged halo\n",             i_column++);
      fprintf(fp_matching_out,"# (%02d): Number of particles in bridged halo's descendant\n",i_column++);
      fprintf(fp_matching_out,"# (%02d): Number of particles in emerged candidate\n",        i_column++);
   }
   else
      fp_matching_out=fopen(filename_matching_out,"a");
   int i_halo;
   for(i_halo=0;i_halo<n_halos_i;i_halo++){
      // Print bridge info
      if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGED)){
         int j_halo;
         for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;j_halo++){
            bridge_info *bridge;
            bridge=&(halos_i[i_halo].bridges[j_halo]);
            fprintf(fp_matching_out,"%4d %7d %7d   %4d %7d %7d   %3d   %4d %7d %7d   %10.4f %10.4f   %7d %7d %7d\n",
                    i_read,
                    i_halo,
                    halos_i[i_halo].id,
                    set_match_snapshot(&(halos_i[i_halo].descendant)),
                    set_match_index(   &(halos_i[i_halo].descendant)),
                    set_match_id(      &(halos_i[i_halo].descendant)),
                    j_halo,
                    set_match_snapshot(bridge),
                    set_match_index(   bridge),
                    set_match_id(      bridge),
                    set_match_score(&(halos_i[i_halo].descendant)),
                    set_match_score(bridge),
                    halos_i[i_halo].n_particles,
                    set_match_n_particles(&(halos_i[i_halo].descendant)),
                    set_match_n_particles(bridge));
         }
      }
   }
   fclose(fp_matching_out);
   SID_log("Done.",SID_LOG_CLOSE);
}

