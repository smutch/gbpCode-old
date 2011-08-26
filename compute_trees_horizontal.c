#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void switch_horizontal_IDs_recursive(tree_horizontal_info **halos,int n_search,int j_file,int j_halo,int id_1,int id_2);
void switch_horizontal_IDs_recursive(tree_horizontal_info **halos,int n_search,int j_file,int j_halo,int id_1,int id_2){
  int current_file;
  int current_index;
  int current_file_temp;
  int current_index_temp;

  // Swap IDs here
  if(halos[j_file%n_search][j_halo].id==id_1)
    halos[j_file%n_search][j_halo].id=id_2;
  else if(halos[j_file%n_search][j_halo].id==id_2)
    halos[j_file%n_search][j_halo].id=id_1;

  // Walk the tree
  current_file =halos[j_file%n_search][j_halo].file_first_progenitor;
  current_index=halos[j_file%n_search][j_halo].index_first_progenitor;
  while(current_file>=0 && current_index>=0){
    switch_horizontal_IDs_recursive(halos,n_search,current_file,current_index,id_1,id_2);
    current_file_temp =current_file;
    current_index_temp=current_index;
    current_file      =halos[current_file_temp%n_search][current_index_temp].file_next_progenitor;
    current_index     =halos[current_file_temp%n_search][current_index_temp].index_next_progenitor;
  }
}

void compute_trees_horizontal_stats(tree_horizontal_info *halos,int n_halos,int n_halos_max,tree_horizontal_stats_info *stats);
void compute_trees_horizontal_stats(tree_horizontal_info *halos,int n_halos,int n_halos_max,tree_horizontal_stats_info *stats){
   int         n_invalid    =0;
   int         n_unprocessed=0;
   int         emerged_diff =0;
   int         i_halo;
   int         j_halo;
   
   stats->n_halos                     =n_halos;
   stats->n_simple                    =0;
   stats->n_mergers                   =0;
   stats->n_strayed                   =0;
   stats->n_sputtered                 =0;
   stats->n_dropped                   =0;
   stats->n_bridged                   =0;
   stats->n_bridge_progenitors        =0;
   stats->n_emerged                   =0;
   stats->n_emerged_lost              =0;
   stats->n_emerged_progenitors       =0;
   stats->max_strayed_size            =0;
   stats->max_sputtered_size          =0;
   stats->max_dropped_size            =0;
   stats->max_bridged_size            =0;
   stats->max_bridge_progenitor_size  =0;
   stats->max_emerged_size            =0;
   stats->max_emerged_lost_size       =0;
   stats->max_emerged_found_diff      =0;
   stats->max_emerged_found_diff_size =0;
   stats->max_emerged_progenitor_size =0;
   stats->max_id                      =0;
   
   for(i_halo=0;i_halo<n_halos_max;i_halo++){

      // Find maximum id
      stats->max_id=MAX(stats->max_id,halos[i_halo].id);

      // Compute statistcs for strays
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_SIMPLE))
         stats->n_simple++;

      // Compute statistcs for mergers
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_MERGER))
         stats->n_mergers++;

      // Compute statistcs for strays
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_STRAYED)){           
         stats->max_strayed_size=MAX(stats->max_strayed_size,halos[i_halo].n_particles);
         stats->n_strayed++;
      }

      // Compute statistcs for sputtered halos
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_SPUTTERED)){           
         stats->max_sputtered_size=MAX(stats->max_sputtered_size,halos[i_halo].n_particles);
         stats->n_sputtered++;
      }

      // Compute statistcs for dropped halos
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_DROPPED)){
         stats->max_dropped_size=MAX(stats->max_dropped_size,halos[i_halo].n_particles);
         stats->n_dropped++;
      }

      // Compute statistcs for bridged and emerged halos
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_BRIDGED)){
         stats->max_bridged_size=MAX(stats->max_bridged_size,halos[i_halo].n_particles);
         stats->n_bridged++;
      }
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_BRIDGE_PROGENITOR)){
         stats->max_bridge_progenitor_size=MAX(stats->max_bridge_progenitor_size,halos[i_halo].n_particles);
         stats->n_bridge_progenitors++;            
         if(!check_mode_for_flag(halos[i_halo].type,TREE_CASE_BRIDGE_DEFAULT)){
            stats->max_emerged_progenitor_size=MAX(stats->max_emerged_progenitor_size,halos[i_halo].n_particles);
            stats->n_emerged_progenitors++;
            if(halos[i_halo].n_particles>halos[i_halo].descendant.n_particles)
               emerged_diff=halos[i_halo].n_particles-halos[i_halo].descendant.n_particles;
            else
               emerged_diff=halos[i_halo].descendant.n_particles-halos[i_halo].n_particles;
            if(emerged_diff>stats->max_emerged_found_diff){
               stats->max_emerged_found_diff     =emerged_diff;
               stats->max_emerged_found_diff_size=halos[i_halo].n_particles;
            }
         }
      }
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_EMERGED)){
         stats->max_emerged_size=MAX(stats->max_emerged_size,halos[i_halo].n_particles);
         stats->n_emerged++;
         if(!check_mode_for_flag(halos[i_halo].type,TREE_CASE_FOUND)){
            stats->max_emerged_lost_size=MAX(stats->max_emerged_lost_size,halos[i_halo].n_particles);
            stats->n_emerged_lost++;
         }
      }

      // Count n_unprocessed and n_invalid so we can perform a couple sanity checks
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_INVALID)){
         n_invalid++;
         if((halos[i_halo].type-TREE_CASE_INVALID)!=0)
            SID_trap_error("An out-of-bounds halo has been manipulated",ERROR_LOGIC);
      }
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_UNPROCESSED))
         n_unprocessed++;
   }
   if(n_halos!=(n_halos_max-n_invalid))
      SID_trap_error("There is an incorrect number of out-of-bounds halos (i.e. %d!=%d)",ERROR_LOGIC,n_halos,(n_halos_max-n_invalid));
   if(n_unprocessed!=0)
      SID_trap_error("A number of halos (%d) have not been marked as processed",ERROR_LOGIC,n_unprocessed);
}

void write_trees_horizontal(tree_horizontal_info **groups,   int n_groups,    int n_groups_max,   int *n_particles_groups,
                            tree_horizontal_info **subgroups,int n_subgroups, int n_subgroups_max,int *n_particles_subgroups,
                            int   **n_subgroups_group,
                            int     n_halos_max,
                            int     max_tree_id_subgroup,
                            int     max_tree_id_group,
                            int     i_write,
                            int     j_write,
                            int     l_write,
                            int     n_search,
                            char   *filename_cat_root_in,
                            char   *filename_output_dir,
                            double *a_list,
                            cosmo_info **cosmo,
                            int     n_k_match);
void write_trees_horizontal(tree_horizontal_info **groups,   int n_groups,    int n_groups_max,   int *n_particles_groups,
                            tree_horizontal_info **subgroups,int n_subgroups, int n_subgroups_max,int *n_particles_subgroups,
                            int   **n_subgroups_group,
                            int     n_halos_max,
                            int     max_tree_id_subgroup,
                            int     max_tree_id_group,
                            int     i_write,
                            int     j_write,
                            int     l_write,
                            int     n_search,
                            char   *filename_cat_root_in,
                            char   *filename_output_dir,
                            double *a_list,
                            cosmo_info **cosmo,
                            int     n_k_match){
   char        filename_output_dir_horizontal[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_trees[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_groups[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_groups_stats[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_groups_properties[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_subgroups[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_subgroups_stats[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_subgroups_properties[MAX_FILENAME_LENGTH];
   char        filename_output_file_root[MAX_FILENAME_LENGTH];
   char        filename_matches_out[MAX_FILENAME_LENGTH];
   char        filename_groups[MAX_FILENAME_LENGTH];
   char        filename_subgroups[MAX_FILENAME_LENGTH];
   char        filename_subgroup_properties_in[MAX_FILENAME_LENGTH];
   char        filename_group_properties_in[MAX_FILENAME_LENGTH];
   char        filename_subgroup_properties_out[MAX_FILENAME_LENGTH];
   char        filename_group_properties_out[MAX_FILENAME_LENGTH];
   char        filename_log[MAX_FILENAME_LENGTH];
   char        filename_mergers_out[MAX_FILENAME_LENGTH];
   char        filename_strayed_out[MAX_FILENAME_LENGTH];
   char        filename_sputtered_out[MAX_FILENAME_LENGTH];
   char        filename_dropped_out[MAX_FILENAME_LENGTH];
   char        filename_bridged_out[MAX_FILENAME_LENGTH];
   char        filename_bridged_progenitor_out[MAX_FILENAME_LENGTH];
   char        filename_emerged_progenitor_out[MAX_FILENAME_LENGTH];
   char        filename_emerged_found_out[MAX_FILENAME_LENGTH];
   char        filename_emerged_unfound_out[MAX_FILENAME_LENGTH];
   char       *filename_output_dir_stats;
   FILE       *fp_mergers_out;
   FILE       *fp_strayed_out;
   FILE       *fp_sputtered_out;
   FILE       *fp_dropped_out;
   FILE       *fp_bridged_out;
   FILE       *fp_bridged_progenitor_out;
   FILE       *fp_emerged_progenitor_out;
   FILE       *fp_emerged_found_out;
   FILE       *fp_emerged_unfound_out;
   int         i_halo;
   char        group_text_prefix[5];
   SID_fp      fp_matches_out;
   SID_fp      fp_group_properties_out;
   SID_fp      fp_subgroup_properties_out;
   FILE       *fp_group_properties_in;
   FILE       *fp_subgroup_properties_in;
   FILE       *fp;
   int         k_match;
   tree_horizontal_info  *halos;
   tree_horizontal_info **halos_all;
   int         n_halos;
   int         j_halo;
   int         k_halo;
   int        *n_particles_parent;
   int         file_offset;
   int         i_subgroup;
   int         j_subgroup;
   int         i_group;
   int         j_group;
   int         i_k_match;
   int         j_k_match;
   int         i_column;
   int         n_particles_residual;
   int         n_p_largest_main;
   int         n_particles_residual_main;
   double      dt_descendant;
   char       *line=NULL;
   int         line_length=0;
   halo_info                  *properties;
   tree_horizontal_stats_info  stats;

   SID_log("Writing results for snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,j_write);
   sprintf(filename_output_dir_horizontal,                     "%s/horizontal",filename_output_dir);
   sprintf(filename_output_dir_horizontal_trees,               "%s/trees",     filename_output_dir_horizontal);
   sprintf(filename_output_dir_horizontal_groups,              "%s/groups",    filename_output_dir_horizontal);
   sprintf(filename_output_dir_horizontal_groups_stats,        "%s/stats",     filename_output_dir_horizontal_groups);
   sprintf(filename_output_dir_horizontal_groups_properties,   "%s/properties",filename_output_dir_horizontal_groups);
   sprintf(filename_output_dir_horizontal_subgroups,           "%s/subgroups", filename_output_dir_horizontal);
   sprintf(filename_output_dir_horizontal_subgroups_stats,     "%s/stats",     filename_output_dir_horizontal_subgroups);
   sprintf(filename_output_dir_horizontal_subgroups_properties,"%s/properties",filename_output_dir_horizontal_subgroups);
   mkdir(filename_output_dir,                                02755);
   mkdir(filename_output_dir_horizontal,                     02755);
   mkdir(filename_output_dir_horizontal_trees,               02755);
   mkdir(filename_output_dir_horizontal_groups,              02755);
   mkdir(filename_output_dir_horizontal_groups_stats,        02755);
   mkdir(filename_output_dir_horizontal_groups_properties,   02755);
   mkdir(filename_output_dir_horizontal_subgroups,           02755);
   mkdir(filename_output_dir_horizontal_subgroups_stats,     02755);
   mkdir(filename_output_dir_horizontal_subgroups_properties,02755);
   strcpy(filename_output_file_root,filename_output_dir);
   strip_path(filename_output_file_root);
   sprintf(filename_matches_out,            "%s/%s.trees_horizontal_%d",                    filename_output_dir_horizontal_trees,            filename_output_file_root,j_write);
   sprintf(filename_group_properties_out,   "%s/%s.trees_horizontal_groups_properties_%d",  filename_output_dir_horizontal_groups_properties,filename_output_file_root,j_write);
   sprintf(filename_subgroup_properties_out,"%s/%s.trees_horizontal_subgroups_properties_%d",filename_output_dir_horizontal_subgroups_properties,filename_output_file_root,j_write);
   sprintf(filename_group_properties_in,    "%s_%03d.catalog_groups_properties",             filename_cat_root_in,j_write);
   sprintf(filename_subgroup_properties_in, "%s_%03d.catalog_subgroups_properties",          filename_cat_root_in,j_write);
   /*
   SID_log("filename_matches_out            ={%s}",SID_LOG_COMMENT,filename_matches_out);
   SID_log("filename_group_properties_out   ={%s}",SID_LOG_COMMENT,filename_group_properties_out);
   SID_log("filename_subgroup_properties_out={%s}",SID_LOG_COMMENT,filename_subgroup_properties_out);
   SID_log("filename_group_properties_in    ={%s}",SID_LOG_COMMENT,filename_group_properties_in);
   SID_log("filename_subgroup_properties_in ={%s}",SID_LOG_COMMENT,filename_subgroup_properties_in);
   */
   SID_fopen(filename_matches_out,            "w",&fp_matches_out);
   SID_fopen(filename_group_properties_out,   "w",&fp_group_properties_out);
   SID_fopen(filename_subgroup_properties_out,"w",&fp_subgroup_properties_out);
   SID_fwrite(&(n_groups),            sizeof(int),1,&fp_matches_out);
   SID_fwrite(&(n_subgroups),         sizeof(int),1,&fp_matches_out);
   SID_fwrite(&(n_halos_max),         sizeof(int),1,&fp_matches_out);
   SID_fwrite(&(max_tree_id_subgroup),sizeof(int),1,&fp_matches_out);
   SID_fwrite(&(max_tree_id_group),   sizeof(int),1,&fp_matches_out);
   properties=(halo_info *)SID_calloc(sizeof(halo_info)); // valgrind throws an error if we don't init the unused bits with calloc
   if(n_groups>0){
     fp_group_properties_in   =fopen(filename_group_properties_in,   "r");
     fp_subgroup_properties_in=fopen(filename_subgroup_properties_in,"r");
     fseeko(fp_group_properties_in,   4*sizeof(int),SEEK_CUR);
     fseeko(fp_subgroup_properties_in,4*sizeof(int),SEEK_CUR);
     for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
       // compute_trees_verticle wants the file offsets to be -ve for the roots, +ve everywhere else
       file_offset=groups[i_write%n_search][i_group].descendant.file-i_write;
       SID_fwrite(&(groups[i_write%n_search][i_group].id),           sizeof(int),1,&fp_matches_out);
       SID_fwrite(&(groups[i_write%n_search][i_group].descendant.id),sizeof(int),1,&fp_matches_out);
       SID_fwrite(&(groups[i_write%n_search][i_group].tree_id),      sizeof(int),1,&fp_matches_out);
       SID_fwrite(&(file_offset),                                    sizeof(int),1,&fp_matches_out);
       SID_fwrite(&(n_subgroups_group[i_write%n_search][i_group]),   sizeof(int),1,&fp_matches_out);
       read_group_properties(fp_group_properties_in,properties,i_group,j_write);
       if(groups[i_write%n_search][i_group].id>=0)
         SID_fwrite(properties,sizeof(halo_info),1,&fp_group_properties_out);
       for(j_subgroup=0;j_subgroup<n_subgroups_group[i_write%n_search][i_group];j_subgroup++,i_subgroup++){
         // compute_trees_verticle wants the file offsets to be -ve for the roots, +ve everywhere else
         file_offset=subgroups[i_write%n_search][i_subgroup].descendant.file-i_write;
         SID_fwrite(&(subgroups[i_write%n_search][i_subgroup].id),           sizeof(int),1,&fp_matches_out);
         SID_fwrite(&(subgroups[i_write%n_search][i_subgroup].descendant.id),sizeof(int),1,&fp_matches_out);
         SID_fwrite(&(subgroups[i_write%n_search][i_subgroup].tree_id),      sizeof(int),1,&fp_matches_out);
         SID_fwrite(&(file_offset),                                          sizeof(int),1,&fp_matches_out);
         read_group_properties(fp_subgroup_properties_in,properties,i_subgroup,j_write);
         if(subgroups[i_write%n_search][i_subgroup].id>=0)
           SID_fwrite(properties,sizeof(halo_info),1,&fp_subgroup_properties_out);
       }
     }
     fclose(fp_group_properties_in);
     fclose(fp_subgroup_properties_in);
   }
   SID_free(SID_FARG properties);
   SID_fclose(&fp_matches_out);
   SID_fclose(&fp_group_properties_out);
   SID_fclose(&fp_subgroup_properties_out);   
   
   // Write statistics to ascii files
   sprintf(filename_log,"%s/%s.log",filename_output_dir_horizontal,filename_output_file_root);
   for(i_k_match=0;i_k_match<n_k_match;i_k_match++){
      switch(i_k_match){
         case 0:
         compute_trees_horizontal_stats(subgroups[i_write%n_search],n_subgroups,n_subgroups_max,&stats);
         break;
         case 1:
         compute_trees_horizontal_stats(groups[i_write%n_search],   n_groups,   n_groups_max,   &stats);
         break;
      }
      if(i_k_match==0 && l_write==0){
         fp=fopen(filename_log,"w");
         i_column=1;
         fprintf(fp,"# (%02d): Expansion factor (a)\n",i_column++);
         fprintf(fp,"# (%02d): Snapshot filenumber\n", i_column++);
         for(j_k_match=0;j_k_match<n_k_match;j_k_match++){
            switch(j_k_match){
               case 0:
               sprintf(group_text_prefix,"sub");
               break;
               case 1:
               sprintf(group_text_prefix,"");
               break;
            }            
            fprintf(fp,"# (%02d): Maximum %sgroup ID\n",                    i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of %sgroups\n",                         i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of simple       %sgroups\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of merging      %sgroups\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of strayed      %sgroups\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of sputtering   %sgroups\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of dropped      %sgroups\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of bridged      %sgroups\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of bridged      %sgroups progenitors\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of emerged      %sgroups\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of emerged lost %sgroups\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of emerged      %sgroups progenitors\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest strayed   %sgroup\n",             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest sputtered %sgroup\n",             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest dropped   %sgroup\n",             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest bridged   %sgroup\n",             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest bridged   %sgroup progenitor\n",  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest emerged   %sgroup\n",             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest emerged   %sgroup lost\n",        i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest emerged   %sgroup found diff\n",       i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest emerged   %sgroup found diff size\n",  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest emerged   %sgroup progenitor size\n",  i_column++,group_text_prefix);
         }
         fclose(fp);
      }
      fp=fopen(filename_log,"a");
      if(i_k_match==0)
         fprintf(fp,"%le %4d",a_list[l_write],j_write);
      fprintf(fp," %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d",
         stats.max_id,
         stats.n_halos,         
         stats.n_simple,         
         stats.n_mergers,         
         stats.n_strayed,
         stats.n_sputtered,
         stats.n_dropped,
         stats.n_bridged,
         stats.n_bridge_progenitors,
         stats.n_emerged,
         stats.n_emerged_lost,
         stats.n_emerged_progenitors,
         stats.max_strayed_size,
         stats.max_sputtered_size,
         stats.max_dropped_size,
         stats.max_bridged_size,
         stats.max_bridge_progenitor_size,
         stats.max_emerged_size,
         stats.max_emerged_lost_size,
         stats.max_emerged_found_diff,
         stats.max_emerged_found_diff_size,
         stats.max_emerged_progenitor_size);
      if(i_k_match==n_k_match-1)
         fprintf(fp,"\n");
      fclose(fp);

      // Write special case information
      n_particles_parent=(int *)SID_malloc(sizeof(int)*MAX(n_groups,n_subgroups));
      for(k_match=0;k_match<n_k_match;k_match++){

         // Initialize a bunch of stuff depending on whether
         //   we are processing groups or subgroups
         switch(k_match){
            case 0:
               filename_output_dir_stats=filename_output_dir_horizontal_subgroups_stats;
               sprintf(group_text_prefix,"sub");
               halos               =subgroups[i_write%n_search];
               halos_all           =subgroups;
               n_halos             =n_subgroups;
               for(i_halo=0,k_halo=0;i_halo<n_groups;i_halo++){
                  for(j_halo=0;j_halo<n_subgroups_group[i_write%n_search][i_halo];j_halo++,k_halo++)
                     n_particles_parent[k_halo]=n_particles_groups[i_halo];
               }
               break;
            case 1:
               filename_output_dir_stats=filename_output_dir_horizontal_groups_stats;
               sprintf(group_text_prefix,"");
               halos               =groups[i_write%n_search];
               halos_all           =groups;
               n_halos             =n_groups;
               n_particles_parent  =(int *)SID_malloc(sizeof(int)*n_halos);
               for(i_halo=0;i_halo<n_halos;i_halo++)
                  n_particles_parent[i_halo]=n_particles_groups[i_halo];
               break;
         }
         sprintf(filename_mergers_out,           "%s/%s.%sgroups_mergers",            filename_output_dir_stats,filename_output_file_root,group_text_prefix);
         sprintf(filename_strayed_out,           "%s/%s.%sgroups_strays",             filename_output_dir_stats,filename_output_file_root,group_text_prefix);
         sprintf(filename_sputtered_out,         "%s/%s.%sgroups_sputters",           filename_output_dir_stats,filename_output_file_root,group_text_prefix);
         sprintf(filename_dropped_out,           "%s/%s.%sgroups_drops",              filename_output_dir_stats,filename_output_file_root,group_text_prefix);
         sprintf(filename_bridged_out,           "%s/%s.%sgroups_bridges",            filename_output_dir_stats,filename_output_file_root,group_text_prefix);      
         sprintf(filename_bridged_progenitor_out,"%s/%s.%sgroups_bridged_progenitors",filename_output_dir_stats,filename_output_file_root,group_text_prefix);      
         sprintf(filename_emerged_progenitor_out,"%s/%s.%sgroups_emerged_progenitors",filename_output_dir_stats,filename_output_file_root,group_text_prefix);      
         sprintf(filename_emerged_found_out,     "%s/%s.%sgroups_emerged_found",      filename_output_dir_stats,filename_output_file_root,group_text_prefix);      
         sprintf(filename_emerged_unfound_out,   "%s/%s.%sgroups_emerged_unfound",    filename_output_dir_stats,filename_output_file_root,group_text_prefix);      
         if(l_write==0){
            fp_mergers_out           =fopen(filename_mergers_out,           "w");
            fp_strayed_out           =fopen(filename_strayed_out,           "w");
            fp_sputtered_out         =fopen(filename_sputtered_out,         "w");
            fp_dropped_out           =fopen(filename_dropped_out,           "w");
            fp_bridged_out           =fopen(filename_bridged_out,           "w");
            fp_bridged_progenitor_out=fopen(filename_bridged_progenitor_out,"w");
            fp_emerged_progenitor_out=fopen(filename_emerged_progenitor_out,"w");
            fp_emerged_found_out     =fopen(filename_emerged_found_out,     "w");
            fp_emerged_unfound_out   =fopen(filename_emerged_unfound_out,   "w");
         }
         else{
            fp_mergers_out           =fopen(filename_mergers_out,           "a");
            fp_strayed_out           =fopen(filename_strayed_out,           "a");
            fp_sputtered_out         =fopen(filename_sputtered_out,         "a");
            fp_dropped_out           =fopen(filename_dropped_out,           "a");
            fp_bridged_out           =fopen(filename_bridged_out,           "a");
            fp_bridged_progenitor_out=fopen(filename_bridged_progenitor_out,"a");
            fp_emerged_progenitor_out=fopen(filename_emerged_progenitor_out,"a");
            fp_emerged_found_out     =fopen(filename_emerged_found_out,     "a");
            fp_emerged_unfound_out   =fopen(filename_emerged_unfound_out,   "a");
         }

         // Loop over each halo
         for(i_halo=0;i_halo<n_halos;i_halo++){

            // Compute the time between the halo and it's descendant
            if(halos[i_halo].descendant.file>i_write){ // This isn't true for first snapshot for instance
              if(l_write<(halos[i_halo].descendant.file-i_write))
                 SID_trap_error("Unreasonable descendant file offset: %d %d %d",ERROR_LOGIC,l_write,halos[i_halo].descendant.file,i_write);
              dt_descendant=deltat_a(cosmo,a_list[l_write],a_list[l_write-(halos[i_halo].descendant.file-i_write)])/S_PER_YEAR;
            }
            else
              dt_descendant=0.;

            // Write mergers
            fp=fp_mergers_out;
            if(l_write==0 && i_halo==0){
               i_column=1;
               fprintf(fp,"# (%02d): %sgroup expansion factor\n",                                   i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): time between %sgroup and its descendant [yrs]\n",               i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot number\n",                                    i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot index\n",                                     i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup id\n",                                                 i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup match type\n",                                         i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup\n",                         i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",            i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n",             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup merger's main progenitor\n",i_column++,group_text_prefix);
            }
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_MERGER)){
               fprintf(fp,"%10.3le %10.3le %5d %7d %7d %8d %3d %8d %8d %8d %8d\n",
                  a_list[l_write],
                  dt_descendant,
                  j_write,
                  i_halo,
                  halos[i_halo].id,
                  halos[i_halo].type,
                  halos[i_halo].descendant.file-i_write,
                  halos[i_halo].n_particles,
                  halos[i_halo].descendant.n_particles,
                  n_particles_parent[i_halo],
                  halos[i_halo].size_main_progenitor);
            }

            // Write strays
            fp=fp_strayed_out;
            if(l_write==0 && i_halo==0){
               i_column=1;
               fprintf(fp,"# (%02d): %sgroup expansion factor\n",                       i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot number\n",                        i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot index\n",                         i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup match type\n",                             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup\n",             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n", i_column++,group_text_prefix);
            }
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_STRAYED)){
               fprintf(fp,"%10.3le %5d %7d %8d %8d %8d\n",
                  a_list[l_write],
                  j_write,
                  i_halo,
                  halos[i_halo].type,
                  halos[i_halo].n_particles,
                  n_particles_parent[i_halo]);
            }

            // Write sputtered halos
            fp=fp_sputtered_out;
            if(l_write==0 && i_halo==0){
               i_column=1;
               fprintf(fp,"# (%02d): %sgroup expansion factor\n",                       i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): time between %sgroup and its descendant [yrs]\n",   i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot number\n",                        i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot index\n",                         i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup match type\n",                             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                 i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup\n",             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n", i_column++,group_text_prefix);
            }
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_SPUTTERED)){
               fprintf(fp,"%10.3le %10.3le %5d %7d %8d %3d %8d %8d %8d\n",
                  a_list[l_write],
                  dt_descendant,
                  j_write,
                  i_halo,
                  halos[i_halo].type,
                  halos[i_halo].descendant.file-i_write,
                  halos[i_halo].n_particles,
                  halos[i_halo].descendant.n_particles,
                  n_particles_parent[i_halo]);
            }

            // Write dropped halos
            fp=fp_dropped_out;
            if(l_write==0 && i_halo==0){
               i_column=1;
               fprintf(fp,"# (%02d): %sgroup expansion factor\n",                       i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): time between %sgroup and its descendant [yrs]\n",   i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot number\n",                        i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot index\n",                         i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup id\n",                                     i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup match type\n",                             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                 i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup\n",             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n", i_column++,group_text_prefix);
            }
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_DROPPED)){
               fprintf(fp,"%10.3le %10.3le %5d %7d %7d %8d %3d %8d %8d %8d\n",
                  a_list[l_write],
                  dt_descendant,
                  j_write,
                  i_halo,
                  halos[i_halo].id,
                  halos[i_halo].type,
                  halos[i_halo].descendant.file-i_write,
                  halos[i_halo].n_particles,
                  halos[i_halo].descendant.n_particles,
                  n_particles_parent[i_halo]);
            }

            // Write bridged and emerged halos
            fp=fp_bridged_out;
            if(l_write==0 && i_halo==0){
               i_column=1;
               fprintf(fp,"# (%02d): %sgroup expansion factor\n",                           i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): time between %sgroup and its descendant [yrs]\n",       i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot number\n",                            i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot index\n",                             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup id\n",                                         i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup match type\n",                                 i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                     i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup\n",                 i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",    i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n",     i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of %sgroups emerging from this bridge\n",       i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the largest emerged %sgroup\n", i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the largest emerged %sgroup's main progenitor\n",  i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the other   emerged %sgroups\n",                   i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the other   emerged %sgroups's main progenitors\n",i_column++,group_text_prefix);
            }
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_BRIDGED)){
               n_p_largest_main=halos_all[halos[i_halo].bridges[0].file%n_search][halos[i_halo].bridges[0].index].size_main_progenitor;
               for(j_halo=1,n_particles_residual=0,n_particles_residual_main=0;j_halo<halos[i_halo].n_bridges;j_halo++){
                 n_particles_residual     +=halos[i_halo].bridges[j_halo].n_particles;
                 n_particles_residual_main+=halos_all[halos[i_halo].bridges[j_halo].file%n_search][halos[i_halo].bridges[j_halo].index].size_main_progenitor;
               }
               fprintf(fp,"%10.3le %10.3le %5d %7d %7d %8d %3d %8d %8d %8d %8d %8d %8d %8d %8d\n",
                  a_list[l_write],
                  dt_descendant,
                  j_write,
                  i_halo,
                  halos[i_halo].id,
                  halos[i_halo].type,
                  halos[i_halo].descendant.file-i_write,
                  halos[i_halo].n_particles,
                  halos[i_halo].descendant.n_particles,
                  n_particles_parent[i_halo],
                  halos[i_halo].n_bridges,
                  halos[i_halo].bridges[0].n_particles,
                  n_p_largest_main,
                  n_particles_residual,
                  n_particles_residual_main);
            }
            fp=fp_emerged_progenitor_out;
            if(l_write==0 && i_halo==0){
               i_column=1;
               fprintf(fp,"# (%02d): %sgroup expansion factor\n",                       i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): time between %sgroup and its descendant [yrs]\n",   i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot number\n",                        i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot index\n",                         i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup id\n",                                     i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup match type\n",                             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                 i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup\n",             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n", i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): id of bridge matched to by %sgroup\n",             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): snapshot of bridge matched to by %sgroup\n",       i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): index    of bridge matched to by %sgroup\n",          i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's bridge\n",    i_column++,group_text_prefix);
            }
            fp=fp_bridged_progenitor_out;
            if(l_write==0 && i_halo==0){
               i_column=1;
               fprintf(fp,"# (%02d): %sgroup expansion factor\n",                       i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): time between %sgroup and its descendant [yrs]\n",   i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot number\n",                        i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot index\n",                         i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup id\n",                                     i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup match type\n",                             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                 i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup\n",             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n", i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): id of bridge matched to by %sgroup\n",             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): snapshot of bridge matched to by %sgroup\n",       i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): index    of bridge matched to by %sgroup\n",          i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's bridge\n",    i_column++,group_text_prefix);
            }
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_BRIDGE_PROGENITOR)){
               if(!check_mode_for_flag(halos[i_halo].type,TREE_CASE_BRIDGE_DEFAULT)){
                  fp=fp_emerged_progenitor_out;
                  fprintf(fp,"%10.3le %10.3le %5d %7d %7d %8d %3d %8d %8d %8d %7d %4d %8d %8d\n",
                     a_list[l_write],
                     dt_descendant,
                     j_write,
                     i_halo,
                     halos[i_halo].id,
                     halos[i_halo].type,
                     halos[i_halo].descendant.file-i_write,
                     halos[i_halo].n_particles,
                     halos[i_halo].descendant.n_particles,
                     n_particles_parent[i_halo],
                     halos[i_halo].bridge_match->id,
                     halos[i_halo].bridge_match->file,
                     halos[i_halo].bridge_match->index,
                     halos[i_halo].bridge_match->n_particles);
               }
               else{
                  fp=fp_bridged_progenitor_out;
                  fprintf(fp,"%10.3le %10.3le %5d %7d %7d %8d %3d %8d %8d %8d %7d %4d %8d %8d\n",
                     a_list[l_write],
                     dt_descendant,
                     j_write,
                     i_halo,
                     halos[i_halo].id,
                     halos[i_halo].type,
                     halos[i_halo].descendant.file-i_write,
                     halos[i_halo].n_particles,
                     halos[i_halo].descendant.n_particles,
                     n_particles_parent[i_halo],
                     halos[i_halo].bridge_match->id,
                     halos[i_halo].bridge_match->file,
                     halos[i_halo].bridge_match->index,
                     halos[i_halo].bridge_match->n_particles);
               }
            }
            fp=fp_emerged_unfound_out;
            if(l_write==0 && i_halo==0){
               i_column=1;
               fprintf(fp,"# (%02d): %sgroup expansion factor\n",                       i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): time between %sgroup and its descendant [yrs]\n",   i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot number\n",                        i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot index\n",                         i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup id\n",                                     i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup match type\n",                             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                 i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup\n",             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n", i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): id of bridge back-matched to by %sgroup\n",                 i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): snapshot of bridge back-matched to by %sgroup\n",           i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): index    of bridge back-matched to by %sgroup\n",              i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's back-matched bridge\n",i_column++,group_text_prefix);
            }
            fp=fp_emerged_found_out;
            if(l_write==0 && i_halo==0){
               i_column=1;
               fprintf(fp,"# (%02d): %sgroup expansion factor\n",                       i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): time between %sgroup and its descendant [yrs]\n",   i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot number\n",                        i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup snapshot index\n",                         i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup id\n",                                     i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup match type\n",                             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                 i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup\n",             i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n", i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): id of bridge back-matched to by %sgroup\n",                 i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): snapshot of bridge matched to by %sgroup\n",                i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): index    of bridge back-matched to by %sgroup\n",              i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's back-matched bridge\n",i_column++,group_text_prefix);
               fprintf(fp,"# (%02d): number of particles in the %sgroup's main progenitor\n",    i_column++,group_text_prefix);
            }
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_EMERGED)){
               if(!check_mode_for_flag(halos[i_halo].type,TREE_CASE_FOUND)){
                  fp=fp_emerged_unfound_out;
                  fprintf(fp,"%10.3le %10.3le %5d %7d %7d %8d %3d %8d %8d %8d %7d %4d %8d %8d\n",
                     a_list[l_write],
                     dt_descendant,
                     j_write,
                     i_halo,
                     halos[i_halo].id,
                     halos[i_halo].type,
                     halos[i_halo].descendant.file-i_write,
                     halos[i_halo].n_particles,
                     halos[i_halo].descendant.n_particles,
                     n_particles_parent[i_halo],
                     halos[i_halo].bridge_backmatch->id,
                     halos[i_halo].bridge_backmatch->file,
                     halos[i_halo].bridge_backmatch->index,
                     halos[i_halo].bridge_backmatch->n_particles);
               }
               else{
                  fp=fp_emerged_found_out;
                  fprintf(fp,"%10.3le %10.3le %5d %7d %7d %8d %3d %8d %8d %8d %7d %4d %8d %8d %8d\n",
                     a_list[l_write],
                     dt_descendant,
                     j_write,
                     i_halo,
                     halos[i_halo].id,
                     halos[i_halo].type,
                     halos[i_halo].descendant.file-i_write,
                     halos[i_halo].n_particles,
                     halos[i_halo].descendant.n_particles,
                     n_particles_parent[i_halo],
                     halos[i_halo].bridge_backmatch->id,
                     halos[i_halo].bridge_backmatch->file,
                     halos[i_halo].bridge_backmatch->index,
                     halos[i_halo].bridge_backmatch->n_particles,
                     halos[i_halo].size_main_progenitor);
               }
            }
         }
         fclose(fp_mergers_out);
         fclose(fp_strayed_out);
         fclose(fp_sputtered_out);
         fclose(fp_dropped_out);
         fclose(fp_bridged_out);
         fclose(fp_bridged_progenitor_out);
         fclose(fp_emerged_progenitor_out);
         fclose(fp_emerged_found_out);
         fclose(fp_emerged_unfound_out);
         SID_free(SID_FARG n_particles_parent);
      }
   }

   SID_log("Done.",SID_LOG_CLOSE);

}

void read_matches_local(char    *filename_root_matches,
                       int      i_read,
                       int      j_read,
                       int      mode,
                       int     *n_groups_i,
                       int     *n_groups_j,
                       int     *n_particles_i,
                       int     *n_particles_j,
                       int     *n_sub_group_i,
                       int     *n_sub_group_j,
                       int     *match_ids,
                       float   *match_score,
                       size_t  *match_index);
void read_matches_local(char    *filename_root_matches,
                       int      i_read,
                       int      j_read,
                       int      mode,
                       int     *n_groups_i,
                       int     *n_groups_j,
                       int     *n_particles_i,
                       int     *n_particles_j,
                       int     *n_sub_group_i,
                       int     *n_sub_group_j,
                       int     *match_ids,
                       float   *match_score,
                       size_t  *match_index){
   char   group_text_prefix[5];
   char   filename_in[MAX_FILENAME_LENGTH];
   SID_fp fp_in;
   int k_read;
   int l_read;
   int i_read_stop;
   int i_read_start;
   int n_search;
   int n_files;
   int n_matches;
   size_t offset;
   int flag_continue;
   int i_read_file;
   int j_read_file;
   int n_groups_file;
   int n_groups_file_1;
   int n_groups_file_2;
   int n_groups;
   int n_groups_i_file;
   int n_groups_j_file;
   
   switch(mode){
      case MATCH_SUBGROUPS:
      sprintf(group_text_prefix,"sub");
      break;
      case MATCH_GROUPS:
      sprintf(group_text_prefix,"");
      break;
   }
   sprintf(filename_in,"%s.%sgroup_matches",filename_root_matches,group_text_prefix);

   if(i_read==j_read)
     SID_trap_error("i_read=j_read in read_matches",ERROR_LOGIC);

   SID_fopen(filename_in,"r",&fp_in);
   SID_fread(&i_read_start,sizeof(int),1,&fp_in);
   SID_fread(&i_read_stop, sizeof(int),1,&fp_in);
   SID_fread(&n_search,    sizeof(int),1,&fp_in);
   SID_fread(&n_files,     sizeof(int),1,&fp_in);
   for(k_read=0,n_groups_i_file=-1,n_groups_j_file=-1;k_read<n_files;k_read++){
      SID_fread(&l_read,  sizeof(int),1,&fp_in);
      SID_fread(&n_groups,sizeof(int),1,&fp_in);
      if(i_read==l_read){
         n_groups_i_file=n_groups;
         SID_fread(n_particles_i,sizeof(int),n_groups,&fp_in);
         if(mode==MATCH_GROUPS){
            if(n_sub_group_i!=NULL)
               SID_fread(n_sub_group_i,sizeof(int),n_groups,&fp_in);
            else
               SID_fseek(&fp_in,sizeof(int),n_groups,SID_SEEK_CUR);
         }
      }
      else if(j_read==l_read){
         n_groups_j_file=n_groups;
         SID_fread(n_particles_j,sizeof(int),n_groups,&fp_in);
         if(mode==MATCH_GROUPS){
            if(n_sub_group_j!=NULL)
               SID_fread(n_sub_group_j,sizeof(int),n_groups,&fp_in);
            else
               SID_fseek(&fp_in,sizeof(int),n_groups,SID_SEEK_CUR);
         }
      }
      else{
         SID_fseek(&fp_in,sizeof(int),n_groups,SID_SEEK_CUR);
         if(mode==MATCH_GROUPS)
            SID_fseek(&fp_in,sizeof(int),n_groups,SID_SEEK_CUR);
      }
   }

   // Sanity check
   if(n_groups_i_file<0) SID_trap_error("File 1 (%d) group count not validly set (%d)",ERROR_LOGIC,i_read,n_groups_i_file);
   if(n_groups_j_file<0) SID_trap_error("File 2 (%d) group count not validly set (%d)",ERROR_LOGIC,j_read,n_groups_j_file);

   // Find the match we are looking for and read the offset to it's matching data
   SID_fread(&n_matches,   sizeof(int),1,&fp_in);
   for(k_read=0,offset=0,flag_continue=TRUE;k_read<n_matches && flag_continue;k_read++){
      SID_fread(&i_read_file,sizeof(int),   1,&fp_in);
      SID_fread(&j_read_file,sizeof(int),   1,&fp_in);
      SID_fread(&offset,     sizeof(size_t),1,&fp_in);
      if(flag_continue){
         if(i_read_file==i_read && j_read_file==j_read){
            (*n_groups_i)=n_groups_i_file;
            (*n_groups_j)=n_groups_j_file;
            flag_continue=FALSE;
         }
      }
   }

   // Sanity check
   if(flag_continue)
      SID_trap_error("Requested matching combination (%d->%d) not present in the matching file.",ERROR_LOGIC,i_read,j_read);

   // Offset to the matching data
   SID_fseek(&fp_in,1,offset,SID_SEEK_SET);

   // Sanity check
   SID_fread(&i_read_file,    sizeof(int),1,&fp_in);
   SID_fread(&j_read_file,    sizeof(int),1,&fp_in);
   SID_fread(&n_groups_file_1,sizeof(int),1,&fp_in);
   SID_fread(&n_groups_file_2,sizeof(int),1,&fp_in);
   if(i_read_file==i_read && j_read_file==j_read){
      // Read matching data
      SID_fread(match_ids,  sizeof(int),   (*n_groups_i),&fp_in);
      SID_fread(match_index,sizeof(size_t),(*n_groups_i),&fp_in);
      SID_fread(match_score,sizeof(float), (*n_groups_i),&fp_in);
      SID_fclose(&fp_in);
   }
   else{
      SID_fclose(&fp_in);
      SID_trap_error("Error encountered in the matching file",ERROR_LOGIC);      
   }
}

// Set halo_i[i_halo] so it points to halo_j[j_halo]
void set_halo_and_descendant(tree_horizontal_info **halos,
                             int                    i_file,
                             int                    i_halo,
                             int                    j_file,
                             int                    j_halo,
                             float                  score,
                             int                   *max_id,
                             int                    n_search);
void set_halo_and_descendant(tree_horizontal_info **halos,
                             int                    i_file,
                             int                    i_halo,
                             int                    j_file,
                             int                    j_halo,
                             float                  score,
                             int                   *max_id,
                             int                    n_search){
   tree_horizontal_info *halos_i;
   tree_horizontal_info *halos_j;
   int                   file_offset;
   int                   k_file;
   int                   k_index;
   int                   k_file_temp;
   int                   k_index_temp;
   int                   k_file_main;
   int                   k_index_main;
   int                   k_size_main;
   int                   flag_process;
   int                   n_p_diff_old;
   int                   n_p_diff_new;

   // Process the inputs a bit
   halos_i    =halos[i_file%n_search];
   halos_j    =halos[j_file%n_search];
   file_offset=j_file-i_file;
   if(file_offset==0)
      SID_trap_error("A zero file offset has been requested.  It should be -ve for roots and +ve otherwise.",ERROR_LOGIC);
      
   // Set non-bridged halos or finalize bridge matches (ie. set defaults for bridge progenitors not matched to emerged halos)
   if(!check_mode_for_flag(halos_j[j_halo].type,TREE_CASE_BRIDGED)                       || 
       check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED) ||
       check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_FINALIZE)){

      // If we are processing a bridge progenitor, only accept
      //   an emerged halo match if it is a good one ...
      flag_process=TRUE;
      if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED) && 
         !check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_FINALIZE)){
        // ... the score must be good ...
        if(halos_i[i_halo].bridge_match->score>0.75*score)
          flag_process=FALSE;
        // ... the change in halo size must be better ...
        n_p_diff_old=IABS(halos_i[i_halo].bridge_match->n_particles-halos_i[i_halo].n_particles);
        n_p_diff_new=IABS(halos_j[j_halo].n_particles-halos_i[i_halo].n_particles);
        if(n_p_diff_new>=n_p_diff_old)
          flag_process=FALSE;
        // ... it must not be the bridged halo at an earlier snapshot ...
        if(halos_i[i_halo].bridge_match->id==halos_j[j_halo].id)
          flag_process=FALSE;
         k_file =j_file;
         k_index=j_halo;
         while(k_file>=0 && k_index>=0 && flag_process){
             if(k_file==halos_i[i_halo].bridge_match->file && k_index==halos_i[i_halo].bridge_match->index)
               flag_process=FALSE;
             k_file_temp =k_file;
             k_index_temp=k_index;
             k_file =halos[k_file_temp%n_search][k_index_temp].descendant.file;
             k_index=halos[k_file_temp%n_search][k_index_temp].descendant.index;
         }
      }
      if(flag_process){

         // Set halo and tree id
         halos_i[i_halo].tree_id=halos_j[j_halo].tree_id;
         if(halos_j[j_halo].n_progenitors>0){
            // Don't create a new ID if we are matching to a
            //   strayed or sputtered halo.
            if(halos_j[j_halo].id>=0)
              halos_i[i_halo].id=(*max_id)++;
            else
              halos_i[i_halo].id=halos_j[j_halo].id;
         }
         else
            halos_i[i_halo].id=halos_j[j_halo].id;

         // Set progenitor info 
         halos_j[j_halo].n_progenitors++;
         if(halos_j[j_halo].n_progenitors==1){
            halos_j[j_halo].file_first_progenitor =i_file;
            halos_j[j_halo].index_first_progenitor=i_halo;
         }
         else{
            halos[(halos_j[j_halo].file_last_progenitor)%n_search][(halos_j[j_halo].index_last_progenitor)].file_next_progenitor =i_file;
            halos[(halos_j[j_halo].file_last_progenitor)%n_search][(halos_j[j_halo].index_last_progenitor)].index_next_progenitor=i_halo;
         }
         halos_j[j_halo].file_last_progenitor =i_file;
         halos_j[j_halo].index_last_progenitor=i_halo;

         // Determine the main progenitor
         k_file      = halos_j[j_halo].file_first_progenitor;
         k_index     = halos_j[j_halo].index_first_progenitor;
         k_file_main = k_file;
         k_index_main= k_index;
         k_size_main = halos[k_file_main%n_search][k_index_main].n_particles;
         k_file_temp =k_file;
         k_index_temp=k_index;
         k_file =halos[k_file_temp%n_search][k_index_temp].file_next_progenitor;
         k_index=halos[k_file_temp%n_search][k_index_temp].index_next_progenitor;
         while(k_file>=0 && k_index>=0){
           if(halos[k_file%n_search][k_index].n_particles>k_size_main){
             // Switch IDs so that the correct halo inherits the main progenitor's ID
             //    This change in IDs needs to be propagated all the way down both halo's trees 
             switch_horizontal_IDs_recursive(halos,n_search,j_file,j_halo,halos[k_file%n_search][k_index].id,halos[j_file%n_search][j_halo].id);
             k_file_main =k_file;
             k_index_main=k_index;
             k_size_main =halos[k_file_main%n_search][k_index_main].n_particles;
           }
           k_file_temp =k_file;
           k_index_temp=k_index;
           k_file =halos[k_file_temp%n_search][k_index_temp].file_next_progenitor;
           k_index=halos[k_file_temp%n_search][k_index_temp].index_next_progenitor;
         }

         // Initialize all the progenitors to be mergers and set main progenitor size
         k_file =halos_j[j_halo].file_first_progenitor;
         k_index=halos_j[j_halo].index_first_progenitor;
         while(k_file>=0 && k_index>=0){
           halos[k_file%n_search][k_index].type                |=TREE_CASE_MERGER;
           halos[k_file%n_search][k_index].size_main_progenitor =k_size_main;
           k_file_temp =k_file;
           k_index_temp=k_index;
           k_file =halos[k_file_temp%n_search][k_index_temp].file_next_progenitor;
           k_index=halos[k_file_temp%n_search][k_index_temp].index_next_progenitor;
         }

         // Set type for main progenitor (make sure it's merger flag is switched off too)
         halos[k_file_main%n_search][k_index_main].type|=TREE_CASE_MAIN_PROGENITOR;
         halos[k_file_main%n_search][k_index_main].type&=(~TREE_CASE_MERGER);

         // Set descendant info
         halos_i[i_halo].descendant.id         =halos_j[j_halo].id;
         halos_i[i_halo].descendant.file       =j_file;
         halos_i[i_halo].descendant.index      =j_halo;
         halos_i[i_halo].descendant.score      =score;
         halos_i[i_halo].descendant.n_particles=halos_j[j_halo].n_particles;

         // Set match-type flags ...                   
         //   If we've matched to a halo with a well-defined id then make sure it isn't marked as sputtering ...
         if(halos_i[i_halo].id>=0)
            halos_i[i_halo].type=halos_i[i_halo].type&=(~TREE_CASE_SPUTTERED);

         // ... else we've matched to a halo which *does not* have a valid id.  Mark it as
         //     a sputtering halo.
         else
            halos_i[i_halo].type|=TREE_CASE_SPUTTERED;
         
         // Set flag for simple matches and dropped halos
         if(file_offset==1)
            halos_i[i_halo].type|=TREE_CASE_SIMPLE;
         else{
            halos_i[i_halo].type|=TREE_CASE_DROPPED;
            halos_j[j_halo].type|=TREE_CASE_FOUND;
         }

         // Mark the halo as processed
         halos_i[i_halo].type&=(~(TREE_CASE_UNPROCESSED|TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED|TREE_CASE_BRIDGE_FINALIZE));
      }
   }
   // ... else we've matched to a bridge.  Here we just set the info needed
   //     to connect this halo with the bridge it's been matched to.  We'll
   //     use that info to search the bridge's emergent halos and if we fail
   //     to identify to one of those, we'll finalize this match as the default.
   else{
      if(halos_i[i_halo].bridge_match==NULL){
         halos_i[i_halo].bridge_match             =(match_info *)SID_calloc(sizeof(match_info));
         halos_i[i_halo].bridge_match->id         =halos_j[j_halo].id;
         halos_i[i_halo].bridge_match->file       =j_file;
         halos_i[i_halo].bridge_match->score      =score;
         halos_i[i_halo].bridge_match->n_particles=halos_j[j_halo].n_particles;
         halos_i[i_halo].bridge_match->index      =j_halo;
         halos_i[i_halo].type|=(TREE_CASE_BRIDGE_PROGENITOR|TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED);
         halos_i[i_halo].type&=(~TREE_CASE_UNPROCESSED);
      }
   }
}

void compute_trees_horizontal(char   *filename_halo_root_in,
                              char   *filename_cat_root_in,
                              char   *filename_root_matches,
                              char   *filename_output_dir,
                              double *a_list,
                              cosmo_info **cosmo,
                              int     i_read_start,
                              int     i_read_stop,
                              int     i_read_step,
                              int     n_search,
                              int    *flag_clean){
  char        group_text_prefix[5];
  FILE       *fp;
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
  int         l_write;
  int         j_file_1;
  int         j_file_2;
  int         i_read;
  int         j_read;
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
  int        *my_descendant;
  int       **n_particles;
  int       **n_particles_groups;
  int       **n_particles_subgroups;
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
  int         i_find,n_find;
  int         flag_continue;
  int         flag_drop;
  int        *match_id=NULL;
  int        *search_id=NULL;
  int         n_progenitors_max;
  int         i_search;
  int         flag_dropped;
  int         flag_first;
  int         n_particles_max;
  int         trunk_index;
  int         biggest_stray;
  int         biggest_stray_drop;
  int         biggest_stray_bridge;
  int         biggest_sputter;
  int        *n_groups=NULL;
  int        *n_subgroups=NULL;
  int         max_tree_id_group;
  int         max_tree_id_subgroup;
  int         max_tree_id;
  int       **n_subgroups_group=NULL;
  int        *n_subgroups_group_1=NULL;
  size_t    **sort_id=NULL;
  size_t    **sort_group_id=NULL;
  size_t    **sort_subgroup_id=NULL;
  size_t     *match_index=NULL;
  size_t     *bridge_index=NULL;
  size_t     *search_index=NULL;
  float      *match_score=NULL;
  int        *bridge_keep=NULL;
  int         flag_match_subgroups;
  int         flag_keep_strays=FALSE;
  int         n_k_match=2;
  int         n_snap;
  
  tree_horizontal_info **subgroups;
  tree_horizontal_info **groups;
  tree_horizontal_info **halos;
  tree_horizontal_info  *halos_i;
  match_info            *bridges;
  match_info            *bridge;

  int  n_files;
  int  n_subgroups_max;
  int  n_groups_max;
  int *n_halos;
  int  n_halos_max;
  int  n_halos_i;
  int  i_halo;
  int      n_halos_1_matches;
  int      n_halos_2_matches;
  int     j_halo;
  int     k_halo;
  int     l_halo;

  int     n_list;
  int     k_file;
  int     k_index;
  int     k_file_temp;
  int     k_index_temp;
  int     main_progenitor_index;
  float   main_progenitor_score;
  
  int     n_strayed;
  int     n_sputtered;
  int     n_bridged;
  int     n_dropped;
  int     n_bridged_systems;
  int     max_strayed_size;
  int     max_sputtered_size;
  int     max_dropped_size;
  int     max_emerged_size;
  int     n_bridge_emerged;
  
  tree_horizontal_stats_info stats;
  
  SID_log("Constructing horizontal merger trees for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_stop,i_read_start,i_read_step);

  if(n_search<1)
    SID_trap_error("n_search=%d but must be at least 1",ERROR_LOGIC);

  // Validate existing matching files &/or perfrom matching
  compute_trees_matches(filename_halo_root_in,
                        filename_root_matches,
                        i_read_start,
                        i_read_stop,
                        i_read_step,
                        &n_files,
                        &n_subgroups,
                        &n_groups,
                        n_search);

  // We need these for allocating arrays
  calc_max(n_subgroups,&n_subgroups_max,n_files,SID_INT,CALC_MODE_DEFAULT);
  calc_max(n_groups,   &n_groups_max,   n_files,SID_INT,CALC_MODE_DEFAULT);
  n_halos_max=MAX(n_subgroups_max,n_groups_max);

  // We need indices for the current and last i_file as well
  n_search+=2; 
     
  // Initialize arrays
  SID_log("Creating arrays...",SID_LOG_OPEN);
  n_particles_groups   =(int   **)SID_malloc(sizeof(int *) *n_search);
  n_particles_subgroups=(int   **)SID_malloc(sizeof(int *) *n_search);
  match_id             =(int    *)SID_malloc(sizeof(int)   *n_halos_max);
  match_score          =(float  *)SID_malloc(sizeof(float) *n_halos_max);
  match_index          =(size_t *)SID_malloc(sizeof(size_t)*n_halos_max);
  subgroups            =(tree_horizontal_info **)SID_malloc(sizeof(tree_horizontal_info *)*n_search);
  groups               =(tree_horizontal_info **)SID_malloc(sizeof(tree_horizontal_info *)*n_search);
  n_subgroups_group    =(int                  **)SID_malloc(sizeof(int                  *)*n_search);
  for(i_search=0;i_search<n_search;i_search++){
     n_particles_groups[i_search]   =(int                  *)SID_malloc(sizeof(int)                 *n_groups_max);       
     n_particles_subgroups[i_search]=(int                  *)SID_malloc(sizeof(int)                 *n_subgroups_max);       
     subgroups[i_search]            =(tree_horizontal_info *)SID_calloc(sizeof(tree_horizontal_info)*n_subgroups_max);
     groups[i_search]               =(tree_horizontal_info *)SID_calloc(sizeof(tree_horizontal_info)*n_groups_max);       
     n_subgroups_group[i_search]    =(int                  *)SID_malloc(sizeof(int)                 *n_groups_max);       
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Process the first file separately
  //   (just give everything ids from a running index) ...
  SID_log("Initializing tree roots...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Initialize everything to a 1:1 simple match
  read_matches_local(filename_root_matches,
                     i_read_stop,i_read_stop-i_read_step,
                     MATCH_SUBGROUPS,
                     &n_halos_1_matches,
                     &n_halos_2_matches,
                     n_particles_subgroups[0],
                     n_particles_subgroups[1],
                     NULL,
                     NULL,
                     match_id,
                     match_score,
                     match_index);
  for(i_search=1;i_search<n_search;i_search++)
    memcpy(n_particles_subgroups[i_search],n_particles_subgroups[0],n_halos_1_matches*sizeof(int));
                     
  for(i_halo=0,max_id_subgroup=0,max_tree_id_subgroup=0;i_halo<n_subgroups_max;i_halo++){
     for(i_search=0;i_search<n_search;i_search++){
        subgroups[i_search][i_halo].id              =-1;
        subgroups[i_search][i_halo].tree_id         =-1; // the resulting file offset must be -ve for roots 
        subgroups[i_search][i_halo].descendant.id   =-1;
        subgroups[i_search][i_halo].descendant.file =-1;
        subgroups[i_search][i_halo].descendant.index=-1;
        subgroups[i_search][i_halo].descendant.score= 0.;
        subgroups[i_search][i_halo].descendant.n_particles=0;
        subgroups[i_search][i_halo].n_progenitors         = 0;
        subgroups[i_search][i_halo].file_first_progenitor =-1;
        subgroups[i_search][i_halo].index_first_progenitor=-1;
        subgroups[i_search][i_halo].file_last_progenitor  =-1;
        subgroups[i_search][i_halo].index_last_progenitor =-1;
        subgroups[i_search][i_halo].file_next_progenitor  =-1;
        subgroups[i_search][i_halo].index_next_progenitor =-1;
        subgroups[i_search][i_halo].size_main_progenitor  = 0;
        subgroups[i_search][i_halo].n_particles     = 0;
        subgroups[i_search][i_halo].type            = TREE_CASE_INVALID;
        subgroups[i_search][i_halo].n_bridges       = 0;
        subgroups[i_search][i_halo].bridges         = NULL;
        subgroups[i_search][i_halo].bridge_match    = NULL;
        subgroups[i_search][i_halo].bridge_backmatch= NULL;
        if(i_halo<n_halos_1_matches){
           subgroups[i_search][i_halo].id                    =max_id_subgroup++;
           subgroups[i_search][i_halo].tree_id               =max_tree_id_subgroup++;
           subgroups[i_search][i_halo].descendant.id         =subgroups[i_search][i_halo].id;
           subgroups[i_search][i_halo].descendant.file       =i_read_stop;
           subgroups[i_search][i_halo].descendant.index      =i_halo;
           subgroups[i_search][i_halo].descendant.score      = 0.;
           subgroups[i_search][i_halo].descendant.n_particles=n_particles_subgroups[i_search][i_halo];
           subgroups[i_search][i_halo].n_particles           =n_particles_subgroups[i_search][i_halo];
           subgroups[i_search][i_halo].type                  =TREE_CASE_SIMPLE|TREE_CASE_MAIN_PROGENITOR;
        }
     }
     subgroups[i_read_stop%n_search][i_halo].n_progenitors=0;
  }

  // Initialize everything to a 1:1 simple match
  read_matches_local(filename_root_matches,
                     i_read_stop,i_read_stop-i_read_step,
                     MATCH_GROUPS,
                     &n_halos_1_matches,
                     &n_halos_2_matches,
                     n_particles_groups[0],
                     n_particles_groups[1],
                     n_subgroups_group[0],
                     n_subgroups_group[1],
                     match_id,
                     match_score,
                     match_index);
  for(i_search=1;i_search<n_search;i_search++){
    memcpy(n_particles_groups[i_search],n_particles_groups[0],n_halos_1_matches*sizeof(int));
    memcpy(n_subgroups_group[i_search], n_subgroups_group[0], n_halos_1_matches*sizeof(int));
  }

  for(i_halo=0,max_id_group=0,max_tree_id_group=0;i_halo<n_groups_max;i_halo++){
     for(i_search=0;i_search<n_search;i_search++){
        groups[i_search][i_halo].id              =-1;
        groups[i_search][i_halo].tree_id         =-1;
        groups[i_search][i_halo].descendant.id   =-1;
        groups[i_search][i_halo].descendant.file =-1; // the resulting file offset must be -ve for roots
        groups[i_search][i_halo].descendant.index=-1;
        groups[i_search][i_halo].descendant.n_particles=0;
        groups[i_search][i_halo].descendant.score     = 0.;
        groups[i_search][i_halo].n_progenitors        = 0;
        groups[i_search][i_halo].file_first_progenitor =-1;
        groups[i_search][i_halo].index_first_progenitor=-1;
        groups[i_search][i_halo].file_last_progenitor  =-1;
        groups[i_search][i_halo].index_last_progenitor =-1;
        groups[i_search][i_halo].file_next_progenitor  =-1;
        groups[i_search][i_halo].index_next_progenitor =-1;
        groups[i_search][i_halo].size_main_progenitor  = 0;
        groups[i_search][i_halo].n_particles     = 0;
        groups[i_search][i_halo].type            = TREE_CASE_INVALID;
        groups[i_search][i_halo].n_bridges       = 0;
        groups[i_search][i_halo].bridges         = NULL;
        groups[i_search][i_halo].bridge_match    = NULL;
        groups[i_search][i_halo].bridge_backmatch= NULL;
        if(i_halo<n_halos_1_matches){
           groups[i_search][i_halo].id                    =max_id_group;
           groups[i_search][i_halo].tree_id               =max_tree_id_group;
           groups[i_search][i_halo].descendant.id         =groups[i_search][i_halo].id;
           groups[i_search][i_halo].descendant.file       =i_read_stop;
           groups[i_search][i_halo].descendant.index      =i_halo;
           groups[i_search][i_halo].descendant.score      = 0.;
           groups[i_search][i_halo].descendant.n_particles=n_particles_groups[i_search][i_halo];
           groups[i_search][i_halo].n_particles           =n_particles_groups[i_search][i_halo];
           groups[i_search][i_halo].type                  =TREE_CASE_SIMPLE|TREE_CASE_MAIN_PROGENITOR;
           max_id_group++;
           max_tree_id_group++;
        }
     }
     groups[i_read_stop%n_search][i_halo].n_progenitors=0;
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // The first snapshot is done now (set to defaults as the roots of trees) ... now loop over all other snapshots ...
  //   There are a bunch of counters at work here.  Because we aren't necessarily using every 
  //     snapshot (if i_read_step>1), we need counters to keep track of which snapshots we
  //     are working with (i_read_*,j_read_*, etc), counters to keep track of which
  //     files's we're dealing with as far as the trees are concerned (i_file_*,j_file_*,etc), and
  //     counters to keep track of which files are being/have been written (i_write_*,j_write_* etc).
  //     We can't write files right away because previously processed snapshots can be changed
  //     when we deal with dropped and bridged halos.
  for(i_read   =i_read_stop-i_read_step,
        i_file =i_read_stop-1, 
        j_file =1,             
        i_write=i_read_stop,      
        j_write=i_read_stop,
        l_write=0;      
      i_read>=i_read_start;
      i_read-=i_read_step,    
         i_file--, 
         j_file++){   

    SID_log("Processing snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read);

    // Loop twice (1st to process subgroups, 2nd to process groups)
    for(k_match=0;k_match<n_k_match;k_match++){

       // Initialize a bunch of stuff depending on whether
       //   we are processing groups or subgroups
       switch(k_match){
          case 0:
          sprintf(group_text_prefix,"sub");
          flag_match_subgroups=MATCH_SUBGROUPS;
          halos               =subgroups;
          n_halos             =n_subgroups;
          n_halos_max         =n_subgroups_max;
          max_id              =max_id_subgroup;
          n_particles         =n_particles_subgroups;
          break;
          case 1:
          sprintf(group_text_prefix,"");
          flag_match_subgroups=MATCH_GROUPS;
          halos               =groups;
          n_halos             =n_groups;
          n_halos_max         =n_groups_max;
          max_id              =max_id_group;
          n_particles         =n_particles_groups;
          break;
       }
       halos_i  =halos[i_file%n_search];
       n_halos_i=n_halos[j_file];
       SID_log("Processing %d %sgroups...",SID_LOG_OPEN|SID_LOG_TIMER,n_halos_i,group_text_prefix);

       // Initialize tree pointer-arrays with unique and negative dummy values
       for(i_halo=0;i_halo<n_halos_max;i_halo++){
          halos_i[i_halo].id                    =-1;
          halos_i[i_halo].tree_id               =-1;
          halos_i[i_halo].descendant.id         =-1;
          halos_i[i_halo].descendant.file       =-1; // the resulting file offset must be -ve for roots
          halos_i[i_halo].descendant.index      =-1;
          halos_i[i_halo].descendant.n_particles= 0;
          halos_i[i_halo].descendant.score      = 0.;
          halos_i[i_halo].n_progenitors         = 0;
          halos_i[i_halo].file_first_progenitor =-1;
          halos_i[i_halo].index_first_progenitor=-1;
          halos_i[i_halo].file_last_progenitor  =-1;
          halos_i[i_halo].index_last_progenitor =-1;
          halos_i[i_halo].file_next_progenitor  =-1;
          halos_i[i_halo].index_next_progenitor =-1;
          halos_i[i_halo].size_main_progenitor  = 0;
          halos_i[i_halo].n_particles           = 0;
          halos_i[i_halo].n_bridges             = 0;
          SID_free(SID_FARG halos_i[i_halo].bridges);
          SID_free(SID_FARG halos_i[i_halo].bridge_match);
          SID_free(SID_FARG halos_i[i_halo].bridge_backmatch);
          if(i_halo<n_halos_i)
             halos_i[i_halo].type=TREE_CASE_UNPROCESSED;
          else
             halos_i[i_halo].type=TREE_CASE_INVALID;
       }
       
       // Use back-matching to identify bridged halos ...
       SID_log("Identifying bridge candidates from back-matching...",SID_LOG_OPEN|SID_LOG_TIMER);
       SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
       //    ... first, do an initial count of matches.  This will not be a list of unique halos
       //        though, since the same halos are likely to appear in repeated snapshots.
       for(j_file_1  =i_file+1,
             j_file_2=i_file,
             j_read_1=i_read+i_read_step,
             j_read_2=i_read,
             i_search=0;
           j_read_1<=i_read_stop && i_search<(n_search-2);
           j_file_1++,
             j_read_1+=i_read_step,
             i_search++){

          SID_log("Counting matches between files %d->%d...",SID_LOG_OPEN,j_read_1,j_read_2);

          // Read back-matching
          read_matches_local(filename_root_matches,
                             j_read_1,j_read_2,
                             flag_match_subgroups,
                             &n_halos_1_matches,
                             &n_halos_2_matches,
                             n_particles[j_file_1%n_search],
                             n_particles[j_file_2%n_search],
                             NULL,
                             NULL,
                             match_id,
                             match_score,
                             match_index);

          // Store halo sizes
          if(i_search==0){
             for(i_halo=0;i_halo<n_halos_2_matches;i_halo++)
                halos[j_read_2%n_search][i_halo].n_particles=n_particles[j_file_2%n_search][i_halo];
          }

          // Perform initial back-match count
          for(i_halo=0;i_halo<n_halos_i;i_halo++){
             j_halo=find_index_int(match_id,i_halo,n_halos_1_matches,match_index);
             while(match_id[match_index[j_halo]]==i_halo && j_halo<(n_halos_1_matches-1)){
                halos_i[i_halo].n_bridges++;
                j_halo++;
             }
             if(match_id[match_index[j_halo]]==i_halo && j_halo==(n_halos_1_matches-1)){
                halos_i[i_halo].n_bridges++;
                j_halo++;
             }
          }
          SID_log("Done.",SID_LOG_CLOSE);
       }
       
       //    ... second, do a conservative allocation using the non-unique counts and reset the counter.
       for(i_halo=0;i_halo<n_halos_i;i_halo++){
          if((halos_i[i_halo].n_bridges)>0)
             (halos_i[i_halo].bridges)=(match_info *)SID_calloc(sizeof(match_info)*(halos_i[i_halo].n_bridges));
          else
             (halos_i[i_halo].bridges)=NULL;
          halos_i[i_halo].n_bridges =0;
       }

       //    ... third, assemble the list of unique halos.
       for(j_file_1  =i_file+1,
             j_file_2=i_file,
             j_read_1=i_read+i_read_step,
             j_read_2=i_read,
             i_search=0;
           j_read_1<=i_read_stop && i_search<(n_search-2);
           j_file_1++,
             j_read_1+=i_read_step,
             i_search++){

          SID_log("Finding unique matches between files %d->%d...",SID_LOG_OPEN,j_read_1,j_read_2);
          
          // Read back-matching
          read_matches_local(filename_root_matches,
                             j_read_1,j_read_2,
                             flag_match_subgroups,
                             &n_halos_1_matches,
                             &n_halos_2_matches,
                             n_particles[j_file_1%n_search],
                             n_particles[j_file_2%n_search],
                             NULL,
                             NULL,
                             match_id,
                             match_score,
                             match_index);
                             
          // For all the halos in i_file_1 with back-matches ...
          for(i_halo=0;i_halo<n_halos_i;i_halo++){
             if((halos_i[i_halo].bridges)!=NULL){
                bridges=halos_i[i_halo].bridges;
                j_halo =find_index_int(match_id,i_halo,n_halos_1_matches,match_index);
                while(match_id[match_index[j_halo]]==i_halo && j_halo<(n_halos_1_matches-1)){
                   // Check to see if this halo is already in the list ...
                   for(k_halo=0,flag_continue=TRUE;k_halo<halos_i[i_halo].n_bridges && flag_continue;k_halo++){
                      if(bridges[k_halo].id==halos[j_file_1%n_search][match_index[j_halo]].id)
                         flag_continue=FALSE;
                   }
                   // ... if not, add it
                   if(flag_continue){
                      bridges[halos_i[i_halo].n_bridges].id         =halos[j_file_1%n_search][match_index[j_halo]].id;
                      bridges[halos_i[i_halo].n_bridges].file       =j_file_1;
                      bridges[halos_i[i_halo].n_bridges].index      =match_index[j_halo];
                      bridges[halos_i[i_halo].n_bridges].score      =match_score[match_index[j_halo]];
                      bridges[halos_i[i_halo].n_bridges].n_particles=halos[j_file_1%n_search][match_index[j_halo]].n_particles;
                      (halos_i[i_halo].n_bridges)++;
                   }
                   j_halo++;
                }
                if(match_id[match_index[j_halo]]==i_halo && j_halo==(n_halos_1_matches-1)){
                   // Check to see if this halo is already in the list ...
                   for(k_halo=0,flag_continue=TRUE;k_halo<halos_i[i_halo].n_bridges && flag_continue;k_halo++){
                      if(bridges[k_halo].id==halos[j_file_1%n_search][match_index[j_halo]].id)
                         flag_continue=FALSE;
                   }
                   // ... if not, add it
                   if(flag_continue){
                      bridges[halos_i[i_halo].n_bridges].id         =halos[j_file_1%n_search][match_index[j_halo]].id;
                      bridges[halos_i[i_halo].n_bridges].file       =j_file_1;
                      bridges[halos_i[i_halo].n_bridges].index      =match_index[j_halo];
                      bridges[halos_i[i_halo].n_bridges].score      =match_score[match_index[j_halo]];
                      bridges[halos_i[i_halo].n_bridges].n_particles=halos[j_file_1%n_search][match_index[j_halo]].n_particles;
                      (halos_i[i_halo].n_bridges)++;
                   }
                   j_halo++;
                }
             }
          }
          SID_log("Done.",SID_LOG_CLOSE);
       }

       // ... lastly, remove bridge descendants and finalize the list ...
       for(i_halo=0;i_halo<n_halos_i;i_halo++){
          if((halos_i[i_halo].n_bridges)>1){ 

             // We may need to remove several halos from the list.  This array will keep track of this.
             bridge_keep=(int *)SID_malloc(sizeof(int)*halos_i[i_halo].n_bridges);

             // Reorder the bridges by their score/size.  We make a temporary copy of the list to do this.
             bridges=(match_info *)SID_calloc(sizeof(match_info)*(halos_i[i_halo].n_bridges));
             for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;j_halo++){
                bridge=&(halos_i[i_halo].bridges[j_halo]);
                memcpy(&(bridges[j_halo]),bridge,sizeof(match_info));
                match_score[j_halo]=bridge->score;
                bridge_keep[j_halo]=TRUE;
             }
             merge_sort((void *)match_score,(size_t)(halos_i[i_halo].n_bridges),&bridge_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

             // Remove any mutual descendants from the list
             //   (since they have their own IDs, this is 
             //    needed to avoid calling them emerged halos)
             for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;j_halo++){
                bridge  = &(bridges[bridge_index[j_halo]]);
                k_file  = halos[bridge->file%n_search][bridge->index].descendant.file;
                k_index = halos[bridge->file%n_search][bridge->index].descendant.index;
                while(k_file>=0 && k_index>=0 && k_file<MIN(i_read_stop,i_file+n_search)){
                   for(k_halo=0;k_halo<halos_i[i_halo].n_bridges;k_halo++){
                      if(k_halo!=j_halo && bridge_keep[k_halo]){
                         bridge  = &(bridges[bridge_index[k_halo]]);
                         if(bridge->file==k_file && bridge->index==k_index)
                            bridge_keep[j_halo]=FALSE;
                      }
                   }
                   k_file_temp =k_file;
                   k_index_temp=k_index;
                   k_file      =halos[k_file_temp%n_search][k_index_temp].descendant.file;
                   k_index     =halos[k_file_temp%n_search][k_index_temp].descendant.index;
                }
             }

             // Since we may have trimmed the list, recount the number remaining
             n_list=halos_i[i_halo].n_bridges;
             for(j_halo=n_list-1,halos_i[i_halo].n_bridges=0;j_halo>=0;j_halo--){
               if(bridge_keep[j_halo])
                  halos_i[i_halo].n_bridges++;
             }

             // We've removed some halos and may not actually be a bridged halo anymore.  Clean-up if so.
             if(halos_i[i_halo].n_bridges<1){
                halos_i[i_halo].type&=(~TREE_CASE_BRIDGED);
                SID_free(SID_FARG halos_i[i_halo].bridges);
                halos_i[i_halo].n_bridges=0;
             }
             else{
                halos_i[i_halo].type|=TREE_CASE_BRIDGED;
                n_bridged++;

                // Because we've overallocated previously, reallocate the bridge list here to save RAM.
                SID_free(SID_FARG halos_i[i_halo].bridges);
                (halos_i[i_halo].bridges)=(match_info *)SID_calloc(sizeof(match_info)*(halos_i[i_halo].n_bridges));

                // Copy the sorted temporary list to the permanent list,  leaving out the main progenitor.
                for(j_halo=n_list-1,l_halo=0;j_halo>=0;j_halo--){
                   if(bridge_keep[j_halo]){
                      memcpy(&(halos_i[i_halo].bridges[l_halo]),&(bridges[bridge_index[j_halo]]),sizeof(match_info));
                      halos[(bridges[bridge_index[j_halo]].file)%n_search][bridges[bridge_index[j_halo]].index].type|=TREE_CASE_EMERGED;
                      halos[(bridges[bridge_index[j_halo]].file)%n_search][bridges[bridge_index[j_halo]].index].bridge_backmatch=(match_info *)SID_calloc(sizeof(match_info));
                      halos[(bridges[bridge_index[j_halo]].file)%n_search][bridges[bridge_index[j_halo]].index].bridge_backmatch->id         =halos_i[i_halo].id;
                      halos[(bridges[bridge_index[j_halo]].file)%n_search][bridges[bridge_index[j_halo]].index].bridge_backmatch->file       =i_file;
                      halos[(bridges[bridge_index[j_halo]].file)%n_search][bridges[bridge_index[j_halo]].index].bridge_backmatch->score      =match_score[bridge_index[j_halo]];
                      halos[(bridges[bridge_index[j_halo]].file)%n_search][bridges[bridge_index[j_halo]].index].bridge_backmatch->n_particles=halos_i[i_halo].n_particles;
                      halos[(bridges[bridge_index[j_halo]].file)%n_search][bridges[bridge_index[j_halo]].index].bridge_backmatch->index      =i_halo;
                      l_halo++;
                   }
                }
             }

             // Clean-up
             SID_free(SID_FARG bridge_keep);
             SID_free(SID_FARG bridge_index);
             SID_free(SID_FARG bridges);
          }
          // This halo is not a bridge.  Perform cleaning.
          else{
             halos_i[i_halo].type&=(~TREE_CASE_BRIDGED);
             SID_free(SID_FARG halos_i[i_halo].bridges);
             halos_i[i_halo].n_bridges=0;
          }
       }
       SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
       SID_log("Done.",SID_LOG_CLOSE);

       // Perform forward-matching
       SID_log("Constructing progenitors from forward-matching...",SID_LOG_OPEN|SID_LOG_TIMER);
       for(j_file_1  =i_file,
             j_file_2=i_file+1,
             j_read_1=i_read,
             j_read_2=i_read+i_read_step,
             i_search=0;
           j_read_2<=i_read_stop && i_search<(n_search-2);
           j_file_2++,
             j_read_2+=i_read_step,
             i_search++){

          // Read forward-matching
          read_matches_local(filename_root_matches,
                             j_read_1,j_read_2,
                             flag_match_subgroups,
                             &n_halos_1_matches,
                             &n_halos_2_matches,
                             n_particles[j_file_1%n_search],
                             n_particles[j_file_2%n_search],
                             n_subgroups_group[j_file_1%n_search],
                             n_subgroups_group[j_file_2%n_search],
                             match_id,
                             match_score,
                             match_index);
                       
          // Perform matching for all the halos in i_file_1.  This loop should deal completely with
          //   all simple matches and dropped halos.  It also identifies matches to bridges, which
          //   require special treatment in the loop that follows (to look for matches to emergent halos)
          //   and at the end of the loop over j_read/i_search to finalize those not associated with emerged halos.
          for(i_halo=0;i_halo<n_halos_1_matches;i_halo++){
             // If this halo hasn't been processed during earlier searches ...
             if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_UNPROCESSED)){
                my_descendant_index=match_id[i_halo];
                // If this halo has been matched to something in i_file_2 ...
                if(my_descendant_index>=0)
                   set_halo_and_descendant(halos,
                                           i_file,
                                           i_halo,
                                           j_file_2,
                                           my_descendant_index,
                                           match_score[i_halo],
                                           &max_id,
                                           n_search);
             }
          }

          // Try to match "progenitors" of bridges to the emergent halos identified in their bridge-lists
          for(i_halo=0,n_drop=0;i_halo<n_halos_1_matches;i_halo++){
             if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED)){
                // Loop over all the emergent halos identified with the bridge that i_halo has been matched to
                if(halos_i[i_halo].bridge_match==NULL)
                  SID_trap_error("Bridge match not defined during emerged halo search.",ERROR_LOGIC);
                bridges=halos[(halos_i[i_halo].bridge_match->file)%n_search][halos_i[i_halo].bridge_match->index].bridges;
                if(bridges==NULL)
                  SID_trap_error("Bridges not defined during emerged halo search.",ERROR_LOGIC);
                for(k_halo=0;k_halo<halos[(halos_i[i_halo].bridge_match->file)%n_search][halos_i[i_halo].bridge_match->index].n_bridges;k_halo++){
                   if(bridges[k_halo].file==j_file_2 && match_id[i_halo]==bridges[k_halo].index){
                      set_halo_and_descendant(halos,
                                              i_file,
                                              i_halo,
                                              bridges[k_halo].file,
                                              bridges[k_halo].index,
                                              match_score[i_halo],
                                              &max_id,
                                              n_search);
                   }
                }
             }
          }
       }
      
       // Finalize matches to bridges that haven't been associated with any emerged halos
       for(i_halo=0,n_drop=0;i_halo<n_halos_1_matches;i_halo++){
          if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED)){
             // The descendant info is already set.  Just increment it's progenitor counter and set this halo's info.
             //   Leave the target flag untouched so we can later identify which BRIDGE_PROGENITORS were found
             halos_i[i_halo].type|=(TREE_CASE_BRIDGE_FINALIZE|TREE_CASE_BRIDGE_DEFAULT);
             set_halo_and_descendant(halos,
                                     i_file,
                                     i_halo,
                                     halos_i[i_halo].bridge_match->file,
                                     halos_i[i_halo].bridge_match->index,
                                     halos_i[i_halo].bridge_match->score,
                                     &max_id,
                                     n_search);
          }
       }
      
       // Assign flags for halos not successfully processed.  They must be strays
       for(i_halo=0;i_halo<n_halos_i;i_halo++){
          if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_UNPROCESSED)){
             halos_i[i_halo].type|=TREE_CASE_STRAYED;
             halos_i[i_halo].type&=(~TREE_CASE_UNPROCESSED);
          }
       }
       SID_log("Done.",SID_LOG_CLOSE);

       // Now that we have assigned all the IDs for the halos in the active snapshot,
       //   we need to remove all of their descendants from their bridge list so that
       //   real (rather than bridged) matches to bridged halos are properly applied
       //   when the default action takes place
       for(i_halo=0;i_halo<n_halos_i;i_halo++){
          // Check all bridged halos ...
          if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGED)){
             n_list=halos_i[i_halo].n_bridges;
             // ... and check all of their descendants ...
             k_file  = halos_i[i_halo].descendant.file;
             k_index = halos_i[i_halo].descendant.index;
             while(k_file>=0 && k_index>=0 && k_file<MIN(i_read_stop,i_file+n_search)){
                // Compare against all of the halo's bridges
                for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;){
                   bridge=&(halos_i[i_halo].bridges[j_halo]);
                   // If any of them match, remove them from the list ...
                   if(bridge->file==k_file && bridge->index==k_index){
                      // ... do this by decrementing the counter
                      halos_i[i_halo].n_bridges--;
                      // ... and sliding all the halos down ...
                      for(k_halo=j_halo;k_halo<halos_i[i_halo].n_bridges;k_halo++)
                         memcpy(&(halos_i[i_halo].bridges[k_halo]),&(halos_i[i_halo].bridges[k_halo+1]),sizeof(match_info));
                   }
                   // We only need to increment the counter if
                   //   we don't find a match
                   else
                     j_halo++;
                }
                k_file_temp =k_file;
                k_index_temp=k_index;
                k_file      =halos[k_file_temp%n_search][k_index_temp].descendant.file;
                k_index     =halos[k_file_temp%n_search][k_index_temp].descendant.index;
             }

             // Since we may have removed items, we might not have a bridged halo any more.
             if(halos_i[i_halo].n_bridges<1){
                halos_i[i_halo].type&=(~TREE_CASE_BRIDGED);
                SID_free(SID_FARG halos_i[i_halo].bridges);
                halos_i[i_halo].n_bridges=0;
             }
          }

          // Print bridge info
          /*
          if(k_match==1 && halos_i[i_halo].n_bridges>1){
             if(i_halo==0)
                printf("\n");
             printf("i_file=%8d i_halo=%8d n_p_i=%8d j_file=%8d j_halo=%8d n_p_j=%8d type=%8d score=%10.4f\n",
                     i_file,i_halo,
                     halos_i[i_halo].n_particles,
                     halos_i[i_halo].descendant.file,
                     halos_i[i_halo].descendant.index,
                     halos_i[i_halo].descendant.n_particles,
                     halos_i[i_halo].type,
                     halos_i[i_halo].descendant.score);
             for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;j_halo++){
                bridge=&(halos_i[i_halo].bridges[j_halo]);
                printf("  b_file=%8d b_halo=%8d n_p_b=%8d b_type=%8d b_score=%10.4f\n",
                        bridge->file,
                        bridge->index,
                        halos[(bridge->file)%n_search][bridge->index].n_particles,
                        halos[bridge->file%n_search][bridge->index].type,bridge->score);
             }
          }
          */
       }

       // Write statistics to the log 
       //   n.b.: This is only an estimate in some cases, since subsequent snapshots may alter this snapshot.  
       //         See the written log file for the most accurate numbers.
       compute_trees_horizontal_stats(halos_i,n_halos_i,n_halos_max,&stats);
       SID_log("Results:",SID_LOG_OPEN);
       SID_log("# of halos              =%-7d",SID_LOG_COMMENT,stats.n_halos);
       SID_log("# of simple matches     =%-7d (%d mergers)",SID_LOG_COMMENT,stats.n_simple,stats.n_mergers);
       if(stats.n_strayed>0)
          SID_log("# of strayed halos      =%-7d (largest=%d particles)",SID_LOG_COMMENT,stats.n_strayed,stats.max_strayed_size);
       else
          SID_log("# of strayed halos      =%-7d",SID_LOG_COMMENT,stats.n_strayed);
       if(stats.n_sputtered>0)
          SID_log("# of sputtering halos   =%-7d (largest=%d particles)",SID_LOG_COMMENT,stats.n_sputtered,stats.max_sputtered_size);
       else
          SID_log("# of sputtering halos   =%-7d",SID_LOG_COMMENT,stats.n_sputtered);
       if(stats.n_dropped>0)
          SID_log("# of dropped halos      =%-7d (largest=%d particles)",SID_LOG_COMMENT,stats.n_dropped,stats.max_dropped_size);
       else
          SID_log("# of dropped halos      =%-7d",SID_LOG_COMMENT,stats.n_dropped);
       if(stats.n_bridged>0)
          SID_log("# of bridged halos      =%-7d (largest=%d particles)",SID_LOG_COMMENT,stats.n_bridged,stats.max_bridged_size);
       else
          SID_log("# of bridged halos      =%-7d",SID_LOG_COMMENT,stats.n_bridged);
       if(stats.n_bridge_progenitors>0)
          SID_log("# of bridge progenitors =%-7d (largest=%d particles)",SID_LOG_COMMENT,stats.n_bridge_progenitors,stats.max_bridge_progenitor_size);
       else
          SID_log("# of bridge progenitors =%-7d",SID_LOG_COMMENT,stats.n_bridge_progenitors);
       if(stats.n_emerged_progenitors>0){
          SID_log("# of emerged progenitors=%-7d (largest=%d particles; max. diff=%d/%d)",SID_LOG_COMMENT,
             stats.n_emerged_progenitors,stats.max_emerged_progenitor_size,stats.max_emerged_found_diff,stats.max_emerged_found_diff_size);
       }
       else
          SID_log("# of emerged progenitors=%-7d",SID_LOG_COMMENT,stats.n_emerged_progenitors);      
       SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);
      
       // Update things changed in this k_match iteration
       switch(k_match){
          case 0:
            max_id_subgroup=max_id;
            break;
          case 1:
            max_id_group=max_id;
            break;
       }      
       SID_log("Done.",SID_LOG_CLOSE);
    } // k_match
    
    // Write tree info once a few files have been processed
    //   and no more dropped groups need to be given ids
    if(j_file>(n_search-2)){
       write_trees_horizontal(groups,   n_groups[l_write],   n_groups_max,   n_particles_groups[i_write%n_search],
                              subgroups,n_subgroups[l_write],n_subgroups_max,n_particles_subgroups[i_write%n_search],
                              n_subgroups_group,
                              n_halos_max,
                              max_tree_id_subgroup,
                              max_tree_id_group,
                              i_write,
                              j_write,
                              l_write,
                              n_search,
                              filename_cat_root_in,
                              filename_output_dir,
                              a_list,
                              cosmo,
                              n_k_match);
       i_write--;
       l_write++;
       j_write-=i_read_step;
    }
    
    SID_log("Done.",SID_LOG_CLOSE);
  } // loop over snaps

  // Write the remaining snapshots
  for(;j_write>=i_read_start;i_write--,j_write-=i_read_step,l_write++)
     write_trees_horizontal(groups,   n_groups[l_write],   n_groups_max,   n_particles_groups[i_write%n_search],
                            subgroups,n_subgroups[l_write],n_subgroups_max,n_particles_subgroups[i_write%n_search],
                            n_subgroups_group,
                            n_halos_max,
                            max_tree_id_subgroup,
                            max_tree_id_group,
                            i_write,
                            j_write,
                            l_write,
                            n_search,
                            filename_cat_root_in,
                            filename_output_dir,
                            a_list,
                            cosmo,
                            n_k_match);
  
  SID_log("Freeing arrays...",SID_LOG_OPEN);
  SID_free(SID_FARG match_id);
  SID_free(SID_FARG match_score);
  SID_free(SID_FARG match_index);
  for(i_search=0;i_search<n_search;i_search++){
     SID_free(SID_FARG n_particles_groups[i_search]);
     SID_free(SID_FARG n_particles_subgroups[i_search]);
     SID_free(SID_FARG n_subgroups_group[i_search]);
     n_halos_max=n_subgroups_max;
     halos      =subgroups;
     for(i_halo=0;i_halo<n_halos_max;i_halo++){
        halos=subgroups;
        SID_free(SID_FARG halos[i_search][i_halo].bridges);
        SID_free(SID_FARG halos[i_search][i_halo].bridge_match);
        SID_free(SID_FARG halos[i_search][i_halo].bridge_backmatch);
     }
     n_halos_max=n_groups_max;
     halos      =groups;
     for(i_halo=0;i_halo<n_halos_max;i_halo++){
        halos=subgroups;
        SID_free(SID_FARG halos[i_search][i_halo].bridges);
        SID_free(SID_FARG halos[i_search][i_halo].bridge_match);
        SID_free(SID_FARG halos[i_search][i_halo].bridge_backmatch);
     }
     SID_free(SID_FARG subgroups[i_search]);
     SID_free(SID_FARG groups[i_search]);
  }
  SID_free(SID_FARG n_particles_groups);
  SID_free(SID_FARG n_particles_subgroups);
  SID_free(SID_FARG subgroups);
  SID_free(SID_FARG groups);
  SID_free(SID_FARG n_subgroups_group);
  SID_log("Done.",SID_LOG_CLOSE);

  // Set flag_clean=TRUE so that the output files generated here
  //  will be treated as temporary in the next step (if called)
  (*flag_clean)=TRUE;

  SID_log("Done.",SID_LOG_CLOSE);
}
