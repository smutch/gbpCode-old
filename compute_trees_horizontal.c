#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void change_horizontal_ID_recursive(tree_horizontal_info *halo,int id_1,int id_2);
void change_horizontal_ID_recursive(tree_horizontal_info *halo,int id_1,int id_2){
  tree_horizontal_info *current;

  // Change IDs here
  if(halo->id==id_1)
    halo->id=id_2;

  // Walk the tree
  current=halo->first_progenitor.halo;
  while(current!=NULL){
    change_horizontal_ID_recursive(current,id_1,id_2);
    current=current->next_progenitor.halo;
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
   stats->n_fragmented                =0;
   stats->n_emerged_progenitors       =0;
   stats->max_strayed_size            =0;
   stats->max_sputtered_size          =0;
   stats->max_dropped_size            =0;
   stats->max_bridged_size            =0;
   stats->max_bridge_progenitor_size  =0;
   stats->max_emerged_size            =0;
   stats->max_fragmented_size         =0;
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
         }
      }
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_EMERGED)){
         if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_FOUND)){
            stats->max_emerged_size=MAX(stats->max_emerged_size,halos[i_halo].n_particles);
            stats->n_emerged++;
         }
         else{
            stats->max_fragmented_size=MAX(stats->max_fragmented_size,halos[i_halo].n_particles);
            stats->n_fragmented++;
         }
      }

      // Count n_unprocessed and n_invalid so we can perform a couple sanity checks
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_INVALID)){
         n_invalid++;
         if((halos[i_halo].type-TREE_CASE_INVALID)!=0)
            SID_trap_error("An out-of-bounds halo has been manipulated (type=%d)",ERROR_LOGIC,halos[i_halo].type);
      }
      if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_UNPROCESSED))
         n_unprocessed++;
   }
   if(n_halos!=(n_halos_max-n_invalid))
      SID_trap_error("There is an incorrect number of out-of-bounds halos (i.e. %d!=%d)",ERROR_LOGIC,n_halos,(n_halos_max-n_invalid));
   if(n_unprocessed!=0)
      SID_trap_error("A number of halos (%d) have not been marked as processed",ERROR_LOGIC,n_unprocessed);
}

int match_id_local(match_info *match);
int match_id_local(match_info *match){
  if(match->halo!=NULL)
    return(match->halo->id);
  else
    return(-1);
}

int match_file_local(match_info *match);
int match_file_local(match_info *match){
  if(match->halo!=NULL)
    return(match->halo->file);
  else
    return(-1);
}

int match_type_local(match_info *match);
int match_type_local(match_info *match){
  if(match->halo!=NULL)
    return(match->halo->type);
  else
    return(-1);
}

int match_n_particles_local(match_info *match);
int match_n_particles_local(match_info *match){
  if(match->halo!=NULL)
    return(match->halo->n_particles);
  else
    return(-1);
}

int match_index_local(match_info *match);
int match_index_local(match_info *match){
  if(match->halo!=NULL)
    return(match->halo->index);
  else
    return(-1);
}

float match_score_local(match_info *match);
float match_score_local(match_info *match){
  if(match->halo!=NULL)
    return(match->score);
  else
    return(0.);
}

void write_trees_horizontal(tree_horizontal_info **groups,   int n_groups,    int n_groups_max,   
                            tree_horizontal_info **subgroups,int n_subgroups, int n_subgroups_max,
                            int   **n_subgroups_group,
                            int     n_halos_max,
                            int     max_tree_id_subgroup,
                            int     max_tree_id_group,
                            int     i_write,
                            int     j_write,
                            int     l_write,
                            int     i_file_start,
                            int     n_wrap,
                            char   *filename_cat_root_in,
                            char   *filename_output_dir,
                            double *a_list,
                            cosmo_info **cosmo,
                            int     n_k_match);
void write_trees_horizontal(tree_horizontal_info **groups,   int n_groups,    int n_groups_max,   
                            tree_horizontal_info **subgroups,int n_subgroups, int n_subgroups_max,
                            int   **n_subgroups_group,
                            int     n_halos_max,
                            int     max_tree_id_subgroup,
                            int     max_tree_id_group,
                            int     i_write,
                            int     j_write,
                            int     l_write,
                            int     n_wrap,
                            int     i_file_start,
                            char   *filename_cat_root_in,
                            char   *filename_output_dir,
                            double *a_list,
                            cosmo_info **cosmo,
                            int     n_k_match){
   char        filename_output_dir_horizontal[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_trees[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_groups[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_groups_stats[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_groups_matching[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_groups_properties[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_subgroups[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_subgroups_stats[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_subgroups_matching[MAX_FILENAME_LENGTH];
   char        filename_output_dir_horizontal_subgroups_properties[MAX_FILENAME_LENGTH];
   char        filename_output_file_root[MAX_FILENAME_LENGTH];
   char        filename_matches_out[MAX_FILENAME_LENGTH];
   char        filename_matching_out[MAX_FILENAME_LENGTH];
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
   char        filename_emerged_out[MAX_FILENAME_LENGTH];
   char        filename_fragmented_out[MAX_FILENAME_LENGTH];
   char       *filename_output_dir_stats;
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
   SID_fp      fp_group_properties_out;
   SID_fp      fp_subgroup_properties_out;
   FILE       *fp_group_properties_in;
   FILE       *fp_subgroup_properties_in;
   FILE       *fp;
   tree_horizontal_info  *halos;
   tree_horizontal_info **halos_all;
   int         n_halos;
   int         n_found;
   int         j_halo;
   int         k_halo;
   int         file_offset;
   int         i_subgroup;
   int         j_subgroup;
   int         i_group;
   int         j_group;
   int         i_k_match;
   int         j_k_match;
   int         i_column;
   int         n_particles_emerged;
   int         n_particles_emerged_main;
   int         n_particles_fragmented;
   int         n_p_largest_main;
   int         n_p_largest_i;
   int         n_p_largest_index;
   double      dt_descendant;
   double      dt_progenitor;
   char       *line=NULL;
   int         line_length=0;
   halo_info                  *properties;
   tree_horizontal_stats_info  stats;
   int                         desc_id;

   SID_log("Writing results for snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,j_write);
   sprintf(filename_output_dir_horizontal,                     "%s/horizontal",filename_output_dir);
   sprintf(filename_output_dir_horizontal_trees,               "%s/trees",     filename_output_dir_horizontal);
   sprintf(filename_output_dir_horizontal_groups,              "%s/groups",    filename_output_dir_horizontal);
   sprintf(filename_output_dir_horizontal_groups_stats,        "%s/stats",     filename_output_dir_horizontal_groups);
   sprintf(filename_output_dir_horizontal_groups_properties,   "%s/properties",filename_output_dir_horizontal_groups);
   sprintf(filename_output_dir_horizontal_groups_matching,     "%s/matching",  filename_output_dir_horizontal_groups);
   sprintf(filename_output_dir_horizontal_subgroups,           "%s/subgroups", filename_output_dir_horizontal);
   sprintf(filename_output_dir_horizontal_subgroups_stats,     "%s/stats",     filename_output_dir_horizontal_subgroups);
   sprintf(filename_output_dir_horizontal_subgroups_properties,"%s/properties",filename_output_dir_horizontal_subgroups);
   sprintf(filename_output_dir_horizontal_subgroups_matching,  "%s/matching",  filename_output_dir_horizontal_subgroups);
   mkdir(filename_output_dir,                                02755);
   mkdir(filename_output_dir_horizontal,                     02755);
   mkdir(filename_output_dir_horizontal_trees,               02755);
   mkdir(filename_output_dir_horizontal_groups,              02755);
   mkdir(filename_output_dir_horizontal_groups_stats,        02755);
   mkdir(filename_output_dir_horizontal_groups_matching,     02755);
   mkdir(filename_output_dir_horizontal_groups_properties,   02755);
   mkdir(filename_output_dir_horizontal_subgroups,           02755);
   mkdir(filename_output_dir_horizontal_subgroups_stats,     02755);
   mkdir(filename_output_dir_horizontal_subgroups_matching,  02755);
   mkdir(filename_output_dir_horizontal_subgroups_properties,02755);
   strcpy(filename_output_file_root,filename_output_dir);
   strip_path(filename_output_file_root);
   sprintf(filename_group_properties_out,   "%s/%s.trees_horizontal_groups_properties_%d",   filename_output_dir_horizontal_groups_properties,   
                                                                                             filename_output_file_root,j_write);
   sprintf(filename_subgroup_properties_out,"%s/%s.trees_horizontal_subgroups_properties_%d",filename_output_dir_horizontal_subgroups_properties,
                                                                                             filename_output_file_root,j_write);
   sprintf(filename_group_properties_in,    "%s_%03d.catalog_groups_properties",             filename_cat_root_in,j_write);
   sprintf(filename_subgroup_properties_in, "%s_%03d.catalog_subgroups_properties",          filename_cat_root_in,j_write);
   sprintf(filename_matches_out,            "%s/%s.trees_horizontal_%d",                     filename_output_dir_horizontal_trees,filename_output_file_root,j_write);

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
       if(i_write<i_file_start)
         file_offset=match_file_local(&(groups[i_write%n_wrap][i_group].descendant))-i_write;
       else
         file_offset=-1;
       desc_id    =match_id_local(&(groups[i_write%n_wrap][i_group].descendant));
       SID_fwrite(&(groups[i_write%n_wrap][i_group].id),        sizeof(int),1,&fp_matches_out);
       SID_fwrite(&(desc_id),                                   sizeof(int),1,&fp_matches_out);
       SID_fwrite(&(groups[i_write%n_wrap][i_group].tree_id),   sizeof(int),1,&fp_matches_out);
       SID_fwrite(&(file_offset),                               sizeof(int),1,&fp_matches_out);
       SID_fwrite(&(n_subgroups_group[i_write%n_wrap][i_group]),sizeof(int),1,&fp_matches_out);
       read_group_properties(fp_group_properties_in,properties,i_group,j_write);
       if(groups[i_write%n_wrap][i_group].id>=0)
         SID_fwrite(properties,sizeof(halo_info),1,&fp_group_properties_out);
       for(j_subgroup=0;j_subgroup<n_subgroups_group[i_write%n_wrap][i_group];j_subgroup++,i_subgroup++){
         // compute_trees_verticle wants the file offsets to be -ve for the roots, +ve everywhere else
         if(i_write<i_file_start)
           file_offset=match_file_local(&(subgroups[i_write%n_wrap][i_subgroup].descendant))-i_write;
         else
           file_offset=-1;
         desc_id    =match_id_local(&(subgroups[i_write%n_wrap][i_subgroup].descendant));
         SID_fwrite(&(subgroups[i_write%n_wrap][i_subgroup].id),     sizeof(int),1,&fp_matches_out);
         SID_fwrite(&(desc_id),                                      sizeof(int),1,&fp_matches_out);
         SID_fwrite(&(subgroups[i_write%n_wrap][i_subgroup].tree_id),sizeof(int),1,&fp_matches_out);
         SID_fwrite(&(file_offset),                                  sizeof(int),1,&fp_matches_out);
         read_group_properties(fp_subgroup_properties_in,properties,i_subgroup,j_write);
         if(subgroups[i_write%n_wrap][i_subgroup].id>=0)
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
      // Initialize a bunch of stuff depending on whether
      //   we are processing groups or subgroups
      switch(i_k_match){
         case 0:
         compute_trees_horizontal_stats(subgroups[i_write%n_wrap],n_subgroups,n_subgroups_max,&stats);
         filename_output_dir_stats=filename_output_dir_horizontal_subgroups_stats;
         halos    =subgroups[i_write%n_wrap];
         halos_all=subgroups;
         n_halos  =n_subgroups;
         break;
         case 1:
         compute_trees_horizontal_stats(groups[i_write%n_wrap],   n_groups,   n_groups_max,   &stats);
         filename_output_dir_stats=filename_output_dir_horizontal_groups_stats;
         halos    =groups[i_write%n_wrap];
         halos_all=groups;
         n_halos  =n_groups;
         break;
      }

      // Write snapshot summary statistics
      if(i_k_match==0 && l_write==0){
         fp=fopen(filename_log,"w");
         i_column=1;
         fprintf(fp,"# (%02d): Expansion factor (a)\n",   i_column++);
         fprintf(fp,"# (%02d): Snapshot filenumber\n",    i_column++);
         fprintf(fp,"# (%02d): Snapshot interval [yrs]\n",i_column++);
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
            fprintf(fp,"# (%02d): # of simple        %sgroups\n",           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of merging       %sgroups\n",           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of strayed       %sgroups\n",           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of sputtering    %sgroups\n",           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of dropped       %sgroups\n",           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of bridged       %sgroups\n",           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of bridged       %sgroup progenitors\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of emerged       %sgroups\n",           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of fragmented    %sgroups\n",           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): # of emerged       %sgroup progenitors\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest strayed    %sgroup\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest sputtered  %sgroup\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest dropped    %sgroup\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest bridged    %sgroup\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest emerged    %sgroup\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): Largest fragmented %sgroup\n",            i_column++,group_text_prefix);
         }
         fclose(fp);
      }
      fp=fopen(filename_log,"a");
      if(i_k_match==0){
         if(l_write>0)
           fprintf(fp,"%le %4d %10.4lf",a_list[l_write],j_write,deltat_a(cosmo,a_list[l_write],a_list[l_write-1])/S_PER_YEAR);
         else
           fprintf(fp,"%le %4d %10.4lf",a_list[l_write],j_write,deltat_a(cosmo,a_list[l_write+1],a_list[l_write])/S_PER_YEAR);
      }
      fprintf(fp," %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d",
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
              stats.n_fragmented,
              stats.n_emerged_progenitors,
              stats.max_strayed_size,
              stats.max_sputtered_size,
              stats.max_dropped_size,
              stats.max_bridged_size,
              stats.max_emerged_size,
              stats.max_fragmented_size);
      if(i_k_match==n_k_match-1)
         fprintf(fp,"\n");
      fclose(fp);

      switch(i_k_match){
         case 0:
         sprintf(group_text_prefix,"sub");
         break;
         case 1:
         sprintf(group_text_prefix,"");
         break;
      }

      // Write matching and special case information
      sprintf(filename_matching_out,          "%s/%s.progenitors_%sgroups_%d",     filename_output_dir_horizontal_subgroups_matching,
                                                                                   filename_output_file_root,group_text_prefix,j_write);
      sprintf(filename_mergers_out,           "%s/%s.%sgroups_mergers",            filename_output_dir_stats,filename_output_file_root,group_text_prefix);
      sprintf(filename_strayed_out,           "%s/%s.%sgroups_strays",             filename_output_dir_stats,filename_output_file_root,group_text_prefix);
      sprintf(filename_sputtered_out,         "%s/%s.%sgroups_sputters",           filename_output_dir_stats,filename_output_file_root,group_text_prefix);
      sprintf(filename_dropped_out,           "%s/%s.%sgroups_drops",              filename_output_dir_stats,filename_output_file_root,group_text_prefix);
      sprintf(filename_bridged_out,           "%s/%s.%sgroups_bridges",            filename_output_dir_stats,filename_output_file_root,group_text_prefix);      
      sprintf(filename_emerged_out,           "%s/%s.%sgroups_emerged",            filename_output_dir_stats,filename_output_file_root,group_text_prefix);      
      sprintf(filename_fragmented_out,        "%s/%s.%sgroups_fragmented",         filename_output_dir_stats,filename_output_file_root,group_text_prefix);      
      fp_matching_out=fopen(filename_matching_out,"w");
      fp             =fp_matching_out;
      i_column=1;
      fprintf(fp,"# (%02d): Progenitor number\n",                i_column++,group_text_prefix);
      fprintf(fp,"# (%02d): Tree ID\n",                          i_column++,group_text_prefix);
      fprintf(fp,"# (%02d): Halo file\n",                        i_column++,group_text_prefix);
      fprintf(fp,"# (%02d): Halo index\n",                       i_column++,group_text_prefix);
      fprintf(fp,"# (%02d): Halo ID\n",                          i_column++,group_text_prefix);
      fprintf(fp,"# (%02d): Progenitor file\n",                  i_column++,group_text_prefix);
      fprintf(fp,"# (%02d): Progenitor index\n",                 i_column++,group_text_prefix);
      fprintf(fp,"# (%02d): Progenitor ID\n",                    i_column++,group_text_prefix);
      fprintf(fp,"# (%02d): Descendant match type\n",            i_column++,group_text_prefix);
      fprintf(fp,"# (%02d): Progenitor match type\n",            i_column++,group_text_prefix);
      fprintf(fp,"# (%02d): Progenitor match score\n",           i_column++,group_text_prefix);
      fprintf(fp,"# (%02d): Number of particles in halo\n",      i_column++,group_text_prefix);
      fprintf(fp,"# (%02d): Number of particles in progenitor\n",i_column++,group_text_prefix);
      if(l_write==0){
         fp_mergers_out           =fopen(filename_mergers_out,           "w");
         fp_strayed_out           =fopen(filename_strayed_out,           "w");
         fp_sputtered_out         =fopen(filename_sputtered_out,         "w");
         fp_dropped_out           =fopen(filename_dropped_out,           "w");
         fp_bridged_out           =fopen(filename_bridged_out,           "w");
         fp_emerged_out     =fopen(filename_emerged_out,     "w");
         fp_fragmented_out   =fopen(filename_fragmented_out,   "w");
      }
      else{
         fp_mergers_out           =fopen(filename_mergers_out,           "a");
         fp_strayed_out           =fopen(filename_strayed_out,           "a");
         fp_sputtered_out         =fopen(filename_sputtered_out,         "a");
         fp_dropped_out           =fopen(filename_dropped_out,           "a");
         fp_bridged_out           =fopen(filename_bridged_out,           "a");
         fp_emerged_out     =fopen(filename_emerged_out,     "a");
         fp_fragmented_out   =fopen(filename_fragmented_out,   "a");
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
            if(halos[i_halo].first_progenitor.halo->file<i_write){ // This isn't true for first snapshot for instance
              dt_progenitor=deltat_a(cosmo,a_list[l_write+(i_write-halos[i_halo].first_progenitor.halo->file)],a_list[l_write])/S_PER_YEAR;
            }
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
            fprintf(fp_matching_out,"%3d %7d   %4d %7d %7d   %4d %7d %7d   %8d %8d   %10.4f   %7d %7d\n",
                    j_halo++,
                    halos[i_halo].tree_id,
                    j_write,
                    i_halo,
                    halos[i_halo].id,
                    match_file_local(current_match),
                    match_index_local(current_match),
                    match_id_local(current_match),
                    halos[i_halo].type,
                    match_type_local(current_match),
                    match_score_local(current_match),
                    halos[i_halo].n_particles,
                    match_n_particles_local(current_match));
            current_match=&(current_match->halo->next_progenitor);
         }

         // Write mergers
         fp=fp_mergers_out;
         if(l_write==0 && i_halo==0){
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                                   i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                                    i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                                     i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup id\n",                                                 i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup match type\n",                                         i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",                         i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n",             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup merger's main progenitor\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup's matching score to it's descendant\n",                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup's matching score to it's main progenitor\n",           i_column++,group_text_prefix);
         }
         if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_MERGER))
            fprintf(fp,"%10.3le %8d %8d %8d %8d %8d %8d %8d %8d %8d %10.3f %10.3f\n",
                    a_list[l_write],
                    j_write,
                    i_halo,
                    halos[i_halo].id,
                    halos[i_halo].type,
                    match_file_local(&(halos[i_halo].descendant))-i_write,
                    halos[i_halo].n_particles,
                    halos[i_halo].n_particles_parent,
                    match_n_particles_local(&(halos[i_halo].descendant)),
                    match_n_particles_local(&(halos[i_halo].descendant.halo->first_progenitor)),
                    match_score_local(&(halos[i_halo].descendant)),
                    match_score_local(&(halos[i_halo].first_progenitor)));

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
         if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_STRAYED))
            fprintf(fp,"%10.3le %8d %8d %8d %8d %8d\n",
               a_list[l_write],
               j_write,
               i_halo,
               halos[i_halo].type,
               halos[i_halo].n_particles,
               halos[i_halo].n_particles_parent);

         // Write sputtered halos
         fp=fp_sputtered_out;
         if(l_write==0 && i_halo==0){
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                       i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                        i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                         i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup match type\n",                             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                 i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n", i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",i_column++,group_text_prefix);
         }
         if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_SPUTTERED))
            fprintf(fp,"%10.3le %8d %8d %8d %8d %8d %8d %8d\n",
               a_list[l_write],
               j_write,
               i_halo,
               halos[i_halo].type,
               match_file_local(&(halos[i_halo].descendant))-i_write,
               halos[i_halo].n_particles,
               halos[i_halo].n_particles_parent,
               match_n_particles_local(&(halos[i_halo].descendant)));

         // Write dropped halos
         fp=fp_dropped_out;
         if(l_write==0 && i_halo==0){
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup->descendant interval [yrs]\n",                  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                              i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup id\n",                                          i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup match type\n",                                  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                      i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",                  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n",      i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",     i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's main progenitor\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's descendant\n",            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's main progenitor\n",       i_column++,group_text_prefix);
         }
         if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_DROPPED))
            fprintf(fp,"%10.3le %10.3le %8d %8d %8d %8d %8d %8d %8d %8d %8d %10.3f %10.3f\n",
               a_list[l_write],
               dt_descendant,
               j_write,
               i_halo,
               halos[i_halo].id,
               halos[i_halo].type,
               match_file_local(&(halos[i_halo].descendant))-i_write,
               halos[i_halo].n_particles,
               halos[i_halo].n_particles_parent,
               match_n_particles_local(&(halos[i_halo].descendant)),
               match_n_particles_local(&(halos[i_halo].first_progenitor)),
               match_score_local(&(halos[i_halo].descendant)),
               match_score_local(&(halos[i_halo].first_progenitor)));

         // Write bridged halos
         fp=fp_bridged_out;
         if(l_write==0 && i_halo==0){
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                                                      i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                                                       i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                                                        i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup id\n",                                                                    i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup match type\n",                                                            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of back-matched %sgroups identified with this bridge\n",                   i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of back-matched %sgroups identified as emerging halos\n",                  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",                                            i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n",                                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's main progenitor\n",                          i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the largest    back-matched %sgroup\n",                    i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the largest    back-matched %sgroup's main progenitor\n",  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the emerged    back-matched %sgroups\n",                   i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the emerged    back-matched %sgroups's main progenitors\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the fragmented back-matched %sgroups\n",                   i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the fragmented back-matched %sgroups's main progenitors\n",i_column++,group_text_prefix);
         }
         if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_BRIDGED)){
            // Find the largest back-matched halo
            j_halo=0;
            n_p_largest_i=match_n_particles_local(&(halos_all[match_file_local(&(halos[i_halo].bridges[j_halo]))%n_wrap]
                                                             [match_index_local(&(halos[i_halo].bridges[j_halo]))].first_progenitor));
            n_p_largest_main =n_p_largest_i;
            n_p_largest_index=j_halo;
            for(j_halo=1,n_p_largest_main=0,n_p_largest_index=0;j_halo<halos[i_halo].n_bridges;j_halo++){
               n_p_largest_i=match_n_particles_local(&(halos_all[match_file_local(&(halos[i_halo].bridges[j_halo]))%n_wrap]
                                                                [match_index_local(&(halos[i_halo].bridges[j_halo]))].first_progenitor));
               if(n_p_largest_i>n_p_largest_main){
                  n_p_largest_main =n_p_largest_i;
                  n_p_largest_index=j_halo;
               }
            }

            // Count the number of emerged/fragmented halos and the number of particles
            //   in each set, as well as the number of particles in their main progenitors
            n_found=0;
            n_particles_emerged        =0;
            n_particles_emerged_main   =0;
            n_particles_fragmented     =0;
            for(j_halo=0;j_halo<halos[i_halo].n_bridges;j_halo++){
               if(check_mode_for_flag(match_type_local(&(halos[i_halo].bridges[j_halo])),TREE_CASE_FOUND)){
                  n_particles_emerged     +=match_n_particles_local(&(halos[i_halo].bridges[j_halo]));
                  n_particles_emerged_main+=match_n_particles_local(&(halos_all[match_file_local(&(halos[i_halo].bridges[j_halo]))%n_wrap]
                                                                             [match_index_local(&(halos[i_halo].bridges[j_halo]))].first_progenitor));
                  n_found++;
               }
               else
                  n_particles_fragmented+=match_n_particles_local(&(halos[i_halo].bridges[j_halo]));
            }
            fprintf(fp,"%10.3le %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d\n",
                    a_list[l_write],
                    j_write,
                    i_halo,
                    halos[i_halo].id,
                    halos[i_halo].type,
                    halos[i_halo].n_bridges,
                    n_found,
                    halos[i_halo].n_particles,
                    halos[i_halo].n_particles_parent,
                    match_n_particles_local(&(halos[i_halo].first_progenitor)),
                    match_n_particles_local(&(halos[i_halo].bridges[n_p_largest_index])),
                    n_p_largest_main,
                    n_particles_emerged,
                    n_particles_emerged_main,
                    n_particles_fragmented);
         }

         // Write emerged and fragmented halos
         fp=fp_emerged_out;
         if(l_write==0 && i_halo==0){
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup->progenitor interval [yrs]\n",                      i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                                 i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                                  i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup id\n",                                              i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup match type\n",                                      i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                          i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup bridge snapshot\n",                                 i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup bridge index\n",                                    i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup bridge id\n",                                       i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's bridge\n",                    i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's descendant\n",                i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's main progenitor\n",           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",                      i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",         i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n",          i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's bridge\n",             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's main progenitor\n",    i_column++,group_text_prefix);
         }
         fp=fp_fragmented_out;
         if(l_write==0 && i_halo==0){
            i_column=1;
            fprintf(fp,"# (%02d): %sgroup expansion factor\n",                       i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot number\n",                        i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup snapshot index\n",                         i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup id\n",                                     i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup match type\n",                             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup descendant file offset\n",                 i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup bridge snapshot\n",                        i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup bridge index\n",                           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): %sgroup bridge id\n",                              i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's bridge\n",           i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): match score for the %sgroup's descendant\n",       i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup\n",             i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's descendant\n",i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's FoF group\n", i_column++,group_text_prefix);
            fprintf(fp,"# (%02d): number of particles in the %sgroup's bridge\n",    i_column++,group_text_prefix);
         }
         if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_EMERGED)){
            if(check_mode_for_flag(halos[i_halo].type,TREE_CASE_FOUND)){
               fp=fp_emerged_out;
               fprintf(fp,"%10.3le %10.3le %8d %8d %8d %8d %8d %8d %8d %8d %10.3f %10.3f %10.3f %8d %8d %8d %8d %8d\n",
                       a_list[l_write],
                       dt_progenitor,
                       j_write,
                       i_halo,
                       halos[i_halo].id,
                       halos[i_halo].type,
                       match_file_local(&(halos[i_halo].descendant))-i_write,
                       match_file_local(&(halos[i_halo].bridge_backmatch)),
                       match_index_local(&(halos[i_halo].bridge_backmatch)),
                       match_id_local(&(halos[i_halo].bridge_backmatch)),
                       match_score_local(&(halos[i_halo].bridge_backmatch)),
                       match_score_local(&(halos[i_halo].descendant)),
                       match_score_local(&(halos[i_halo].first_progenitor)),
                       halos[i_halo].n_particles,
                       match_n_particles_local(&(halos[i_halo].descendant)),
                       halos[i_halo].n_particles_parent,
                       match_n_particles_local(&(halos[i_halo].bridge_backmatch)),
                       match_n_particles_local(&(halos[i_halo].first_progenitor)));
            }
            else{
               fp=fp_fragmented_out;
               fprintf(fp,"%10.3le %8d %8d %8d %8d %8d %8d %8d %8d %10.3f %10.3f %8d %8d %8d %8d\n",
                       a_list[l_write],
                       j_write,
                       i_halo,
                       halos[i_halo].id,
                       halos[i_halo].type,
                       match_file_local(&(halos[i_halo].descendant))-i_write,
                       match_file_local(&(halos[i_halo].bridge_backmatch)),
                       match_index_local(&(halos[i_halo].bridge_backmatch)),
                       match_id_local(&(halos[i_halo].bridge_backmatch)),
                       match_score_local(&(halos[i_halo].bridge_backmatch)),
                       match_score_local(&(halos[i_halo].descendant)),
                       halos[i_halo].n_particles,
                       match_n_particles_local(&(halos[i_halo].descendant)),
                       halos[i_halo].n_particles_parent,
                       match_n_particles_local(&(halos[i_halo].bridge_backmatch)));
            }
         }
      }
      fclose(fp_matching_out);
      fclose(fp_mergers_out);
      fclose(fp_strayed_out);
      fclose(fp_sputtered_out);
      fclose(fp_dropped_out);
      fclose(fp_bridged_out);
      fclose(fp_emerged_out);
      fclose(fp_fragmented_out);
   }

   SID_log("Done.",SID_LOG_CLOSE);

}

// Set halo_i[i_halo] so it points to halo_j[j_halo]
void set_halo_and_descendant(tree_horizontal_info **halos,
                             int                    i_file,
                             int                    i_halo,
                             int                    j_file,
                             int                    j_halo,
                             float                  score,
                             int                   *max_id,
                             int                    n_wrap);
void set_halo_and_descendant(tree_horizontal_info **halos,
                             int                    i_file,
                             int                    i_halo,
                             int                    j_file,
                             int                    j_halo,
                             float                  score,
                             int                   *max_id,
                             int                    n_wrap){
   tree_horizontal_info *halos_i;
   tree_horizontal_info *halos_j;
   int                   file_offset;
   int                   k_file;
   int                   l_file;
   int                   k_index;
   int                   k_file_temp;
   int                   k_index_temp;
   int                   k_file_main;
   int                   k_index_main;
   int                   k_size_main;
   int                   flag_process;
   int                   n_p_diff_old;
   int                   n_p_diff_new;

   if(j_file<=i_file)
     SID_trap_error("j_file<=i_file in set_halo_and_descendant().",ERROR_NONE);

   // Process the inputs a bit
   halos_i    =halos[i_file%n_wrap];
   halos_j    =halos[j_file%n_wrap];
   file_offset=j_file-i_file;
   if(file_offset==0)
      SID_trap_error("A zero file offset has been requested.  It should be -ve for roots and +ve otherwise.",ERROR_LOGIC);
      
   // Set non-bridged halos or finalize bridge matches (ie. set defaults for bridge progenitors not matched to emerged halos)
   if(!check_mode_for_flag(halos_j[j_halo].type,TREE_CASE_BRIDGED)                       || 
       check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED) ||
       check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_FINALIZE)){

      // If we are processing a bridge progenitor, only accept
      //   this new match if it meets these criteria ...
      flag_process=TRUE;
      if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED) && 
         !check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_FINALIZE)){
        // ... the score must be good ...
        if(score<0.5*(halos_i[i_halo].bridge_forematch.score))
          flag_process=FALSE;
        // ... the change in halo size must be better ...
        n_p_diff_old=IABS(halos_i[i_halo].bridge_forematch.halo->n_particles-halos_i[i_halo].n_particles);
        n_p_diff_new=IABS(halos_j[j_halo].n_particles-halos_i[i_halo].n_particles);
        if(n_p_diff_new>=n_p_diff_old)
          flag_process=FALSE;
        // ... we must not be matching to a descendant of the initial bridged match ...
        tree_horizontal_info *current;
        current=halos_i[i_halo].bridge_forematch.halo->descendant.halo;
        if(current!=NULL)
           k_file =current->file;
        l_file=k_file;
        while(current!=NULL && k_file>=l_file && k_file<=j_file && flag_process){
           if(current==&(halos_j[j_halo]))
             flag_process=FALSE;
           current=current->descendant.halo;
           l_file=k_file;
           if(current!=NULL)
              k_file =current->file;
        }
      }
      if(flag_process){

         // Set progenitor IDs and pointers ...
         match_info old_progenitor;
         match_info new_progenitor;
         // ... create new progenitor ...
         new_progenitor.halo =&(halos_i[i_halo]);
         new_progenitor.score=score;
         // ... increment counter ...
         halos_j[j_halo].n_progenitors++;
         // ... set tree id ...
         new_progenitor.halo->tree_id=halos_j[j_halo].tree_id;
         // ... create initial progenitor ....
         if(halos_j[j_halo].n_progenitors==1){
            // ... set initial progenitor ...
            memcpy(&(halos_j[j_halo].first_progenitor),&new_progenitor,sizeof(match_info));
            // ... set initial progenitor id ...
            new_progenitor.halo->id=halos_j[j_halo].id;
            // ... set initial progenitor type ...
            new_progenitor.halo->type|=TREE_CASE_MAIN_PROGENITOR;
         }
         else{
            // If we have a higher-score (ie a new main) progenitor, insert it at the 
            //   beginning of the list and swap IDs with the main progenitor so that
            //   the correct halo gets the main progenitor ID and all others get a new one ...
            memcpy(&old_progenitor,&(halos_j[j_halo].first_progenitor),sizeof(match_info));
            if(score>old_progenitor.score){
               // ... set new main progenitor ...
               memcpy(&(halos_j[j_halo].first_progenitor),                      &new_progenitor,sizeof(match_info));
               memcpy(&(halos_j[j_halo].first_progenitor.halo->next_progenitor),&old_progenitor,sizeof(match_info));
               // ... let the new main progenitor inherit the descendant's id ...
               halos_i[i_halo].id=halos_j[j_halo].id;
               // ... and give the old main progenitor the new id ...
               change_horizontal_ID_recursive(old_progenitor.halo,old_progenitor.halo->id,(*max_id)++);
               // ... set new main progenitor type ...
               old_progenitor.halo->type|=  TREE_CASE_MERGER;
               old_progenitor.halo->type&=(~TREE_CASE_MAIN_PROGENITOR);
               new_progenitor.halo->type&=(~TREE_CASE_MERGER);
               new_progenitor.halo->type|=  TREE_CASE_MAIN_PROGENITOR;
            }
            // ... else just add the new halo to the end of the list and create a new ID for it (if this is not a strayed/sputtered halo).
            else{
               // ... set new non-main progenitor ...
               memcpy(&(halos_j[j_halo].last_progenitor.halo->next_progenitor),&new_progenitor,sizeof(match_info));
               // ... set new non-main progenitor id ...
               if(halos_j[j_halo].id>=0)
                  halos_i[i_halo].id=(*max_id)++;
               else
                  halos_i[i_halo].id=halos_j[j_halo].id;
               // ... set new non-main progenitor type ...
               halos_i[i_halo].type|=  TREE_CASE_MERGER;
               halos_i[i_halo].type&=(~TREE_CASE_MAIN_PROGENITOR);
            }
         }
         // ... set last progenitor
         memcpy(&(halos_j[j_halo].last_progenitor),&new_progenitor,sizeof(match_info));

         // Set descendant info
         halos_i[i_halo].descendant.halo =&(halos_j[j_halo]);
         halos_i[i_halo].descendant.score=score;

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
            if(!check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_PROGENITOR))
               halos_i[i_halo].type|=TREE_CASE_DROPPED;
            halos_j[j_halo].type|=TREE_CASE_FOUND;
         }

         // Mark emerged halos as found
         if(check_mode_for_flag(halos_j[j_halo].type,TREE_CASE_EMERGED))
            halos_j[j_halo].type|=TREE_CASE_FOUND;

         // Mark the halo as processed
         halos_i[i_halo].type&=(~(TREE_CASE_UNPROCESSED|TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED|TREE_CASE_BRIDGE_FINALIZE));
      }
   }
   // ... else we've matched to a bridge.  Here we just set the info needed
   //     to connect this halo with the bridge it's been matched to.  We'll
   //     use that info to search the bridge's emergent halos and if we fail
   //     to identify to one of those, we'll finalize this match as the default.
   else if(halos_i[i_halo].bridge_forematch.halo==NULL){
      halos_i[i_halo].bridge_forematch.halo =&(halos_j[j_halo]);
      halos_i[i_halo].bridge_forematch.score =score;
      halos_i[i_halo].type                  |=(TREE_CASE_BRIDGE_PROGENITOR|TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED);
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
                              int     flag_fix_bridges,
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
  int        *n_particles;
  int        *n_particles_groups;
  int        *n_particles_subgroups;
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
  bridge_info           *bridges;
  bridge_info           *bridge;

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
  int     l_file;
  int     k_index;
  int     k_file_temp;
  int     k_index_temp;

  int     n_wrap;
  
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
  int     i_file_start;
  
  tree_horizontal_stats_info stats;
 
  char  filename_output_dir_horizontal[MAX_FILENAME_LENGTH];
  char  filename_output_dir_horizontal_groups[MAX_FILENAME_LENGTH];
  char  filename_output_dir_horizontal_groups_matching[MAX_FILENAME_LENGTH];
  char  filename_output_dir_horizontal_subgroups[MAX_FILENAME_LENGTH];
  char  filename_output_dir_horizontal_subgroups_matching[MAX_FILENAME_LENGTH];
  char  filename_output_file_root[MAX_FILENAME_LENGTH];
  char  filename_matching_out[MAX_FILENAME_LENGTH];
  FILE *fp_matching_out;
  int   i_column;

 
  SID_log("Constructing horizontal merger trees for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_stop,i_read_start,i_read_step);

  if(!flag_fix_bridges)
    SID_log("Bridge-fixing is turned off.",SID_LOG_COMMENT);

  if(n_search<1)
    SID_trap_error("n_search=%d but must be at least 1",ERROR_LOGIC,n_search);

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

  // We need indices that allow us to hold-on to descendants until outputting
  //   and for the current and last i_file as well
  n_wrap=2*n_search+2;
     
  // Initialize arrays
  SID_log("Creating arrays...",SID_LOG_OPEN);
  n_particles_groups     =(int    *)SID_malloc(sizeof(int)   *n_halos_max);
  n_particles_subgroups  =(int    *)SID_malloc(sizeof(int)   *n_halos_max);
  match_id               =(int    *)SID_malloc(sizeof(int)   *n_halos_max);
  match_score            =(float  *)SID_malloc(sizeof(float) *n_halos_max);
  match_index            =(size_t *)SID_malloc(sizeof(size_t)*n_halos_max);
  subgroups              =(tree_horizontal_info **)SID_malloc(sizeof(tree_horizontal_info *)*n_wrap);
  groups                 =(tree_horizontal_info **)SID_malloc(sizeof(tree_horizontal_info *)*n_wrap);
  n_subgroups_group      =(int                  **)SID_malloc(sizeof(int                  *)*n_wrap);
  for(i_search=0;i_search<n_wrap;i_search++){
     subgroups[i_search]            =(tree_horizontal_info *)SID_calloc(sizeof(tree_horizontal_info)*n_subgroups_max);
     groups[i_search]               =(tree_horizontal_info *)SID_calloc(sizeof(tree_horizontal_info)*n_groups_max);       
     n_subgroups_group[i_search]    =(int                  *)SID_malloc(sizeof(int)                 *n_groups_max);       
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Process the first file separately
  //   (just give everything ids from a running index) ...
  SID_log("Initializing tree roots...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Initialize everything to a 1:1 simple match
  read_matches(filename_root_matches,
               i_read_stop,i_read_stop-i_read_step,
               MATCH_GROUPS,
               &n_halos_1_matches,
               &n_halos_2_matches,
               n_particles_groups,
               NULL,
               n_subgroups_group[0],
               n_subgroups_group[1],
               match_id,
               match_score,
               match_index);
  for(i_search=1;i_search<n_wrap;i_search++)
    memcpy(n_subgroups_group[i_search],n_subgroups_group[0],n_halos_1_matches*sizeof(int));

  i_file_start=n_files-1;

  for(i_halo=0,max_id_group=0,max_tree_id_group=0;i_halo<n_groups_max;i_halo++,j_halo++){
     for(i_search=0;i_search<n_wrap;i_search++){
        groups[i_search][i_halo].file                  =   i_file_start; // The resulting file offset must be -ve for tree roots
        groups[i_search][i_halo].index                 =(size_t)i_halo;
        groups[i_search][i_halo].n_bridges             =   0;
        groups[i_search][i_halo].descendant.halo       =NULL;
        groups[i_search][i_halo].descendant.score      =  0.;
        groups[i_search][i_halo].first_progenitor.halo =NULL;
        groups[i_search][i_halo].first_progenitor.score=  0.;
        groups[i_search][i_halo].last_progenitor.halo  =NULL;
        groups[i_search][i_halo].last_progenitor.score =  0.;
        groups[i_search][i_halo].next_progenitor.halo  =NULL;
        groups[i_search][i_halo].next_progenitor.score =  0.;
        groups[i_search][i_halo].bridge_forematch.halo =NULL;
        groups[i_search][i_halo].bridge_forematch.score=  0.;
        groups[i_search][i_halo].bridge_backmatch.halo =NULL;
        groups[i_search][i_halo].bridge_backmatch.score=  0.;
        groups[i_search][i_halo].bridges               =NULL;
        groups[i_search][i_halo].type                  =TREE_CASE_INVALID;
        groups[i_search][i_halo].id                    =-1;
        groups[i_search][i_halo].tree_id               =-1; 
        groups[i_search][i_halo].n_particles           = 0;
        groups[i_search][i_halo].n_particles_parent    = 0;
        groups[i_search][i_halo].n_progenitors         = 0;
        if(i_halo<n_halos_1_matches){
           groups[i_search][i_halo].id                    =max_id_group;
           groups[i_search][i_halo].tree_id               =max_tree_id_group;
           groups[i_search][i_halo].type                  =TREE_CASE_SIMPLE|TREE_CASE_MAIN_PROGENITOR;
           groups[i_search][i_halo].n_particles           =n_particles_groups[i_halo];
           groups[i_search][i_halo].n_particles_parent    =n_particles_groups[i_halo];
           groups[i_search][i_halo].descendant.halo       =&(subgroups[(i_search+1)%n_wrap][i_halo]);
           groups[i_search][i_halo].descendant.score       =1.;
           if(i_search!=(i_file_start%n_wrap)){
              groups[i_search][i_halo].n_progenitors       =1;
              if(i_search>0){
                 groups[i_search][i_halo].first_progenitor.halo =&(groups[i_search-1][i_halo]);
                 groups[i_search][i_halo].last_progenitor.halo  =&(groups[i_search-1][i_halo]);
              }
              else{
                 groups[i_search][i_halo].first_progenitor.halo =&(groups[n_wrap-1][i_halo]);
                 groups[i_search][i_halo].last_progenitor.halo  =&(groups[n_wrap-1][i_halo]);
              }
              groups[i_search][i_halo].first_progenitor.score=1.;
              groups[i_search][i_halo].last_progenitor.score =1.;
           }
        }
     }
     if(i_halo<n_halos_1_matches){
        max_id_group++;
        max_tree_id_group++;
     }
  }

  // Initialize everything to a 1:1 simple match
  read_matches(filename_root_matches,
               i_read_stop,i_read_stop-i_read_step,
               MATCH_SUBGROUPS,
               &n_halos_1_matches,
               &n_halos_2_matches,
               n_particles_subgroups,
               NULL,
               NULL,
               NULL,
               match_id,
               match_score,
               match_index);
                     
  for(i_halo=0,j_halo=0,k_halo=0,max_id_subgroup=0,max_tree_id_subgroup=0;i_halo<n_subgroups_max;i_halo++,j_halo++){
     if(j_halo>n_subgroups_group[0][k_halo] && i_halo<n_halos_1_matches){
       k_halo++;
       j_halo=0;
     }
     for(i_search=0;i_search<n_wrap;i_search++){
        subgroups[i_search][i_halo].file                  =   i_file_start; // The resulting file offset must be -ve for tree roots
        subgroups[i_search][i_halo].index                 =(size_t)i_halo;
        subgroups[i_search][i_halo].n_bridges             =   0;
        subgroups[i_search][i_halo].descendant.halo       =NULL;
        subgroups[i_search][i_halo].descendant.score      =  0.;
        subgroups[i_search][i_halo].first_progenitor.halo =NULL;
        subgroups[i_search][i_halo].first_progenitor.score=  0.;
        subgroups[i_search][i_halo].last_progenitor.halo  =NULL;
        subgroups[i_search][i_halo].last_progenitor.score =  0.;
        subgroups[i_search][i_halo].next_progenitor.halo  =NULL;
        subgroups[i_search][i_halo].next_progenitor.score =  0.;
        subgroups[i_search][i_halo].bridge_forematch.halo =NULL;
        subgroups[i_search][i_halo].bridge_forematch.score=  0.;
        subgroups[i_search][i_halo].bridge_backmatch.halo =NULL;
        subgroups[i_search][i_halo].bridge_backmatch.score=  0.;
        subgroups[i_search][i_halo].bridges               =NULL;
        subgroups[i_search][i_halo].type                  =TREE_CASE_INVALID;
        subgroups[i_search][i_halo].id                    =-1;
        subgroups[i_search][i_halo].tree_id               =-1; 
        subgroups[i_search][i_halo].n_particles           = 0;
        subgroups[i_search][i_halo].n_particles_parent    = 0;
        subgroups[i_search][i_halo].n_progenitors         = 0;
        if(i_halo<n_halos_1_matches){
           subgroups[i_search][i_halo].id                    =max_id_subgroup;
           subgroups[i_search][i_halo].tree_id               =max_tree_id_subgroup;
           subgroups[i_search][i_halo].type                  =TREE_CASE_SIMPLE|TREE_CASE_MAIN_PROGENITOR;
           subgroups[i_search][i_halo].n_particles           =n_particles_subgroups[i_halo];
           subgroups[i_search][i_halo].n_particles_parent    =n_particles_groups[k_halo];
           subgroups[i_search][i_halo].descendant.halo       =&(subgroups[(i_search+1)%n_wrap][i_halo]);
           subgroups[i_search][i_halo].descendant.score       =1.;
           if(i_search!=(i_file_start%n_wrap)){
              subgroups[i_search][i_halo].n_progenitors       =1;
              if(i_search>0){
                 subgroups[i_search][i_halo].first_progenitor.halo =&(subgroups[i_search-1][i_halo]);
                 subgroups[i_search][i_halo].last_progenitor.halo  =&(subgroups[i_search-1][i_halo]);
              }
              else{
                 subgroups[i_search][i_halo].first_progenitor.halo =&(subgroups[n_wrap-1][i_halo]);
                 subgroups[i_search][i_halo].last_progenitor.halo  =&(subgroups[n_wrap-1][i_halo]);
              }
              subgroups[i_search][i_halo].first_progenitor.score=1.;
              subgroups[i_search][i_halo].last_progenitor.score =1.;
           }
        }
     }
     if(i_halo<n_halos_1_matches){
        max_id_subgroup++;
        max_tree_id_subgroup++;
     }
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
        i_file =i_file_start-1, 
        j_file =1,             
        i_write=i_file_start,      
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
       halos_i  =halos[i_file%n_wrap];
       n_halos_i=n_halos[j_file];
       SID_log("Processing %d %sgroups...",SID_LOG_OPEN|SID_LOG_TIMER,n_halos_i,group_text_prefix);

       // Initialize tree pointer-arrays with unique and negative dummy values
       for(i_halo=0;i_halo<n_halos_max;i_halo++){
          halos_i[i_halo].file                  = i_file;
          halos_i[i_halo].index                 = (size_t)i_halo;
          halos_i[i_halo].n_bridges             =   0;
          halos_i[i_halo].descendant.halo       =NULL;
          halos_i[i_halo].descendant.score      =  0.;
          halos_i[i_halo].first_progenitor.halo =NULL;
          halos_i[i_halo].first_progenitor.score=  0.;
          halos_i[i_halo].last_progenitor.halo  =NULL;
          halos_i[i_halo].last_progenitor.score =  0.;
          halos_i[i_halo].next_progenitor.halo  =NULL;
          halos_i[i_halo].next_progenitor.score =  0.;
          halos_i[i_halo].bridge_forematch.halo =NULL;
          halos_i[i_halo].bridge_forematch.score=  0.;
          halos_i[i_halo].bridge_backmatch.halo =NULL;
          halos_i[i_halo].bridge_backmatch.score=  0.;
          SID_free(SID_FARG halos_i[i_halo].bridges);
          if(i_halo<n_halos_i)
             halos_i[i_halo].type=TREE_CASE_UNPROCESSED;
          else
             halos_i[i_halo].type=TREE_CASE_INVALID;
          halos_i[i_halo].id              =-1;
          halos_i[i_halo].tree_id         =-1;
          halos_i[i_halo].n_particles     = 0;
          halos_i[i_halo].n_progenitors   = 0;
       }
       
       // Use back-matching to identify bridged halos ...
       if(flag_fix_bridges){
          SID_log("Identifying bridge candidates from back-matching...",SID_LOG_OPEN|SID_LOG_TIMER);
          SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
          //    ... first, do an initial count of matches.  This will not be a list of unique halos
          //        though, since the same halos are likely to appear in repeated snapshots.
          for(j_file_1  =i_file+1,
                j_file_2=i_file,
                j_read_1=i_read+i_read_step,
                j_read_2=i_read,
                i_search=0;
              j_read_1<=i_read_stop && i_search<n_search;
              j_file_1++,
                j_read_1+=i_read_step,
                i_search++){

             SID_log("Counting matches between files %d->%d...",SID_LOG_OPEN,j_read_1,j_read_2);

             // Read back-matching
             read_matches(filename_root_matches,
                          j_read_1,j_read_2,
                          flag_match_subgroups,
                          &n_halos_1_matches,
                          &n_halos_2_matches,
                          NULL,
                          n_particles,
                          NULL,
                          NULL,
                          match_id,
                          match_score,
                          match_index);

             // Store halo sizes
             if(i_search==0){
                for(i_halo=0;i_halo<n_halos_2_matches;i_halo++)
                   halos[i_file%n_wrap][i_halo].n_particles=n_particles[i_halo];
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
                (halos_i[i_halo].bridges)=(bridge_info *)SID_calloc(sizeof(bridge_info)*(halos_i[i_halo].n_bridges));
             else
                (halos_i[i_halo].bridges)=NULL;
             halos_i[i_halo].n_bridges =0;
          }

          //    ... third, assemble the list of unique back-matched halos.
          for(j_file_1  =i_file+1,
                j_file_2=i_file,
                j_read_1=i_read+i_read_step,
                j_read_2=i_read,
                i_search=0;
              j_read_1<=i_read_stop && i_search<n_search;
              j_file_1++,
                j_read_1+=i_read_step,
                i_search++){

             SID_log("Finding unique matches between files %d->%d...",SID_LOG_OPEN,j_read_1,j_read_2);
          
             // Read back-matching
             read_matches(filename_root_matches,
                          j_read_1,j_read_2,
                          flag_match_subgroups,
                          &n_halos_1_matches,
                          &n_halos_2_matches,
                          NULL,
                          NULL,
                          NULL,
                          NULL,
                          match_id,
                          match_score,
                          match_index);
                             
             // For all the halos in i_file_1 with back-matches ...
             for(i_halo=0;i_halo<n_halos_i;i_halo++){
                if((halos_i[i_halo].bridges)!=NULL){
                   // Scan over the list of halos from snapshot=j_read_1 
                   //   that match this halo in j_read_2 ...
                   bridges=halos_i[i_halo].bridges;
                   j_halo =find_index_int(match_id,i_halo,n_halos_1_matches,match_index);
                   // Loop over all but the last halo in the list ...
                   while(match_id[match_index[j_halo]]==i_halo && j_halo<(n_halos_1_matches-1)){
                      // Check to see if this halo is already in the bridge list ...
                      for(k_halo=0,flag_continue=TRUE;k_halo<halos_i[i_halo].n_bridges && flag_continue;k_halo++){
                         if(bridges[k_halo].halo->id==halos[j_file_1%n_wrap][match_index[j_halo]].id)
                            flag_continue=FALSE;
                      }
                      // ... if not, add it
                      if(flag_continue){
                         bridges[halos_i[i_halo].n_bridges].score=match_score[match_index[j_halo]];
                         bridges[halos_i[i_halo].n_bridges].halo =&(halos[j_file_1%n_wrap][match_index[j_halo]]);
                         (halos_i[i_halo].n_bridges)++;
                      }
                      j_halo++;
                   }
                   // ... then do the last halo in the list ...
                   if(match_id[match_index[j_halo]]==i_halo && j_halo==(n_halos_1_matches-1)){
                      // Check to see if this halo is already in the list ...
                      for(k_halo=0,flag_continue=TRUE;k_halo<halos_i[i_halo].n_bridges && flag_continue;k_halo++){
                         if(bridges[k_halo].halo->id==halos[j_file_1%n_wrap][match_index[j_halo]].id)
                            flag_continue=FALSE;
                      }
                      // ... if not, add it
                      if(flag_continue){
                         bridges[halos_i[i_halo].n_bridges].score=match_score[match_index[j_halo]];
                         bridges[halos_i[i_halo].n_bridges].halo =&(halos[j_file_1%n_wrap][match_index[j_halo]]);
                         (halos_i[i_halo].n_bridges)++;
                      }
                      j_halo++;
                   }
                }
             }
             SID_log("Done.",SID_LOG_CLOSE);
          }

          // ... lastly, reorder the bridges by score, keep only the most immediate bridge descendants and finalize the list ...
          SID_log("Re-ordering bridges...",SID_LOG_OPEN);
          for(i_halo=0;i_halo<n_halos_i;i_halo++){
             if((halos_i[i_halo].n_bridges)>1){ 

                // We may need to remove several halos from the list.  This array will keep track of this.
                bridge_keep=(int *)SID_malloc(sizeof(int)*halos_i[i_halo].n_bridges);

                // Reorder the bridges by their score.  We make a temporary copy of the list 
                //   to do this and initially set all bridges as halos to keep..
                bridges=(bridge_info *)SID_calloc(sizeof(bridge_info)*(halos_i[i_halo].n_bridges));
                for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;j_halo++){
                   bridge=&(halos_i[i_halo].bridges[j_halo]);
                   memcpy(&(bridges[j_halo]),bridge,sizeof(bridge_info));
                   match_score[j_halo]=bridge->score;
                   bridge_keep[j_halo]=TRUE;
                }
                merge_sort((void *)match_score,(size_t)(halos_i[i_halo].n_bridges),&bridge_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

                // Remove any mutual descendants from the list
                //   (since they have their own IDs, this is 
                //    needed to avoid calling them emerged halos)
                for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;j_halo++){
                   bridge = &(bridges[bridge_index[j_halo]]);
                   tree_horizontal_info *current;
                   // ... walk the tree upwards ...
                   current=bridge->halo->descendant.halo;
                   if(current!=NULL)
                      k_file=current->file;
                   l_file=k_file;
                   while(current!=NULL && k_file>=l_file && k_file<MIN(n_files,i_file+(n_search+1))){
                      for(k_halo=0;k_halo<halos_i[i_halo].n_bridges;k_halo++){
                         bridge = &(bridges[bridge_index[k_halo]]);
                         if(bridge->halo==current)
                            bridge_keep[k_halo]=FALSE;
                      }
                      current=current->descendant.halo;
                      l_file=k_file;
                      if(current!=NULL)
                         k_file =current->file;
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
                   (halos_i[i_halo].bridges)=(bridge_info *)SID_calloc(sizeof(bridge_info)*(halos_i[i_halo].n_bridges));

                   // Copy the sorted temporary list to the permanent list.
                   for(j_halo=n_list-1,l_halo=0;j_halo>=0;j_halo--){
                      if(bridge_keep[j_halo]){
                         memcpy(&(halos_i[i_halo].bridges[l_halo]),&(bridges[bridge_index[j_halo]]),sizeof(bridge_info));
                         halos[(bridges[bridge_index[j_halo]].halo->file)%n_wrap][bridges[bridge_index[j_halo]].halo->index].type|=TREE_CASE_EMERGED;
                         if(halos[(bridges[bridge_index[j_halo]].halo->file)%n_wrap][bridges[bridge_index[j_halo]].halo->index].first_progenitor.halo!=NULL)
                            halos[(bridges[bridge_index[j_halo]].halo->file)%n_wrap][bridges[bridge_index[j_halo]].halo->index].type|=TREE_CASE_FOUND;
                         if(halos[(bridges[bridge_index[j_halo]].halo->file)%n_wrap][bridges[bridge_index[j_halo]].halo->index].bridge_backmatch.halo==NULL){
                            halos[(bridges[bridge_index[j_halo]].halo->file)%n_wrap][bridges[bridge_index[j_halo]].halo->index].bridge_backmatch.halo =&(halos_i[i_halo]);
                            halos[(bridges[bridge_index[j_halo]].halo->file)%n_wrap][bridges[bridge_index[j_halo]].halo->index].bridge_backmatch.score=match_score[bridge_index[j_halo]];
                         }
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
          SID_log("Done.",SID_LOG_CLOSE);
          SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
          SID_log("Done.",SID_LOG_CLOSE);
       }

       // Perform forward-matching
       SID_log("Constructing progenitors from forward-matching...",SID_LOG_OPEN|SID_LOG_TIMER);
       for(j_file_1  =i_file,
             j_file_2=i_file+1,
             j_read_1=i_read,
             j_read_2=i_read+i_read_step,
             i_search=0;
           j_read_2<=i_read_stop && i_search<n_search;
           j_file_2++,
             j_read_2+=i_read_step,
             i_search++){

          // Read forward-matching
          read_matches(filename_root_matches,
                       j_read_1,j_read_2,
                       flag_match_subgroups,
                       &n_halos_1_matches,
                       &n_halos_2_matches,
                       n_particles,
                       NULL,
                       n_subgroups_group[j_file_1%n_wrap],
                       n_subgroups_group[j_file_2%n_wrap],
                       match_id,
                       match_score,
                       match_index);

          // Store halo sizes
          if(!flag_fix_bridges){
             if(i_search==0){
                for(i_halo=0;i_halo<n_halos_2_matches;i_halo++)
                   halos[i_file%n_wrap][i_halo].n_particles=n_particles[i_halo];
             }
          }

          
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
                                           n_wrap);
             }
          }

          // Try to match "progenitors" of bridges to the emergent halos identified in their bridge-lists
          for(i_halo=0,n_drop=0;i_halo<n_halos_1_matches;i_halo++){
             if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED)){
                // Loop over all the emergent halos identified with the bridge that i_halo has been matched to
                if(halos_i[i_halo].bridge_forematch.halo==NULL)
                  SID_trap_error("Bridge match not defined during emerged halo search.",ERROR_LOGIC);
                bridges=halos[(halos_i[i_halo].bridge_forematch.halo->file)%n_wrap][halos_i[i_halo].bridge_forematch.halo->index].bridges;
                if(bridges==NULL)
                  SID_trap_error("Bridges not defined during emerged halo search.",ERROR_LOGIC);
                n_list=halos[(halos_i[i_halo].bridge_forematch.halo->file)%n_wrap][halos_i[i_halo].bridge_forematch.halo->index].n_bridges;
                for(k_halo=0;k_halo<n_list;k_halo++){
                   if(bridges[k_halo].halo->file==j_file_2 && match_id[i_halo]==bridges[k_halo].halo->index){
                      set_halo_and_descendant(halos,
                                              i_file,
                                              i_halo,
                                              bridges[k_halo].halo->file,
                                              bridges[k_halo].halo->index,
                                              match_score[i_halo],
                                              &max_id,
                                              n_wrap);
                   }
                }
             }
          }
       }
       SID_log("Done.",SID_LOG_CLOSE);
      
       // Finalize matches to bridges that haven't been associated with any emerged halos
       SID_log("Applying defaults to unprocessed halos...",SID_LOG_OPEN|SID_LOG_TIMER);
       for(i_halo=0,n_drop=0;i_halo<n_halos_1_matches;i_halo++){
          if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED)){
             // The descendant info is already set.  Just increment it's progenitor counter and set this halo's info.
             //   Leave the target flag untouched so we can later identify which BRIDGE_PROGENITORS were found
             halos_i[i_halo].type|=(TREE_CASE_BRIDGE_FINALIZE|TREE_CASE_BRIDGE_DEFAULT);
             set_halo_and_descendant(halos,
                                     i_file,
                                     i_halo,
                                     halos_i[i_halo].bridge_forematch.halo->file,
                                     halos_i[i_halo].bridge_forematch.halo->index,
                                     halos_i[i_halo].bridge_forematch.score,
                                     &max_id,
                                     n_wrap);
          }
       }
      
       // Assign flags for halos not successfully processed.  They must be strays.
       for(i_halo=0;i_halo<n_halos_i;i_halo++){
          if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_UNPROCESSED)){
             halos_i[i_halo].type|=TREE_CASE_STRAYED;
             halos_i[i_halo].type&=(~TREE_CASE_UNPROCESSED);
          }
       }
       SID_log("Done.",SID_LOG_CLOSE);
 
       // Now that we have assigned all the IDs for the halos in the active snapshot,
       //   we need to remove all of their descendants from their bridge list to avoid
       //   matching to the main progenitor's descendants when we should be matching to
       //   it at the current snapshot.
       SID_log("Removing main progenitors from bridged halo lists...",SID_LOG_OPEN|SID_LOG_TIMER);
       sprintf(filename_output_dir_horizontal,                     "%s/horizontal",filename_output_dir);
       sprintf(filename_output_dir_horizontal_groups,              "%s/groups",    filename_output_dir_horizontal);
       sprintf(filename_output_dir_horizontal_groups_matching,     "%s/matching",  filename_output_dir_horizontal_groups);
       sprintf(filename_output_dir_horizontal_subgroups,           "%s/subgroups", filename_output_dir_horizontal);
       sprintf(filename_output_dir_horizontal_subgroups_matching,  "%s/matching",  filename_output_dir_horizontal_subgroups);
       mkdir(filename_output_dir,                              02755);
       mkdir(filename_output_dir_horizontal,                   02755);
       mkdir(filename_output_dir_horizontal_groups,            02755);
       mkdir(filename_output_dir_horizontal_groups_matching,   02755);
       mkdir(filename_output_dir_horizontal_subgroups,         02755);
       mkdir(filename_output_dir_horizontal_subgroups_matching,02755);
       strcpy(filename_output_file_root,filename_output_dir);
       strip_path(filename_output_file_root);
       if(k_match==0)
          sprintf(filename_matching_out,"%s/%s.bridges_%sgroups_%d",filename_output_dir_horizontal_subgroups_matching,filename_output_file_root,group_text_prefix,i_read);
       else
          sprintf(filename_matching_out,"%s/%s.bridges_%sgroups_%d",filename_output_dir_horizontal_groups_matching,filename_output_file_root,group_text_prefix,i_read);
       fp_matching_out=fopen(filename_matching_out,"w");
       i_column=1;
       fprintf(fp_matching_out,"# (%02d): Bridge number\n",                    i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Halo file\n",                        i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Halo index\n",                       i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Halo ID\n",                          i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Tree ID\n",                          i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Descendant file\n",                  i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Descendant index\n",                 i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Descendant ID\n",                    i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Bridge file\n",                      i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Bridge index\n",                     i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Bridge ID\n",                        i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Descendant match type\n",            i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Bridge     match type\n",            i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Descendant match score\n",           i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Bridge     match score\n",           i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Number of particles in halo\n",      i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Number of particles in descendant\n",i_column++,group_text_prefix);
       fprintf(fp_matching_out,"# (%02d): Number of particles in bridge\n",    i_column++,group_text_prefix);
       for(i_halo=0;i_halo<n_halos_i;i_halo++){
          // Check all bridged halos ...
          if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGED)){
             n_list=halos_i[i_halo].n_bridges;
             // ... and check all of their descendants ...
             tree_horizontal_info *current;
             current=halos_i[i_halo].descendant.halo;
           
             if(current!=NULL)
                k_file =current->file;
             l_file=k_file;
             while(current!=NULL && k_file>=l_file && k_file<MIN(n_files,i_file+(n_search+1))){
                for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;){
                   bridge = &(halos_i[i_halo].bridges[j_halo]);
                   if(bridge->halo==current){
                      // ... remove "emerged" flag from halo ...
                      current->type&=(~(TREE_CASE_EMERGED|TREE_CASE_FOUND));
                      // ... do this by decrementing the counter
                      halos_i[i_halo].n_bridges--;
                      // ... and sliding all the halos down ...
                      for(k_halo=j_halo;k_halo<halos_i[i_halo].n_bridges;k_halo++)
                         memcpy(&(halos_i[i_halo].bridges[k_halo]),&(halos_i[i_halo].bridges[k_halo+1]),sizeof(bridge_info));
                   }
                   else
                     j_halo++; // We only need to increment the counter if we don't find a match
                }
                current=current->descendant.halo;
                l_file=k_file;
                if(current!=NULL)
                   k_file =current->file;
             }
            
             // Since we may have removed items, we might not have a bridged halo any more.
             if(halos_i[i_halo].n_bridges<1){
                halos_i[i_halo].type&=(~TREE_CASE_BRIDGED);
                SID_free(SID_FARG halos_i[i_halo].bridges);
                halos_i[i_halo].n_bridges=0;
             }
          }

          // Print bridge info
          if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGED)){
             for(j_halo=0;j_halo<halos_i[i_halo].n_bridges;j_halo++){
                bridge=&(halos_i[i_halo].bridges[j_halo]);
                fprintf(fp_matching_out,"%3d %7d   %4d %7d %7d   %4d %7d %7d   %4d %7d %7d   %8d %8d   %10.4f %10.4f   %7d %7d %7d\n",
                        j_halo,
                        halos_i[i_halo].tree_id,
                        i_read,
                        i_halo,
                        halos_i[i_halo].id,
                        match_file_local(&(halos_i[i_halo].descendant)),
                        match_index_local(&(halos_i[i_halo].descendant)),
                        match_id_local(&(halos_i[i_halo].descendant)),
                        match_file_local(bridge),
                        match_index_local(bridge),
                        match_id_local(bridge),
                        halos_i[i_halo].type,
                        match_type_local(bridge),
                        match_score_local(&(halos_i[i_halo].descendant)),
                        match_score_local(bridge),
                        halos_i[i_halo].n_particles,
                        match_n_particles_local(&(halos_i[i_halo].descendant)),
                        match_n_particles_local(bridge));
             }
          }
          else{
                fprintf(fp_matching_out,"%3d %7d   %4d %7d %7d   %4d %7d %7d   %4d %7d %7d   %8d %8d   %10.4f %10.4f   %7d %7d %7d\n",
                        -1,
                        halos_i[i_halo].tree_id,
                        i_read,
                        i_halo,
                        halos_i[i_halo].id,
                        match_file_local(&(halos_i[i_halo].descendant)),
                        match_index_local(&(halos_i[i_halo].descendant)),
                        match_id_local(&(halos_i[i_halo].descendant)),
                        -1,
                        -1,
                        -1,
                        halos_i[i_halo].type,
                        -1,
                        match_score_local(&(halos_i[i_halo].descendant)),
                        -1.0,
                        halos_i[i_halo].n_particles,
                        match_n_particles_local(&(halos_i[i_halo].descendant)),
                        -1);
          }
       }
       fclose(fp_matching_out);
       SID_log("Done.",SID_LOG_CLOSE);

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
          SID_log("# of emerged progenitors=%-7d (largest=%d particles)",SID_LOG_COMMENT,
             stats.n_emerged_progenitors,stats.max_emerged_progenitor_size);
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
    if(j_file>n_search){
       write_trees_horizontal(groups,   n_groups[l_write],   n_groups_max,   
                              subgroups,n_subgroups[l_write],n_subgroups_max,
                              n_subgroups_group,
                              n_halos_max,
                              max_tree_id_subgroup,
                              max_tree_id_group,
                              i_write,
                              j_write,
                              l_write,
                              n_wrap,
                              i_file_start,
                              filename_cat_root_in,
                              filename_output_dir,
                              a_list,
                              cosmo,
                              n_k_match);
       i_write--;
       l_write++;
       j_write-=i_read_step;
    }
//if(i_read<=828) SID_exit(ERROR_NONE);  

    SID_log("Done.",SID_LOG_CLOSE);
  } // loop over snaps

  // Write the remaining snapshots
  for(;j_write>=i_read_start;i_write--,j_write-=i_read_step,l_write++)
     write_trees_horizontal(groups,   n_groups[l_write],   n_groups_max,   
                            subgroups,n_subgroups[l_write],n_subgroups_max,
                            n_subgroups_group,
                            n_halos_max,
                            max_tree_id_subgroup,
                            max_tree_id_group,
                            i_write,
                            j_write,
                            l_write,
                            n_wrap,
                            i_file_start,
                            filename_cat_root_in,
                            filename_output_dir,
                            a_list,
                            cosmo,
                            n_k_match);
  
  SID_log("Freeing arrays...",SID_LOG_OPEN);
  SID_free(SID_FARG match_id);
  SID_free(SID_FARG match_score);
  SID_free(SID_FARG match_index);
  for(i_search=0;i_search<n_wrap;i_search++){
     // Free subgroup information
     for(i_halo=0;i_halo<n_subgroups_max;i_halo++)
        SID_free(SID_FARG subgroups[i_search][i_halo].bridges);
     SID_free(SID_FARG subgroups[i_search]);

     // Free group information
     SID_free(SID_FARG n_subgroups_group[i_search]);
     for(i_halo=0;i_halo<n_groups_max;i_halo++)
        SID_free(SID_FARG groups[i_search][i_halo].bridges);
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
  //(*flag_clean)=TRUE;

  SID_log("Done.",SID_LOG_CLOSE);
}
