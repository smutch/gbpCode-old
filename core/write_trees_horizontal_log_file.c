#include <gbpTrees_build.h>

void write_trees_horizontal_log_file(char *filename_log,int l_write,int j_write,int i_k_match,int n_k_match,tree_horizontal_stats_info *stats,double *a_list,cosmo_info **cosmo,int flag_init){
   FILE *fp;
   char  group_text_prefix[5];
   if(flag_init){
      int   i_column;
      fp=fopen(filename_log,"w");
      i_column=1;
      fprintf(fp,"# (%02d): Expansion factor (a)\n",   i_column++);
      fprintf(fp,"# (%02d): Snapshot filenumber\n",    i_column++);
      fprintf(fp,"# (%02d): Snapshot interval [yrs]\n",i_column++);
      int j_k_match;
      for(j_k_match=0;j_k_match<n_k_match;j_k_match++){
         switch(j_k_match){
            case 0:
            sprintf(group_text_prefix,"sub");
            break;
            case 1:
            sprintf(group_text_prefix,"");
            break;
         }
         fprintf(fp,"# (%02d): Maximum %sgroup ID\n",                       i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): # of %sgroups\n",                            i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): # of simple        %sgroups\n",              i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): # of merging       %sgroups\n",              i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): # of strayed       %sgroups\n",              i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): # of dropped       %sgroups\n",              i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): # of bridged       %sgroups\n",              i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): # of matches to    %sgroup  bridges\n",      i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): # of emerged       %sgroups\n",              i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): # of fragmented    %sgroups strayed\n",      i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): # of fragmented    %sgroups returned\n",     i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): # of fragmented    %sgroups exchanged\n",    i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): # of matches to    %sgroup  emerged halos\n",i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): Largest strayed    %sgroup\n",               i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): Largest dropped    %sgroup\n",               i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): Largest bridged    %sgroup\n",               i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): Largest emerged    %sgroup\n",               i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): Largest strayed   fragmented %sgroup\n",     i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): Largest returned  fragmented %sgroup\n",     i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): Largest exchanged fragmented %sgroup\n",     i_column++,group_text_prefix);
         fprintf(fp,"# (%02d): Largest emerged    %sgroup  progenitor\n",   i_column++,group_text_prefix);
      }
      fclose(fp);
   }

   switch(i_k_match){
      case 0:
      sprintf(group_text_prefix,"sub");
      break;
      case 1:
      sprintf(group_text_prefix,"");
      break;
   }

   fp=fopen(filename_log,"a");
   if(i_k_match==0){
      if(l_write>0)
        fprintf(fp,"%le %4d %10.4lf",a_list[l_write],j_write,deltat_a(cosmo,a_list[l_write],a_list[l_write-1])/S_PER_YEAR);
      else
        fprintf(fp,"%le %4d %10.4lf",a_list[l_write],j_write,deltat_a(cosmo,a_list[l_write+1],a_list[l_write])/S_PER_YEAR);
   }
   fprintf(fp," %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d %08d",
           stats->max_id,
           stats->n_halos,
           stats->n_simple,
           stats->n_mergers,
           stats->n_strayed,
           stats->n_dropped,
           stats->n_bridged,
           stats->n_bridge_progenitors,
           stats->n_emerged,
           stats->n_fragmented_strayed,
           stats->n_fragmented_returned,
           stats->n_fragmented_exchanged,
           stats->n_emerged_progenitors,
           stats->max_strayed_size,
           stats->max_dropped_size,
           stats->max_bridged_size,
           stats->max_emerged_size,
           stats->max_fragmented_strayed_size,
           stats->max_fragmented_returned_size,
           stats->max_fragmented_exchanged_size,
           stats->max_emerged_progenitor_size);
   if(i_k_match==n_k_match-1)
      fprintf(fp,"\n");
   fclose(fp);
}

