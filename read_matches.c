#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void read_matches(char    *filename_root_matches,
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
   int    k_read;
   int    l_read;
   int    n_search;
   int    n_matches;
   size_t offset;
   int    flag_continue;
   int    i_read_file;
   int    j_read_file;
   int    n_groups_file;
   int    n_groups_file_1;
   int    n_groups_file_2;
   int    n_groups;
   int    n_groups_i_file;
   int    n_groups_j_file;
   
   if(i_read==j_read)
     SID_trap_error("i_read=j_read in read_matches",ERROR_LOGIC);

   switch(mode){
      case MATCH_SUBGROUPS:
      sprintf(group_text_prefix,"sub");
      break;
      case MATCH_GROUPS:
      sprintf(group_text_prefix,"");
      break;
   }

   // Read the needed info from the header file
   int i_read_start;
   int i_read_stop;
   int n_search_total;
   int n_files;
   int i_read_in;
   int n_groups_in;
   int counter=0;
   sprintf(filename_in,"%s.%sgroup_matches_header",filename_root_matches,group_text_prefix);
   SID_fopen(filename_in,"r",&fp_in);
   SID_fread(&i_read_start,  sizeof(int),1,&fp_in);
   SID_fread(&i_read_stop,   sizeof(int),1,&fp_in);
   SID_fread(&n_search_total,sizeof(int),1,&fp_in);
   SID_fread(&n_files,       sizeof(int),1,&fp_in);
   for(i_read=i_read_start;i_read<=i_read_stop && counter<2;i_read++){
      SID_fread(&i_read_in, sizeof(int),1,&fp_in);
      if(i_read==i_read_in){
         SID_fread(n_groups_i,sizeof(int),1,&fp_in);
         if((*n_groups_i)>0){
            SID_fread_ordered(n_particles_i,sizeof(int),(size_t)(*n_groups_i),&fp_in);
            if(mode==MATCH_GROUPS)
               SID_fread_ordered(n_sub_group_i,sizeof(int),(size_t)(*n_groups_i),&fp_in);
         }
         counter++;
      }
      else if(j_read==i_read_in){
         SID_fread(n_groups_j,sizeof(int),1,&fp_in);
         if((*n_groups_j)>0){
            SID_fread_ordered(n_particles_j,sizeof(int),(size_t)(*n_groups_j),&fp_in);
            if(mode==MATCH_GROUPS)
               SID_fread_ordered(n_sub_group_j,sizeof(int),(size_t)(*n_groups_j),&fp_in);
         }
         counter++;
      }
      else{
         SID_fread(&n_groups_in,sizeof(int),1,&fp_in);
         if(n_groups_in>0){
            SID_fskip(sizeof(int),n_groups_in,&fp_in);
            if(mode==MATCH_GROUPS)
               SID_fskip(sizeof(int),n_groups_in,&fp_in);
         }
      }
   }
   SID_fclose(&fp_in);

   // Read the matching file
   sprintf(filename_in,"%s_%03d_%03d.%sgroup_matches",filename_root_matches,i_read,j_read,group_text_prefix);
   SID_fopen(filename_in,"r",&fp_in);
   SID_fread(&i_read_file,sizeof(int),1,&fp_in);
   SID_fread(&j_read_file,sizeof(int),1,&fp_in);
   SID_fread(n_groups_i,  sizeof(int),1,&fp_in);
   SID_fread(n_groups_j,  sizeof(int),1,&fp_in);

   // Read matching data
   SID_fread(match_ids,  sizeof(int),   (*n_groups_i),&fp_in);
   SID_fread(match_index,sizeof(size_t),(*n_groups_i),&fp_in);
   SID_fread(match_score,sizeof(float), (*n_groups_i),&fp_in);
   SID_fclose(&fp_in);
}

