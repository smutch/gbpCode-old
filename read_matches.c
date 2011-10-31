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
         if(n_particles_i!=NULL)
            SID_fread(n_particles_i,sizeof(int),n_groups,&fp_in);
         else
            SID_fseek(&fp_in,sizeof(int),n_groups,SID_SEEK_CUR);
         if(mode==MATCH_GROUPS){
            if(n_sub_group_i!=NULL)
               SID_fread(n_sub_group_i,sizeof(int),n_groups,&fp_in);
            else
               SID_fseek(&fp_in,sizeof(int),n_groups,SID_SEEK_CUR);
         }
      }
      else if(j_read==l_read){
         n_groups_j_file=n_groups;
         if(n_particles_j!=NULL)
            SID_fread(n_particles_j,sizeof(int),n_groups,&fp_in);
         else
            SID_fseek(&fp_in,sizeof(int),n_groups,SID_SEEK_CUR);
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
