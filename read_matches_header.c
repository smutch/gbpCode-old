#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int read_matches_header(char   *filename_root_in,
                        char   *filename_root_out,
                        int     i_read_start,
                        int     i_read_stop,
                        int     i_read_step,
                        int    *n_files_return,
                        int   **n_subgroups_return,
                        int   **n_groups_return,
                        int     n_search){
  char        filename_out[256];
  char        filename_out_dir[256];
  char        filename_out_name[256];
  FILE       *fp_test;
  FILE       *fp_out;
  SID_fp      fp_in;
  int         i_read,k_read,l_read;
  int         flag_go;
  int         i_read_start_file;
  int         i_read_stop_file;
  int         i_read_step_file;
  int         n_search_file;
  int         k_match;
  int         n_search_total;
  int         n_k_match;
  int         flag_match_subgroups;
  char        group_text_prefix[5];
  int         n_matches;
  int         n_files;
  int         j_read;
  char        filename_cat1[256];
  char        filename_cat2[256];
  char        filename_cat1_order[256];
  char        filename_cat2_order[256];
  char        filename_out_dir_snap[256];
  int         n_groups_1;
  int         n_groups_1_local;
  int         n_groups_2;
  int         n_groups_2_local;
  int         i_group;
  int         buffered_count;
  int         buffered_count_local;
  int         j_group;
  int         index_test;
  int         i_rank;
  int        *n_particles;
  int        *n_sub_group;
  int        *match_id   =NULL;
  float      *match_score=NULL;
  char        cat_name_1[20];
  char        cat_name_2[20];
  size_t     *match_rank =NULL;
  size_t     *match_index=NULL;
  size_t      offset;
  plist_info  plist1;
  plist_info  plist2;
  void       *buffer;
  int        *buffer_int;
  size_t     *buffer_size_t;
  float      *buffer_float;
  int         n_buffer_max=1000;
  int         n_buffer;
  int         i_buffer;
  int         j_buffer;
  int         flag_compute_header_subgroups;
  int         flag_compute_header_groups;
  int         flag_compute_header;
  int         flag_sucessful_completion=TRUE;

  SID_log("Reading header information...",SID_LOG_OPEN);

  // Count the number of snapshots we are going to use and
  //    initialize the arrays that are to be returned
  for(i_read=i_read_stop,(*n_files_return)=0;      
      i_read>=i_read_start;
      i_read-=i_read_step) (*n_files_return)++;  
  (*n_subgroups_return)=(int *)SID_malloc(sizeof(int)*(*n_files_return));
  (*n_groups_return)   =(int *)SID_malloc(sizeof(int)*(*n_files_return));

  if(SID.I_am_Master){
     FILE *fp_read_header;
     for(k_match=0;k_match<2;k_match++){
        switch(k_match){
           case 0:
           sprintf(group_text_prefix,"sub");
           flag_compute_header=flag_compute_header_subgroups;
           break;
           case 1:
           sprintf(group_text_prefix,"");
           flag_compute_header=flag_compute_header_groups;
           break;
        }
        for(i_read=i_read_stop,j_read=0;i_read>=i_read_start;j_read++){
           // Open file and skip header           
           if(i_read==i_read_stop){
              sprintf(filename_out,"%s/%sgroup_matches_header.dat",filename_out_dir,group_text_prefix);
              if((fp_read_header=fopen(filename_out,"r"))==NULL)
                 SID_trap_error("Could not open file {%s} when reading header information.",ERROR_IO_OPEN,filename_out);
              fseek(fp_read_header,4*sizeof(int),SEEK_SET);
           }

           // Read-forward for the appropriate number of snapshots
           int k_read;
           for(k_read=0;i_read>=i_read_start && k_read<i_read_step;i_read--,k_read++){
              fseek(fp_read_header,1*sizeof(int),SEEK_CUR);
              fread(&n_groups_1,sizeof(int),1,fp_read_header);
              fseek(fp_read_header,n_groups_1*sizeof(int),SEEK_CUR);
              if(k_match==1)
                 fseek(fp_read_header,n_groups_1*sizeof(int),SEEK_CUR);
              if(k_read==0){
                 switch(k_match){
                    case 0:
                    (*n_subgroups_return)[j_read]=n_groups_1;
                    break;
                    case 1:
                    (*n_groups_return)[j_read]   =n_groups_1;
                    break;
                 }
              }
           }
        }
        fclose(fp_read_header);
     }
  }
  SID_Bcast((*n_subgroups_return),sizeof(int)*(*n_files_return),MASTER_RANK,SID.COMM_WORLD);
  SID_Bcast((*n_groups_return),   sizeof(int)*(*n_files_return),MASTER_RANK,SID.COMM_WORLD);
  SID_log("Done.",SID_LOG_CLOSE);

  return(flag_sucessful_completion);
}

