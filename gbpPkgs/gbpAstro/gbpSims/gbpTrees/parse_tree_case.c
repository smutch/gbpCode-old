#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){
  int     n_search;
  int     i_halo;
  char    filename_root_in[MAX_FILENAME_LENGTH];
  char    filename_in[MAX_FILENAME_LENGTH];
  char    group_text_prefix[4];
  int     n_files;
  int     k_read;
  int     max_n_groups;
  int     l_read;
  int     n_groups;
  int    *n_particles_i;
  int    *n_particles_j;
  int     j_read;
  int     mode;
  int     n_groups_i;
  int     n_groups_j;
  int    *match_ids;
  float  *match_score;
  size_t *match_index;
  int     j_halo;
  int     match;
  int     i_read;
  int     i_read_start;
  int     i_read_stop;
  SID_fp  fp_in;

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user input
  int match_type;
  match_type=atoi(argv[1]);
  SID_log("Parsing tree_case_flag_list for match_type=%d...",SID_LOG_OPEN,match_type);

  // Perform parsing
  int i_parse=0;
  int count  =0;
  int test_restult;
  for(i_parse=0;i_parse<n_tree_case_flag_list;i_parse++){
     if(check_mode_for_flag(match_type,tree_case_flag_list[i_parse])){
        count++;
        if(count==1)
           SID_log("The following tree_case_flag_list are switched ON:",SID_LOG_OPEN);
        SID_log("%s",SID_LOG_COMMENT,tree_case_flag_list_text[i_parse]);
     }
  }
  if(count==0)
     SID_log("No tree_case_flag_list are switched ON",SID_LOG_OPEN);
  SID_log("",SID_LOG_SILENT_CLOSE);
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

