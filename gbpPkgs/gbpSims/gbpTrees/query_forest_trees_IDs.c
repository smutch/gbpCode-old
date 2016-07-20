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
  int     n_files;
  int     k_read;
  int     max_n_groups;
  int     l_read;
  int     n_groups;
  int    *n_particles_i;
  int    *n_particles_j;
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

  // Fetch user inputs
  if(argc!=6)
     SID_trap_error("Invalid Syntax.",ERROR_SYNTAX);
  char    filename_SSimPL_root[MAX_FILENAME_LENGTH];
  char    filename_halos_root[MAX_FILENAME_LENGTH];
  char    filename_trees_root[MAX_FILENAME_LENGTH];
  char    halo_type_prefix_text[8];
  strcpy(filename_SSimPL_root,argv[1]);
  strcpy(filename_halos_root, argv[2]);
  strcpy(filename_trees_root, argv[3]);
  if(!strcmp(argv[4],"groups") || !strcmp(argv[4],"group"))
     mode=MATCH_GROUPS;
  else if(!strcmp(argv[4],"subgroups") || !strcmp(argv[4],"subgroup"))
     mode=MATCH_SUBGROUPS;
  else{
     SID_log("Invalid mode selection {%s}.  Should be 'group' or 'subgroup'.",SID_LOG_COMMENT,argv[4]);
     SID_exit(ERROR_SYNTAX);
  }
  int halo_forest_find=atoi(argv[5]);

  SID_log("Querying forest tree IDs...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Read tree header information
  tree_info *trees;
  char       filename_file_root[MAX_FILENAME_LENGTH];
  sprintf(filename_file_root,"%s/trees/%s",filename_SSimPL_root,filename_trees_root);
  init_trees_read(filename_SSimPL_root,filename_halos_root,filename_trees_root,TREE_READ_DEFAULT,&trees);

  if(mode==MATCH_GROUPS){
     SID_log("The following group tree IDs are in forest No. %d:",SID_LOG_COMMENT,halo_forest_find);
     for(int i_halo=0;i_halo<trees->n_trees_group;i_halo++){
        if(trees->tree2forest_mapping_group[i_halo]==halo_forest_find){
           fprintf(stdout,"%d\n",i_halo);
        }
     }
  }
  else if(mode==MATCH_SUBGROUPS){
     SID_log("The following subgroup tree IDs are in forest No. %d:",SID_LOG_COMMENT,halo_forest_find);
     for(int i_halo=0;i_halo<trees->n_trees_subgroup;i_halo++){
        if(trees->tree2forest_mapping_subgroup[i_halo]==halo_forest_find){
           fprintf(stdout,"%d\n",i_halo);
        }
     }
  }
  else
     SID_trap_error("Invalid mode (%d).",ERROR_LOGIC,mode);

  // Clean-up
  free_trees(&trees);
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}


