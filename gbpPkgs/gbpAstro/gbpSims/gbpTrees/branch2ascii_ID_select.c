#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs -- filenames
  char filename_SSimPL_root[MAX_FILENAME_LENGTH];
  char filename_halos_version[MAX_FILENAME_LENGTH];
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_trees_version[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_root,  argv[1]);
  strcpy(filename_halos_version,argv[2]);
  strcpy(filename_trees_version,argv[3]);

  // Fetch user inputs -- halo ID list
  int n_list=argc-4;
  if(n_list<1)
     SID_trap_error("You have not listed any halo_IDs to report on.",ERROR_SYNTAX);
  int *halo_ID_list=(int *)SID_malloc(sizeof(int)*n_list);
  for(int i_list=0;i_list<n_list;i_list++)
     halo_ID_list[i_list]=atoi(argv[4+i_list]);
  merge_sort(halo_ID_list,n_list,NULL,SID_INT,SORT_INPLACE_ONLY,SORT_COMPUTE_INPLACE);
  SID_log("Selected halo IDs:",SID_LOG_OPEN);
  for(int i_list=0;i_list<n_list;i_list++)
     SID_log(" %d",SID_LOG_CONTINUE,halo_ID_list[i_list]);
  SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);

  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_root,filename_trees_version);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_root,filename_halos_version);

  SID_log("Writing ascii versions of %d main progenitor branch(s)...",SID_LOG_OPEN|SID_LOG_TIMER,n_list);

  // Perform analysis
  tree_info *trees;
  read_trees(filename_SSimPL_root,
             filename_halos_version,
             filename_trees_version,
             TREE_MODE_DEFAULT,
             &trees);

  // Read ancillary data
  read_trees_catalogs(trees,
                      filename_SSimPL_root,
                      filename_halos_version,
                      READ_TREES_CATALOGS_GROUPS|READ_TREES_CATALOGS_SUBGROUPS);

  // Initialize a tree_node array large enough for all the requested halo IDs
  tree_node_info **selected =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_list);
  for(int i_list=0;i_list<n_list;i_list++) selected[i_list]=NULL;

  // Find the root tree_nodes of the given halo IDs
  int n_unfound=n_list;
  for(int i_forest=(trees->n_forests_local-1);i_forest>=0 && n_unfound>0;i_forest--){
     tree_node_info *current_halo=trees->first_in_forest_subgroups[i_forest];
     while(current_halo && n_unfound>0){
        for(int i_list=0;i_list<n_list;i_list++){
           if(current_halo->halo_ID==halo_ID_list[i_list] && selected[i_list]==NULL){
              find_treenode_branch_root(trees,current_halo,&(selected[i_list]));
              if(selected[i_list]!=NULL)
                 n_unfound--;
           }
        }
        current_halo=current_halo->next_in_forest;
     }
  }

  // Write file
  for(int i_list=0;i_list<n_list;i_list++){
     char  filename_out[MAX_FILENAME_LENGTH];
     sprintf(filename_out,"subhalo_ID_%09d_branch.txt",halo_ID_list[i_list]);
     write_tree_branch_ascii(trees,selected[i_list],filename_out,filename_trees_version);
  }

  // Clean-up
  free_trees(&trees);
  SID_free(SID_FARG halo_ID_list);
  SID_free(SID_FARG selected);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

