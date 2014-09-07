#define _MAIN
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

// Structure that will carry the needed information to the select-and-analyze function
typedef struct select_and_analyze_params_local select_and_analyze_params_local;
struct select_and_analyze_params_local{
   char        filename_output_root[MAX_FILENAME_LENGTH];
   trend_info *trends_z;
};

// ** Define the calculation here **
void select_and_analyze_treenodes_fctn_init_local(tree_info *trees,void *params_in,int mode,int i_type);
void select_and_analyze_treenodes_fctn_init_local(tree_info *trees,void *params_in,int mode,int i_type){

  // Create an alias for the passed void pointer
  select_and_analyze_params_local *params=(select_and_analyze_params_local *)params_in;

  // Initialize markers (just one-at-a-time to save RAM)
  if(i_type==0){
     //read_treenode_markers(trees,params->filename_output_root,PRECOMPUTE_TREENODE_MARKER_GROUPS);
     precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
  }
  else{
     //read_treenode_markers(trees,params->filename_output_root,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
     precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
  }

  // Initialize the trend(s) to be populated
  init_trend           (&(params->trends_z));
  init_trend_ordinate  ( (params->trends_z),"z",        trees,init_tree_property_z,     free_tree_property_z,     calc_tree_property_index_z);
  init_trend_coordinate( (params->trends_z),"logM",     trees,init_tree_property_logM,  free_tree_property_logM,  calc_tree_property_index_logM);
  init_trend_coordinate( (params->trends_z),"SSFctn",   trees,init_tree_property_SSFctn,free_tree_property_SSFctn,calc_tree_property_index_SSFctn);
  init_trend_coordinate( (params->trends_z),"tau_form", trees,init_tree_property_tau,   free_tree_property_tau,   calc_tree_property_index_tau_form);
  init_trend_coordinate( (params->trends_z),"tau_3to1", trees,init_tree_property_tau,   free_tree_property_tau,   calc_tree_property_index_tau_3to1);
  init_trend_coordinate( (params->trends_z),"tau_10to1",trees,init_tree_property_tau,   free_tree_property_tau,   calc_tree_property_index_tau_10to1);
}

// ** Perform the calculation here **
void select_and_analyze_treenodes_fctn_analyze_local(tree_info *trees,void *params,int mode,int i_type,int flag_init,tree_node_info *halo);
void select_and_analyze_treenodes_fctn_analyze_local(tree_info *trees,void *params,int mode,int i_type,int flag_init,tree_node_info *halo){
   add_item_to_trend(((select_and_analyze_params_local *)params)->trends_z,halo);
}

// ** Write the results here **
void select_and_analyze_treenodes_fctn_fin_local(tree_info *trees,void *params_in,int mode,int i_type);
void select_and_analyze_treenodes_fctn_fin_local(tree_info *trees,void *params_in,int mode,int i_type){
   // Create an alias for the passed void pointer
   select_and_analyze_params_local *params=(select_and_analyze_params_local *)params_in;

   // Clean-up markers
   char filename_output_root[MAX_FILENAME_LENGTH];
   if(i_type==0){
      sprintf(filename_output_root,"%s_groups",params->filename_output_root);
      free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
   }
   else{
      sprintf(filename_output_root,"%s_subgroups",params->filename_output_root);
      free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
   }

   // Write results
   write_trend_ascii(params->trends_z,filename_output_root);

   // Clean-up trend
   free_trend(&(params->trends_z));
}

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  char filename_SSimPL_dir[MAX_FILENAME_LENGTH];
  char filename_halo_version_root[MAX_FILENAME_LENGTH];
  char filename_trees_name[MAX_FILENAME_LENGTH];
  char filename_output_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_dir,       argv[1]);
  strcpy(filename_halo_version_root,argv[2]);
  strcpy(filename_trees_name,       argv[3]);
  strcpy(filename_output_root,      argv[4]);

  // Set the halo and tree filename roots
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_dir,filename_trees_name);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_dir,filename_halo_version_root);

  SID_log("Generating treenode markers & analysis of merger trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Perform analysis
  tree_info *trees;
  read_trees(filename_trees_root,
             filename_halos_root,
             TREE_MODE_DEFAULT,
             &trees);

  // Read catalogs
  read_trees_catalogs(trees,
                      filename_SSimPL_dir,
                      filename_halo_version_root,
                      READ_TREES_CATALOGS_BOTH);

  // ** PERFORM the calculation here **
  select_and_analyze_params_local params;
  strcpy(params.filename_output_root,filename_output_root);
  select_and_analyze_treenodes_by_snap(trees,&params,SELECT_AND_ANALYZE_BOTH,0,trees->n_snaps,
                                       select_and_analyze_treenodes_fctn_init_local,
                                       select_and_analyze_treenodes_fctn_init_snap_null,
                                       select_and_analyze_treenodes_fctn_select_null,
                                       select_and_analyze_treenodes_fctn_analyze_local,
                                       select_and_analyze_treenodes_fctn_fin_snap_null,
                                       select_and_analyze_treenodes_fctn_fin_local);

  // Clean-up
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

