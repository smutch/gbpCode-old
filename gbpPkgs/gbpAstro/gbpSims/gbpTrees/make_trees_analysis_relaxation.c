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
typedef struct process_trees_params_local process_trees_params_local;
struct process_trees_params_local{
   char         filename_output_root[MAX_FILENAME_LENGTH];
   int          snap_lo_tau_trends;
   int          snap_hi_tau_trends;
   trend_info  *mass_binning;
   trend_info **trends_z;
   trend_info **trends_tau_form;
   trend_info **trends_tau_3to1;
   trend_info **trends_tau_10to1;
};

// ** Define the calculation here **
void process_trees_fctn_init_local(tree_info *trees,void *params_in,int mode,int i_type);
void process_trees_fctn_init_local(tree_info *trees,void *params_in,int mode,int i_type){
  // Create an alias for the passed void pointer
  process_trees_params_local *params=(process_trees_params_local *)params_in;
  // Initialize markers (just one-at-a-time to save RAM)
  if(i_type==0)
     precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
  else
     precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
  // Create the trend which will describe the mass binning
  init_treenode_trend(trees,&(params->mass_binning),"logM_course");
  // Allocate the trends
  int n_M=params->mass_binning->ordinate->hist->n_bins;
  params->trends_z        =(trend_info **)SID_malloc(sizeof(trend_info *)*n_M);
  params->trends_tau_form =(trend_info **)SID_malloc(sizeof(trend_info *)*n_M);
  params->trends_tau_3to1 =(trend_info **)SID_malloc(sizeof(trend_info *)*n_M);
  params->trends_tau_10to1=(trend_info **)SID_malloc(sizeof(trend_info *)*n_M);
  // Initialize the trend(s) to be populated
  for(int i_M=0;i_M<n_M;i_M++){
     // f(z) trends
     init_treenode_trend           (trees,&(params->trends_z[i_M]),"z");
     init_treenode_trend_coordinate(trees, (params->trends_z[i_M]),"xoff");
     init_treenode_trend_coordinate(trees, (params->trends_z[i_M]),"tau_form");
     init_treenode_trend_coordinate(trees, (params->trends_z[i_M]),"tau_3to1");
     init_treenode_trend_coordinate(trees, (params->trends_z[i_M]),"tau_10to1");
     if(i_type==0)
        init_treenode_trend_coordinate(trees, (params->trends_z[i_M]),"SSFctn");
     // f(tau_form) trends
     init_treenode_trend           (trees,&(params->trends_tau_form[i_M]),"tau_form");
     init_treenode_trend_coordinate(trees, (params->trends_tau_form[i_M]),"xoff");
     if(i_type==0)
        init_treenode_trend_coordinate(trees, (params->trends_tau_form[i_M]),"SSFctn");
     // f(tau_3to1) trends
     init_treenode_trend           (trees,&(params->trends_tau_3to1[i_M]),"tau_3to1");
     init_treenode_trend_coordinate(trees, (params->trends_tau_3to1[i_M]),"xoff");
     if(i_type==0)
        init_treenode_trend_coordinate(trees, (params->trends_tau_3to1[i_M]), "SSFctn");
     // f(tau_10to1) trends
     init_treenode_trend           (trees,&(params->trends_tau_10to1[i_M]),"tau_10to1");
     init_treenode_trend_coordinate(trees, (params->trends_tau_10to1[i_M]),"xoff");
     if(i_type==0)
        init_treenode_trend_coordinate(trees, (params->trends_tau_10to1[i_M]),"SSFctn");
  }
}

// ** Perform the calculation here **
void process_trees_fctn_analyze_local(tree_info *trees,void *params,int mode,int i_type,int flag_init,tree_node_info *halo);
void process_trees_fctn_analyze_local(tree_info *trees,void *params,int mode,int i_type,int flag_init,tree_node_info *halo){
   int i_M;
   int i_snap   =halo->snap_tree;
   int i_snap_lo=((process_trees_params_local *)params)->snap_lo_tau_trends;
   int i_snap_hi=((process_trees_params_local *)params)->snap_hi_tau_trends;
   if((i_M=add_item_to_trend(((process_trees_params_local *)params)->mass_binning,GBP_ADD_ITEM_TO_TREND_DEFAULT,halo))>=0){
      add_item_to_trend(((process_trees_params_local *)params)->trends_z[i_M],GBP_ADD_ITEM_TO_TREND_DEFAULT,halo);
      if(i_snap>=i_snap_lo&&i_snap<=i_snap_hi){
         add_item_to_trend(((process_trees_params_local *)params)->trends_tau_form[i_M], GBP_ADD_ITEM_TO_TREND_DEFAULT,halo);
         add_item_to_trend(((process_trees_params_local *)params)->trends_tau_3to1[i_M], GBP_ADD_ITEM_TO_TREND_DEFAULT,halo);
         add_item_to_trend(((process_trees_params_local *)params)->trends_tau_10to1[i_M],GBP_ADD_ITEM_TO_TREND_DEFAULT,halo);
      }
   }
}

// ** Write the results here **
void process_trees_fctn_fin_local(tree_info *trees,void *params_in,int mode,int i_type);
void process_trees_fctn_fin_local(tree_info *trees,void *params_in,int mode,int i_type){
   // Create an alias for the passed void pointer
   process_trees_params_local *params=(process_trees_params_local *)params_in;
   // Clean-up markers
   char filename_output_root_root[MAX_FILENAME_LENGTH];
   if(i_type==0){
      sprintf(filename_output_root_root,"%s_groups",params->filename_output_root);
      free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
   }
   else{
      sprintf(filename_output_root_root,"%s_subgroups",params->filename_output_root);
      free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
   }
   // Write results
   finalize_trend(params->mass_binning); // needed for the check on bin_count (below) to work w/ n_proc>1
   int n_M=params->mass_binning->ordinate->hist->n_bins;
   for(int i_M=0;i_M<n_M;i_M++){
      if(params->mass_binning->ordinate->hist->bin_count[i_M]>50){
         // Set the root of the output filename
         char filename_output_root[MAX_FILENAME_LENGTH];
         char mass_text_lo[64];
         char mass_text_hi[64];
         float_to_text((float)histogram_bin_x_lo(params->mass_binning->ordinate->hist,i_M),1,mass_text_lo);
         float_to_text((float)histogram_bin_x_hi(params->mass_binning->ordinate->hist,i_M),1,mass_text_hi);
         sprintf(filename_output_root,"%s_logM_%s_to_%s",filename_output_root_root,mass_text_lo,mass_text_hi);
         // Write the trends
         write_trend_ascii(params->trends_z[i_M],        filename_output_root);
         write_trend_ascii(params->trends_tau_form[i_M], filename_output_root);
         write_trend_ascii(params->trends_tau_3to1[i_M], filename_output_root);
         write_trend_ascii(params->trends_tau_10to1[i_M],filename_output_root);
      }
      // Clean-up the trends
      free_trend(&(params->trends_z[i_M]));
      free_trend(&(params->trends_tau_form[i_M]));
      free_trend(&(params->trends_tau_3to1[i_M]));
      free_trend(&(params->trends_tau_10to1[i_M]));
   }
   SID_free(SID_FARG params->trends_z);
   SID_free(SID_FARG params->trends_tau_form);
   SID_free(SID_FARG params->trends_tau_3to1);
   SID_free(SID_FARG params->trends_tau_10to1);
   // Write and clean-up the mass-binning trend
   write_trend_ascii(params->mass_binning,filename_output_root_root);
   free_trend(&(params->mass_binning));
}

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  char filename_SSimPL_dir[MAX_FILENAME_LENGTH];
  char filename_halo_version_root[MAX_FILENAME_LENGTH];
  char filename_trees_name[MAX_FILENAME_LENGTH];
  char filename_output_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_dir,       argv[1]);
  strcpy(filename_halo_version_root,argv[2]);
  strcpy(filename_trees_name,       argv[3]);
  strcpy(filename_output_root,      argv[4]);
  double z_lo_tau_trends      =atof(argv[5]);
  double z_hi_tau_trends      =atof(argv[6]);

  // Set the halo and tree filename roots
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_dir,filename_trees_name);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_dir,filename_halo_version_root);

  SID_log("Generating treenode markers & analysis of merger trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Perform analysis
  tree_info *trees;
  read_trees(filename_SSimPL_dir,
             filename_halo_version_root,
             filename_trees_name,
             TREE_MODE_DEFAULT,
             &trees);

  // Read catalogs
  read_trees_catalogs(trees,
                      filename_SSimPL_dir,
                      filename_halo_version_root,
                      READ_TREES_CATALOGS_BOTH);

  // Populate the structure that gets passed to the calculation
  process_trees_params_local params;
  strcpy(params.filename_output_root,filename_output_root);
  params.snap_hi_tau_trends=find_treesnap_z(trees,z_lo_tau_trends); // snap_tree and z_list run in opposite orders
  params.snap_lo_tau_trends=find_treesnap_z(trees,z_hi_tau_trends); // snap_tree and z_list run in opposite orders
  SID_log("Trends with tau will be compiled between (z,snap)=(%lf,%d) and (z,snap)=(%lf,%d).",SID_LOG_COMMENT,
          trees->z_list[params.snap_hi_tau_trends],params.snap_hi_tau_trends,
          trees->z_list[params.snap_lo_tau_trends],params.snap_lo_tau_trends);

  // ** PERFORM the calculation here **
  process_trees_by_snap(trees,&params,PROCESS_TREES_BOTH,0,trees->n_snaps,
                        process_trees_fctn_init_local,
                        process_trees_fctn_init_snap_null,
                        process_trees_fctn_select_null,
                        process_trees_fctn_analyze_local,
                        process_trees_fctn_fin_snap_null,
                        process_trees_fctn_fin_local);

  // Clean-up
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

