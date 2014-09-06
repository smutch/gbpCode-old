#define _MAIN
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

typedef struct treenode_trend_property_info treenode_trend_property_info;
struct treenode_trend_property_info{
   char       name[128];
   int        n_ordinate;
   hist_info *hist;
   int      (*calc_index_function)(tree_info *trees,hist_info *hist,tree_node_info *halo);
};

#define INIT_TREENODE_TREND_DEFAULT 0
typedef struct treenode_trend_info treenode_trend_info;
struct treenode_trend_info{
   int                            mode;
   int                            n_properties;
   char                           filename_output_root[MAX_FILENAME_LENGTH];
   treenode_trend_property_info  *ordinate;
   treenode_trend_property_info **property;
};

int calc_z_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo);
int calc_z_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo){
   return(halo->snap_tree);
}

int calc_logM_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo);
int calc_logM_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo){
   return(calc_histogram_index(hist,take_log10(fetch_treenode_Mvir(trees,halo))));
}

int calc_x_off_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo);
int calc_x_off_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo){
   return(calc_histogram_index(hist,fetch_treenode_x_off(trees,halo)));
}

int calc_SSFctn_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo);
int calc_SSFctn_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo){
   return(calc_histogram_index(hist,fetch_treenode_SSFctn(trees,halo)));
}

int calc_tau_form_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo);
int calc_tau_form_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo){
   tree_markers_info *markers=fetch_treenode_precomputed_markers(trees,halo);
   return(fetch_treenode_snap_tree(trees,markers->half_peak_mass));
}

int calc_tau_3to1_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo);
int calc_tau_3to1_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo){
   tree_markers_info *markers=fetch_treenode_precomputed_markers(trees,halo);
   return(fetch_treenode_snap_tree(trees,markers->merger_33pc_remnant));
}

int calc_tau_10to1_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo);
int calc_tau_10to1_hist_index(tree_info *trees,hist_info *hist,tree_node_info *halo){
   tree_markers_info *markers=fetch_treenode_precomputed_markers(trees,halo);
   return(fetch_treenode_snap_tree(trees,markers->merger_10pc_remnant));
}

void add_treenode_to_trend(tree_info *trees,treenode_trend_info *trend,tree_node_info *halo);
void add_treenode_to_trend(tree_info *trees,treenode_trend_info *trend,tree_node_info *halo){
   int i_ordinate=trend->ordinate->calc_index_function(trees,trend->ordinate->hist,halo);
   if(is_histogram_index_in_range(trend->ordinate->hist,i_ordinate)){
      add_to_histogram_index(trend->ordinate->hist,i_ordinate);
      for(int i_property=0;i_property<trend->n_properties;i_property++){
         treenode_trend_property_info *property_i  =trend->property[i_property];
         hist_info                    *hist_i      =&(property_i->hist[i_ordinate]);
         int                           i_coordinate=property_i->calc_index_function(trees,hist_i,halo);
         add_to_histogram_index(hist_i,i_coordinate);
      }
   }
}

void init_treenode_trend_property(tree_info *trees,treenode_trend_property_info **property,const char *name,int n_ordinate);
void init_treenode_trend_property(tree_info *trees,treenode_trend_property_info **property,const char *name,int n_ordinate){
   (*property)=(treenode_trend_property_info *)SID_malloc(sizeof(treenode_trend_property_info));
   strcpy((*property)->name,name);
   (*property)->n_ordinate=n_ordinate;
 
   if(!strcmp(name,"z")){
      (*property)->hist=(hist_info *)SID_malloc(sizeof(hist_info)*((*property)->n_ordinate));
      for(int i_hist=0;i_hist<(*property)->n_ordinate;i_hist++)
         init_histogram(&((*property)->hist[i_hist]),GBP_HISTOGRAM_IRREGULAR_XLO_DEFINED,trees->z_list,trees->n_snaps);
      (*property)->calc_index_function=calc_z_hist_index;
   }
   else if(!strcmp(name,"tau_form") ||
           !strcmp(name,"tau_3to1") ||
           !strcmp(name,"tau_10to1")){
      (*property)->hist       =(hist_info *)SID_malloc(sizeof(hist_info)*((*property)->n_ordinate));
      double *tau_array=(double *)SID_malloc(sizeof(double)*trees->n_snaps);
      for(int i_hist=0;i_hist<(*property)->n_ordinate;i_hist++){
         for(int i_tau=0;i_tau<trees->n_snaps;i_tau++)
            tau_array[i_tau]=10.*((trees->t_list[i_hist]-trees->t_list[i_tau])/trees->t_list[i_hist]);
         init_histogram(&((*property)->hist[i_hist]),GBP_HISTOGRAM_IRREGULAR_XLO_DEFINED,tau_array,trees->n_snaps);
      }
      SID_free(SID_FARG tau_array);
      if(!strcmp(name,"tau_form"))
         (*property)->calc_index_function=calc_tau_form_hist_index;
      else if(!strcmp(name,"tau_3to1"))
         (*property)->calc_index_function=calc_tau_3to1_hist_index;
      else if(!strcmp(name,"tau_10to1"))
         (*property)->calc_index_function=calc_tau_10to1_hist_index;
   }
   else if(!strcmp(name,"logM")){
      double x_min = 7.0;
      double dx    = 0.5;
      int    n_x   =  14;
      (*property)->hist=(hist_info *)SID_malloc(sizeof(hist_info)*((*property)->n_ordinate));
      for(int i_hist=0;i_hist<(*property)->n_ordinate;i_hist++)
         init_histogram(&((*property)->hist[i_hist]),GBP_HISTOGRAM_FIXED,x_min,dx,n_x);
      (*property)->calc_index_function=calc_logM_hist_index;
   }
   else if(!strcmp(name,"x_off")){
      double x_min = 0.0;
      double dx    =0.01;
      int    n_x   =  50;
      (*property)->hist=(hist_info *)SID_malloc(sizeof(hist_info)*((*property)->n_ordinate));
      for(int i_hist=0;i_hist<(*property)->n_ordinate;i_hist++)
         init_histogram(&((*property)->hist[i_hist]),GBP_HISTOGRAM_FIXED,x_min,dx,n_x);
      (*property)->calc_index_function=calc_x_off_hist_index;
   }
   else if(!strcmp(name,"SSFctn")){
      double x_min = 0.0;
      double dx    =0.05;
      int    n_x   =  20;
      (*property)->hist=(hist_info *)SID_malloc(sizeof(hist_info)*((*property)->n_ordinate));
      for(int i_hist=0;i_hist<(*property)->n_ordinate;i_hist++)
         init_histogram(&((*property)->hist[i_hist]),GBP_HISTOGRAM_FIXED,x_min,dx,n_x);
      (*property)->calc_index_function=calc_SSFctn_hist_index;
   }
   else
      SID_trap_error("Invalid property {%s} requested in init_treenode_trend_property().",ERROR_LOGIC,name);
}

void free_trend_property(treenode_trend_property_info **property);
void free_trend_property(treenode_trend_property_info **property){
   if((*property)->hist!=NULL){
      for(int i_hist=0;i_hist<(*property)->n_ordinate;i_hist++)
         free_histogram(&((*property)->hist[i_hist]));
      SID_free(SID_FARG (*property)->hist);
   }
   SID_free(SID_FARG (*property));
}

void init_treenode_trend(tree_info *trees,treenode_trend_info **trend,int mode,const char *filename_out_root,const char *ordinate,int n_properties,...);
void init_treenode_trend(tree_info *trees,treenode_trend_info **trend,int mode,const char *filename_out_root,const char *ordinate,int n_properties,...){
   va_list  vargs;
   va_start(vargs,n_properties);

   (*trend)              =(treenode_trend_info *)SID_malloc(sizeof(treenode_trend_info));
   (*trend)->mode        =mode;
   (*trend)->n_properties=n_properties;
   (*trend)->property    =(treenode_trend_property_info **)SID_malloc(sizeof(treenode_trend_property_info *)*((*trend)->n_properties));
   init_treenode_trend_property(trees,&((*trend)->ordinate),ordinate,1);
   for(int i_property=0;i_property<((*trend)->n_properties);i_property++){
      const char *property=(const char *)va_arg(vargs,const char *);
      init_treenode_trend_property(trees,&((*trend)->property[i_property]),property,((*trend)->ordinate)->hist->n_bins);
   }
   strcpy((*trend)->filename_output_root,filename_out_root);

   va_end(vargs);
}

void free_treenode_trend(treenode_trend_info **trend);
void free_treenode_trend(treenode_trend_info **trend){
   free_trend_property(&((*trend)->ordinate));
   for(int i_property=0;i_property<((*trend)->n_properties);i_property++)
      free_trend_property(&((*trend)->property[i_property]));
   SID_free(SID_FARG (*trend)->property);
   SID_free(SID_FARG (*trend));
}

void write_treenode_trend_binning_file(treenode_trend_info *trend);
void write_treenode_trend_binning_file(treenode_trend_info *trend){
   SID_log("Writing binning description file...",SID_LOG_OPEN);
   if(SID.I_am_Master){
      char  filename_out[MAX_FILENAME_LENGTH];
      sprintf(filename_out,"%s_%s_bins.txt",trend->filename_output_root,trend->ordinate->name);
      FILE *fp_out=fopen(filename_out,"w");
      int i_column=1;
      fprintf(fp_out,"#  Ordinate %s-binning for {%s}\n",trend->ordinate->name,trend->filename_output_root);
      fprintf(fp_out,"#  Column (%03d): Bin index\n",i_column++);
      fprintf(fp_out,"#         (%03d): %s - lo\n",  i_column++,trend->ordinate->name);
      fprintf(fp_out,"#         (%03d): %s - hi\n",  i_column++,trend->ordinate->name);
      hist_info *hist=trend->ordinate->hist;
      for(int i_bin=0;i_bin<hist->n_bins;i_bin++)
         fprintf(fp_out,"%03d %le %le\n",i_bin,histogram_bin_x_lo(hist,i_bin),histogram_bin_x_hi(hist,i_bin));
      fclose(fp_out);
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

void write_treenode_trend_ascii(treenode_trend_info *trend);
void write_treenode_trend_ascii(treenode_trend_info *trend){

   SID_log("Writing trend to output file...",SID_LOG_OPEN);

   // Write binning description file
   write_treenode_trend_binning_file(trend);

   // Open files and write file headers
   SID_log("Opening trend output files and writing headers...",SID_LOG_OPEN);

   // Set filename and open file
   char filename[MAX_FILENAME_LENGTH];
   sprintf(filename,"%s_%s.txt",trend->filename_output_root,trend->ordinate->name);
   FILE *fp_out=fopen(filename,"w");

   // Write header
   int i_column=1;
   fprintf(fp_out,"#  Column (%03d): Snapshot\n",   i_column++);
   fprintf(fp_out,"#         (%03d): %s bin - lo\n",i_column++,trend->ordinate->name);
   fprintf(fp_out,"#         (%03d): %s bin - hi\n",i_column++,trend->ordinate->name);
   fprintf(fp_out,"#         (%03d): n_halos_all\n",i_column++);
   for(int i_property=0;i_property<trend->n_properties;i_property++){
      fprintf(fp_out,"#         (%03d): n_halos_hist (%s)\n",      i_column++,trend->property[i_property]->name);
      fprintf(fp_out,"#         (%03d): %s\n",                     i_column++,trend->property[i_property]->name);
      fprintf(fp_out,"#         (%03d): %s - 68%% confidence lo\n",i_column++,trend->property[i_property]->name);
      fprintf(fp_out,"#         (%03d): %s - 68%% confidence hi\n",i_column++,trend->property[i_property]->name);
      fprintf(fp_out,"#         (%03d): %s - 95%% confidence lo\n",i_column++,trend->property[i_property]->name);
      fprintf(fp_out,"#         (%03d): %s - 95%% confidence hi\n",i_column++,trend->property[i_property]->name);
   }

   // Finalize the snapshot histograms and write results
   hist_info *hist_ordinate=trend->ordinate->hist;
   for(int i_bin=0;i_bin<hist_ordinate->n_bins;i_bin++){
      fprintf(fp_out,"%3d %le %le %d",
                     i_bin,
                     histogram_bin_x_lo(hist_ordinate,i_bin),
                     histogram_bin_x_hi(hist_ordinate,i_bin),
                     hist_ordinate->bin_count[i_bin]);
      for(int i_property=0;i_property<trend->n_properties;i_property++){
         double x_peak;
         double x_68_lo;
         double x_68_hi;
         double x_95_lo;
         double x_95_hi;
         hist_info *hist_i=&(trend->property[i_property]->hist[i_bin]);
         finalize_histogram(hist_i);
         compute_histogram_range(hist_i,68.,GBP_HISTOGRAM_RANGE_HIST,&x_peak,&x_68_lo,&x_68_hi);
         compute_histogram_range(hist_i,95.,GBP_HISTOGRAM_RANGE_HIST,&x_peak,&x_95_lo,&x_95_hi);
         fprintf(fp_out," %d %le %le %le %le %le",
                        hist_i->count_hist,
                        x_peak,
                        x_68_lo,
                        x_68_hi,
                        x_95_lo,
                        x_95_hi);
      }
      fprintf(fp_out,"\n");
   }
   fclose(fp_out);
   SID_log("Done.",SID_LOG_CLOSE);

   SID_log("Done.",SID_LOG_CLOSE);
}

// Structure that will carry the needed information to the select-and-analyze function
typedef struct select_and_analyze_params_local select_and_analyze_params_local;
struct select_and_analyze_params_local{
   char                 filename_output_root[MAX_FILENAME_LENGTH];
   treenode_trend_info *trends_z;
};

// ** Define the calculation here **
void select_and_analyze_treenodes_fctn_init_local(tree_info *trees,void *params_in,int mode,int i_type);
void select_and_analyze_treenodes_fctn_init_local(tree_info *trees,void *params_in,int mode,int i_type){

  // Create an alias for the passed void pointer
  select_and_analyze_params_local *params=(select_and_analyze_params_local *)params_in;

  // Initialize for analyzing groups or subgroups
  char filename_output_root[MAX_FILENAME_LENGTH];
  if(i_type==0){
     sprintf(filename_output_root,"%s_groups",params->filename_output_root);
     precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
     write_treenode_markers_all(trees,params->filename_output_root,PRECOMPUTE_TREENODE_MARKER_GROUPS);
  }
  else{
     sprintf(filename_output_root,"%s_subgroups",params->filename_output_root);
     precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
     write_treenode_markers_all(trees,params->filename_output_root,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
  }

  // Initialize the trend
  init_treenode_trend(trees,&(params->trends_z),INIT_TREENODE_TREND_DEFAULT,filename_output_root,
                      "z",5,"x_off","SSFctn","tau_form","tau_3to1","tau_10to1");
  
}

// ** Perform the calculation here **
void select_and_analyze_treenodes_fctn_analyze_local(tree_info *trees,void *params,int mode,int i_type,int flag_init,tree_node_info *halo);
void select_and_analyze_treenodes_fctn_analyze_local(tree_info *trees,void *params,int mode,int i_type,int flag_init,tree_node_info *halo){
   add_treenode_to_trend(trees,((select_and_analyze_params_local *)params)->trends_z,halo);
}

// ** Write the results here **
void select_and_analyze_treenodes_fctn_fin_local(tree_info *trees,void *params_in,int mode,int i_type);
void select_and_analyze_treenodes_fctn_fin_local(tree_info *trees,void *params_in,int mode,int i_type){
   // Create an alias for the passed void pointer
   select_and_analyze_params_local *params=(select_and_analyze_params_local *)params_in;

   // Write results and clean-up trend
   write_treenode_trend_ascii(params->trends_z);
   free_treenode_trend(&(params->trends_z));

   // Clean-up markers
   if(i_type==0) free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
   else          free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
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

