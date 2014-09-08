#include <stdio.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

void init_tree_property_z(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_IRREGULAR_XLO_DEFINED;
   tree_info *trees =(tree_info *)(property->params);
   double    *x_lo  =trees->z_list;
   int        n_bins=trees->n_snaps;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double *),&x_lo);
   gbp_add_va_arg(vargs,sizeof(int),     &n_bins);
}
void free_tree_property_z(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){}
int calc_tree_property_index_z(trend_property_info *property,hist_info *hist,void *halo){
   return(((tree_node_info *)halo)->snap_tree);
}

void init_tree_property_logM(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_FIXED;
   double x_min = 7.0;
   double dx    = 0.5;
   int    n_x   =  14;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double),&x_min);
   gbp_add_va_arg(vargs,sizeof(double),&dx);
   gbp_add_va_arg(vargs,sizeof(int),   &n_x);
}
void free_tree_property_logM(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){}
int calc_tree_property_index_logM(trend_property_info *property,hist_info *hist,void *halo_in){
   tree_info      *trees=(tree_info      *)(property->params);
   tree_node_info *halo =(tree_node_info *)(halo_in);
   return(calc_histogram_index(hist,take_log10(fetch_treenode_Mvir(trees,halo))));
}

void init_tree_property_xoff(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_FIXED;
   double x_min = 0.0;
   double dx    =0.01;
   int    n_x   =  50;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double),&x_min);
   gbp_add_va_arg(vargs,sizeof(double),&dx);
   gbp_add_va_arg(vargs,sizeof(int),   &n_x);
}
void free_tree_property_xoff(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){}
int calc_tree_property_index_xoff(trend_property_info *property,hist_info *hist,void *halo_in){
   tree_info      *trees=(tree_info      *)(property->params);
   tree_node_info *halo =(tree_node_info *)(halo_in);
   return(calc_histogram_index(hist,fetch_treenode_x_off(trees,halo)));
}

void init_tree_property_SSFctn(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_FIXED;
   double x_min = 0.0;
   double dx    =0.05;
   int    n_x   =  20;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double),&x_min);
   gbp_add_va_arg(vargs,sizeof(double),&dx);
   gbp_add_va_arg(vargs,sizeof(int),   &n_x);
}
void free_tree_property_SSFctn(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){}
int calc_tree_property_index_SSFctn(trend_property_info *property,hist_info *hist,void *halo_in){
   tree_info      *trees=(tree_info      *)(property->params);
   tree_node_info *halo =(tree_node_info *)(halo_in);
   return(calc_histogram_index(hist,fetch_treenode_SSFctn(trees,halo)));
}

void init_tree_property_tau(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_IRREGULAR_XLO_DEFINED;
   tree_info *trees    =(tree_info *)trees_in;
   double    *tau_array=(double *)SID_malloc(sizeof(double)*trees->n_snaps);
   int        n_bins   =trees->n_snaps;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double),&tau_array);
   gbp_add_va_arg(vargs,sizeof(int),   &n_bins);
   for(int i_tau=0;i_tau<trees->n_snaps;i_tau++)
      tau_array[i_tau]=10.*((trees->t_list[i_hist]-trees->t_list[i_tau])/trees->t_list[i_hist]);
}
void free_tree_property_tau(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){
   double *tau_array;
   gbp_va_start(vargs);
   gbp_fetch_va_arg(vargs,sizeof(double),&tau_array);
   SID_free(SID_FARG tau_array);
}
int calc_tree_property_index_tau_form(trend_property_info *property,hist_info *hist,void *halo_in){
   tree_info         *trees=(tree_info      *)(property->params);
   tree_node_info    *halo =(tree_node_info *)(halo_in);
   tree_markers_info *markers=fetch_treenode_precomputed_markers(trees,halo);
   return(fetch_treenode_snap_tree(trees,markers->half_peak_mass));
}
int calc_tree_property_index_tau_3to1(trend_property_info *property,hist_info *hist,void *halo_in){
   tree_info         *trees=(tree_info      *)(property->params);
   tree_node_info    *halo =(tree_node_info *)(halo_in);
   tree_markers_info *markers=fetch_treenode_precomputed_markers(trees,halo);
   return(fetch_treenode_snap_tree(trees,markers->merger_33pc_remnant));
}
int calc_tree_property_index_tau_10to1(trend_property_info *property,hist_info *hist,void *halo_in){
   tree_info         *trees=(tree_info      *)(property->params);
   tree_node_info    *halo =(tree_node_info *)(halo_in);
   tree_markers_info *markers=fetch_treenode_precomputed_markers(trees,halo);
   return(fetch_treenode_snap_tree(trees,markers->merger_10pc_remnant));
}

