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

void init_tree_property_logM_course(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_FIXED;
   double x_min = 7.0;
   double dx    = 0.5;
   int    n_x   =  18;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double),&x_min);
   gbp_add_va_arg(vargs,sizeof(double),&dx);
   gbp_add_va_arg(vargs,sizeof(int),   &n_x);
}
void init_tree_property_logM(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_FIXED;
   double x_min = 7.0;
   double dx    = 0.1;
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
   int    n_x   = 100;
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
   double dx    =0.01;
   int    n_x   = 100;
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

void init_tree_property_Vir_ratio(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_FIXED;
   double x_min = 0.0;
   double dx    =0.02;
   int    n_x   = 200;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double),&x_min);
   gbp_add_va_arg(vargs,sizeof(double),&dx);
   gbp_add_va_arg(vargs,sizeof(int),   &n_x);
}
void free_tree_property_Vir_ratio(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){}
int calc_tree_property_index_Vir_ratio(trend_property_info *property,hist_info *hist,void *halo_in){
   tree_info      *trees=(tree_info      *)(property->params);
   tree_node_info *halo =(tree_node_info *)(halo_in);
   return(calc_histogram_index(hist,fetch_treenode_Vir_ratio(trees,halo)));
}

void init_tree_property_log_sigma_vx(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_FIXED;
   double x_min = 0.0;
   double dx    =0.01;
   int    n_x   = 500;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double),&x_min);
   gbp_add_va_arg(vargs,sizeof(double),&dx);
   gbp_add_va_arg(vargs,sizeof(int),   &n_x);
}
void free_tree_property_log_sigma_vx(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){}
int calc_tree_property_index_log_sigma_vx(trend_property_info *property,hist_info *hist,void *list_in){
   tree_info          *trees=(tree_info          *)(property->params);
   treenode_list_info *list =(treenode_list_info *)(list_in);
   return(calc_histogram_index(hist,fetch_treenode_list_local_log_sigma_vx(trees,list)));
}

void init_tree_property_tau(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_IRREGULAR_XLO_DEFINED;
   tree_info *trees    =(tree_info *)trees_in;
   double    *tau_array=(double *)SID_malloc(sizeof(double)*trees->n_snaps);
   int        n_bins; 
   int        i_0;
   if(property->is_ordinate){
      n_bins=trees->n_snaps;
      i_0   =trees->n_snaps-1;
   }
   else{
      n_bins=MIN(i_hist+1,trees->n_snaps);
      i_0   =i_hist;
   }
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double),&tau_array);
   gbp_add_va_arg(vargs,sizeof(int),   &n_bins);
   for(int i_tau=0;i_tau<n_bins;i_tau++)
      tau_array[i_tau]=10.*((trees->t_list[i_0]-trees->t_list[i_0-i_tau])/trees->t_list[i_0]);
}
void free_tree_property_tau(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){
   double *tau_array;
   gbp_va_start(vargs);
   gbp_fetch_va_arg(vargs,sizeof(double),&tau_array);
   SID_free(SID_FARG tau_array);
}
int calc_tree_property_index_tau_form(trend_property_info *property,hist_info *hist,void *halo_in){
   int             r_val=-1;
   tree_node_info *halo =(tree_node_info *)(halo_in);
   if(halo!=NULL){
      tree_info         *trees  =(tree_info *)(property->params);
      tree_markers_info *markers=fetch_treenode_precomputed_markers(trees,halo);
      tree_node_info    *marker =markers->half_peak_mass;
      if(marker!=NULL){
         int    halo_snap  =fetch_treenode_snap_tree(trees,halo);
         int    marker_snap=fetch_treenode_snap_tree(trees,marker);
         double tau        =10.*((trees->t_list[halo_snap]-trees->t_list[marker_snap])/trees->t_list[halo_snap]);
         r_val=calc_histogram_index(hist,tau);
      }
   }
   return(r_val);
}
int calc_tree_property_index_tau_3to1(trend_property_info *property,hist_info *hist,void *halo_in){
   int             r_val=-1;
   tree_node_info *halo =(tree_node_info *)(halo_in);
   if(halo!=NULL){
      tree_info         *trees       =(tree_info *)(property->params);
      tree_markers_info *markers_halo=fetch_treenode_precomputed_markers(trees,halo);
      tree_node_info    *merger      =markers_halo->merger_33pc_remnant;
      if(merger!=NULL){
         tree_markers_info *markers_merger=fetch_treenode_precomputed_markers(trees,merger);
         tree_node_info    *marker        =markers_merger->peak_mass;
         if(marker!=NULL){
            int    halo_snap  =fetch_treenode_snap_tree(trees,halo);
            int    marker_snap=fetch_treenode_snap_tree(trees,marker);
            double tau        =10.*((trees->t_list[halo_snap]-trees->t_list[marker_snap])/trees->t_list[halo_snap]);
            r_val=calc_histogram_index(hist,tau);
         }
      }
   }
   return(r_val);
}
int calc_tree_property_index_tau_10to1(trend_property_info *property,hist_info *hist,void *halo_in){
   int             r_val=-1;
   tree_node_info *halo =(tree_node_info *)(halo_in);
   if(halo!=NULL){
      tree_info         *trees       =(tree_info *)(property->params);
      tree_markers_info *markers_halo=fetch_treenode_precomputed_markers(trees,halo);
      tree_node_info    *merger      =markers_halo->merger_10pc_remnant;
      if(merger!=NULL){
         tree_markers_info *markers_merger=fetch_treenode_precomputed_markers(trees,merger);
         tree_node_info    *marker        =markers_merger->peak_mass;
         if(marker!=NULL){
            int    halo_snap  =fetch_treenode_snap_tree(trees,halo);
            int    marker_snap=fetch_treenode_snap_tree(trees,marker);
            double tau        =10.*((trees->t_list[halo_snap]-trees->t_list[marker_snap])/trees->t_list[halo_snap]);
            r_val=calc_histogram_index(hist,tau);
         }
      }
   }
   return(r_val);
}

