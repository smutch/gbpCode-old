#include <stdio.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>

void init_halo_trend_property_z(trend_property_info *property,void *property_data,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_IRREGULAR_XLO_DEFINED;
   halo_trend_info *halo_trend =(halo_trend_info *)(property->params);
   double    *x_lo  =halo_trend->z_list;
   int        n_bins=halo_trend->n_snaps;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double *),&x_lo);
   gbp_add_va_arg(vargs,sizeof(int),     &n_bins);
}
void free_halo_trend_property_z(trend_property_info *property,void *property_data,int i_hist,int *mode,gbp_va_list *vargs){}
int calc_halo_trend_property_index_z(trend_property_info *property,hist_info *hist,void *halo){
   return(((halo_info *)halo)->snapshot);
}

void init_halo_trend_property_logM_FoF(trend_property_info *property,void *property_data,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_FIXED;
   double x_min = 7.0;
   double dx    = 0.1;
   int    n_x   =  80;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double),&x_min);
   gbp_add_va_arg(vargs,sizeof(double),&dx);
   gbp_add_va_arg(vargs,sizeof(int),   &n_x);
}
void free_halo_trend_property_logM_FoF(trend_property_info *property,void *property_data,int i_hist,int *mode,gbp_va_list *vargs){}
int calc_halo_trend_property_index_logM_FoF(trend_property_info *property,hist_info *hist,void *halo_in){
   halo_trend_info *halo_trend_data=(halo_trend_info *)(property->params);
   halo_info       *halo           =(halo_info       *)(halo_in);
   return(calc_histogram_index(hist,take_log10(halo->properties_group->n_particles*halo_trend_data->m_p)));
}

void init_halo_trend_property_SSFctn(trend_property_info *property,void *property_data,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_FIXED;
   double x_min = 0.0;
   double dx    =0.01;
   int    n_x   = 100;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double),&x_min);
   gbp_add_va_arg(vargs,sizeof(double),&dx);
   gbp_add_va_arg(vargs,sizeof(int),   &n_x);
}
void free_halo_trend_property_SSFctn(trend_property_info *property,void *property_data,int i_hist,int *mode,gbp_va_list *vargs){}
int calc_halo_trend_property_index_SSFctn(trend_property_info *property,hist_info *hist,void *halo_in){
   halo_trend_info *halo_trend_data=(halo_trend_info *)(property->params);
   halo_info       *halo           =(halo_info       *)(halo_in);
   double SSFctn=(double)(halo->np_sub-halo->np_sub_largest)/(double)halo->np_sub;
   return(calc_histogram_index(hist,SSFctn));
}

