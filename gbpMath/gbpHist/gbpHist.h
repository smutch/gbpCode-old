#ifndef GBPHIST_AWAKE
#define GBPHIST_AWAKE
#include <gbpLib.h>

#define GBP_ADD_ITEM_TO_TREND_ORDINATE   TTTP01
#define GBP_ADD_ITEM_TO_TREND_COORDINATE TTTP02
#define GBP_ADD_ITEM_TO_TREND_DEFAULT    GBP_ADD_ITEM_TO_TREND_COORDINATE

#define GBP_HISTOGRAM_GBP_VARGS             TTTP01
#define GBP_HISTOGRAM_FIXED                 TTTP02
#define GBP_HISTOGRAM_IRREGULAR_XLO_DEFINED TTTP03
#define GBP_HISTOGRAM_IRREGULAR_REVERSED    TTTP04
#define GBP_HISTOGRAM_RANGE_ALL             TTTP05
#define GBP_HISTOGRAM_RANGE_HIST            TTTP06
#define GBP_HISTOGRAM_DEFAULT               GBP_HISTOGRAM_FIXED|GBP_HISTOGRAM_RANGE_HIST
typedef struct hist_info hist_info;
struct hist_info{
   int     mode;
   int     n_bins;
   double  x_min;
   double  dx;
   double *x_lo;
   double *x_hi;
   int     flag_bin_order_inverted;
   size_t *bin_count;
   size_t  count_hist;
   size_t  count_all;
   int     flag_finalized; // Needed to prevent from finalizing more than once.
};

typedef struct trend_property_info trend_property_info;
struct trend_property_info{
   char   name[128];
   void  *params;
   void (*init_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp);
   void (*free_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp);
   int  (*calc_function)(trend_property_info *property,hist_info *hist,void *params_calc);
   int    n_hist;
   int    is_ordinate;
   hist_info           *hist;
   trend_property_info *next;
};

typedef struct trend_info trend_info;
struct trend_info{
   trend_property_info *ordinate;
   trend_property_info *coordinate_first;
   trend_property_info *coordinate_last;
};

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif

void   init_histogram(hist_info *hist,int mode,...);
void   free_histogram(hist_info *hist);
void   clear_histogram(hist_info *hist);
int    calc_histogram_index(hist_info *hist,double x);
int    is_histogram_index_in_range(hist_info *hist,int index);
void   add_to_histogram(hist_info *hist,double x);
void   add_to_histogram_index(hist_info *hist,int index);
void   finalize_histogram     (hist_info *hist);
void   finalize_trend_property(trend_property_info *property);
void   finalize_trend         (trend_info          *trend);
void   compute_histogram_range(hist_info *hist,double confidence_percent,int mode,double *x_peak,double *x_lo,double *x_hi);
double histogram_bin_x_lo(hist_info *hist,int bin);
double histogram_bin_x_hi(hist_info *hist,int bin);
double histogram_bin_x_mid(hist_info *hist,int bin);

void init_trend(trend_info **trend,
                const char  *name,
                void        *params,
                void (*init_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                void (*free_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                int  (*calc_function)(trend_property_info *property,hist_info *hist,void *params_calc));
void free_trend(trend_info **trend);
void init_trend_property(trend_property_info **property,
                         const char           *name,
                         int                   flag_is_ordinate,
                         int                   n_hist,
                         void                 *params,
                         void (*init_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                         void (*free_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                         int  (*calc_function)(trend_property_info *property,hist_info *hist,void *params_calc));
void init_trend_ordinate(trend_info  *trend,
                         const char  *name,
                         void        *params, 
                         void (*init_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                         void (*free_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                         int  (*calc_function)(trend_property_info *property,hist_info *hist,void *params_calc));
void init_trend_coordinate(trend_info  *trend,
                           const char  *name,
                           void        *params, 
                           void (*init_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                           void (*free_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                           int  (*calc_function)(trend_property_info *property,hist_info *hist,void *params_calc));
void free_trend_property(trend_property_info **property);
int  add_item_to_trend(trend_info *trend,int mode,void *item_ordinate,...);
void write_trend_ascii(trend_info *trend,const char *filename_root);
void write_trend_property_binning_file(trend_property_info *property,const char *filename_output_root);

#ifdef __cplusplus
}
#endif
#endif

