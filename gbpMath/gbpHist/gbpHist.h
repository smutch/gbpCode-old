#ifndef GBPHIST_AWAKE
#define GBPHIST_AWAKE
#include <gbpLib.h>

#define GBP_HISTOGRAM_FIXED                 TTTP01
#define GBP_HISTOGRAM_IRREGULAR_XLO_DEFINED TTTP02
#define GBP_HISTOGRAM_IRREGULAR_REVERSED    TTTP03
#define GBP_HISTOGRAM_RANGE_ALL             TTTP04
#define GBP_HISTOGRAM_RANGE_HIST            TTTP05
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
void   finalize_histogram(hist_info *hist);
void   compute_histogram_range(hist_info *hist,double confidence_percent,int mode,double *x_peak,double *x_lo,double *x_hi);
double histogram_bin_x_lo(hist_info *hist,int bin);
double histogram_bin_x_hi(hist_info *hist,int bin);
double histogram_bin_x_mid(hist_info *hist,int bin);

#ifdef __cplusplus
}
#endif
#endif

