#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHist.h>

void init_histogram(hist_info *hist,int mode,...){
  va_list  vargs;
  va_start(vargs,mode);

  // Set some initial defaults
  hist->flag_bin_order_inverted=FALSE;

  // Set things dependant on mode
  if(check_mode_for_flag(mode,GBP_HISTOGRAM_FIXED)){
     hist->x_min =(double)va_arg(vargs,double);
     hist->dx    =(double)va_arg(vargs,double);
     hist->n_bins=(int   )va_arg(vargs,int);
     hist->x_lo  =NULL;
     hist->x_hi  =NULL;
  }
  else if(check_mode_for_flag(mode,GBP_HISTOGRAM_IRREGULAR_XLO_DEFINED)){
     double *x_lo_in;
     x_lo_in     =(double *)va_arg(vargs,double *);
     hist->n_bins=(int     )va_arg(vargs,int);
     hist->x_lo  =(double *)SID_malloc(sizeof(double)*hist->n_bins);
     hist->x_hi  =NULL;
     memcpy(hist->x_lo,x_lo_in,sizeof(double)*hist->n_bins);
     // Check that the x_min values are sorted properly
     int flag_error=FALSE;
     for(int i_bin=0;i_bin<(hist->n_bins-1);i_bin++){
        if(hist->x_lo[i_bin+1]<hist->x_lo[i_bin])
           flag_error=TRUE;
     }
     // Try reversing axis order if they are not
     if(flag_error){
        flag_error=FALSE;
        for(int i_bin=0;i_bin<(hist->n_bins-1);i_bin++){
           if(hist->x_lo[i_bin]<hist->x_lo[i_bin+1])
              flag_error=TRUE;
        }
        if(!flag_error)
           hist->flag_bin_order_inverted=TRUE;
     }
     // If reversing bin order failed, then this array is not sorted
     if(flag_error)
        SID_trap_error("The x_lo array passed to init_histogram() is not sorted in ascending order",ERROR_LOGIC);
  }
  else
     SID_trap_error("Invalid mode (%d) specified in init_histogram().",ERROR_LOGIC,mode);

  // Set things common to all modes
  hist->mode      =mode;
  hist->count_hist=0;
  hist->count_all =0;
  hist->bin_count =(size_t *)SID_calloc(sizeof(size_t)*hist->n_bins);

  va_end(vargs);
}

