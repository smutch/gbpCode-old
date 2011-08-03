#ifndef GBPSTATS_AWAKE
#define GBPSTATS_AWAKE
#include <gbpLib.h>

#define CALC_MODE_DEFAULT       DEFAULT_MODE
#define CALC_MODE_RETURN_DOUBLE 1

void calc_max(void   *data,
              void   *result,
	      size_t  n_data,
              SID_Datatype type,
              int          mode);
void calc_max_global(void   *data_local,
                     void   *result,
 	    	     size_t  n_data_local,
                     SID_Datatype type,
                     int          mode,
                     SID_Comm    *comm);
void calc_mean(void  *data,
               void  *result,
	       size_t n_data,
               SID_Datatype type,
               int          mode);
void calc_mean_global(void   *data_local,
                      void   *result,
 		      size_t  n_data_local,
                      SID_Datatype type,
                      int          mode,
                      SID_Comm    *comm);
void calc_median(void   *data,
                 void   *result,
		 size_t  n_data,
                 SID_Datatype type,
                 int          mode);
void calc_median_global(void   *data_local,
                        void   *result,
		        size_t  n_data_local,
                        SID_Datatype type,
                        int          mode,
                        SID_Comm    *comm);
void calc_min(void   *data,
              void   *result,
 	      size_t  n_data,
              SID_Datatype type,
              int          mode);
void calc_min_global(void   *data_local,
                     void   *result,
 	    	     size_t  n_data_local,
                     SID_Datatype type,
                     int          mode,
                     SID_Comm    *comm);
void calc_stddev(void   *data,
                 void   *result,
 		 size_t  n_data,
                 SID_Datatype type,
                 int          mode);
void calc_stddev_global(void   *data_local,
                        void   *result,
 		        size_t  n_data_local,
                        SID_Datatype type,
                        int          mode,
                        SID_Comm    *comm);
void calc_sum(void   *data,
              void   *result,
 	      size_t  n_data,
              SID_Datatype type,
              int          mode);
void calc_sum_global(void   *data_local,
                     void   *result,
 	    	     size_t  n_data_local,
                     SID_Datatype type,
                     int          mode,
                     SID_Comm    *comm);
double calc_sep_periodic(double x_1,
                         double y_1,
                         double z_1,
                         double x_2,
                         double y_2,
                         double z_2,
                         double box_size);
#endif
