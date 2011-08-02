#ifndef GBPSTATS_AWAKE
#define GBPSTATS_AWAKE
#include <gbpLib.h>

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
<<<<<<< HEAD
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
=======
>>>>>>> 4e5f3e7f51c761d5faade46f5a9ecca3fa094716
#endif
