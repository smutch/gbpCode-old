#include <gbpLib.h>

double calc_sum_global(void   *data_local,
 	    	       size_t  n_data_local,
                       SID_Datatype  type){
  int     i_data;
  double  sum_local,sum;
  sum_local=calc_sum(data_local,n_data_local,type);
  #if USE_MPI
    MPI_Allreduce(&sum_local,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  #else
    sum=sum_local;
  #endif
  return(sum);
}
