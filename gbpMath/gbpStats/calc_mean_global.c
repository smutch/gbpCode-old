#include <gbpLib.h>
#include <gbpStats.h>

double calc_mean_global(void   *data_local,
 		        size_t  n_data_local,
                        SID_Datatype type){
  double  sum;
  size_t  n_data;
  double  mean;
  if(n_data_local<1)
    mean=0.;
  else{
    sum=calc_sum_global(data_local,n_data_local,type);
    #ifdef USE_MPI
      MPI_Allreduce(&n_data_local,&n_data,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
    #else
      n_data=n_data_local;
    #endif
    mean=sum/(double)n_data;
  }
  return(mean);
}
