#include <gbpLib.h>

double calc_max_global(void   *data_local,
 	    	       size_t  n_data_local,
                       SID_Datatype     type){
  int     i_data;
  double  max_local,max;
  max_local=calc_max(data_local,n_data_local,type);
  #ifdef USE_MPI
    MPI_Allreduce(&max_local,&max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  #else
    max=max_local;
  #endif
  return(max);
}
