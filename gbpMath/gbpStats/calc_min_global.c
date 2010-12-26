#include <gbpLib.h>
#include <gbpStats.h>

double calc_min_global(void   *data_local,
 	    	       size_t  n_data_local,
                       SID_Datatype   type){
  int     i_data;
  double  min_local,min;
  min_local=calc_min(data_local,n_data_local,type);
  #ifdef USE_MPI
    MPI_Allreduce(&min_local,&min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  #else
    min=min_local;
  #endif
  return(min);
}
