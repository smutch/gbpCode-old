#include <math.h>
#include <gbpLib.h>
#include <gbpStats.h>

double calc_stddev_global(void   *data_local,
 		          size_t  n_data_local,
                          SID_Datatype  type){
  size_t  i_data;
  size_t  n_data;
  double  accumulator_local,accumulator;
  double  mean,stddev;
  if(n_data_local<1)
    accumulator_local=0.;
  else{
    mean=calc_mean_global(data_local,n_data_local,type);
    for(i_data=0,accumulator_local=0.;i_data<n_data_local;i_data++){
        if(type==SID_DOUBLE)
          accumulator_local+=pow((double)((double *)data_local)[i_data]-mean,2.);
        else if(type==SID_FLOAT)
          accumulator_local+=pow((double)((float  *)data_local)[i_data]-mean,2.);
        else if (type==SID_INT)
          accumulator_local+=pow((double)((int    *)data_local)[i_data]-mean,2.);
        else if(type==SID_SIZE_T)
          accumulator_local+=pow((double)((size_t *)data_local)[i_data]-mean,2.);
        else
          SID_trap_error("Unknown variable type in calc_stddev",ERROR_LOGIC);
    }
  }
  #if USE_MPI
    MPI_Allreduce(&accumulator_local,&accumulator,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&n_data_local,     &n_data,     1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
  #else
    accumulator=accumulator_local;
    n_data     =n_data_local;
  #endif
  stddev=sqrt(accumulator/(double)n_data);
  return(stddev);
}
