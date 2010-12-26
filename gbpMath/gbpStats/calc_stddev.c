#include <math.h>
#include <gbpLib.h>
#include <gbpStats.h>

double calc_stddev(void   *data,
 		   size_t  n_data,
                   SID_Datatype type){
  int     i_data;
  double  stddev;
  double  mean;
  if(n_data<1)
    stddev=0.;
  else{
    mean=calc_mean(data,n_data,type);
    for(i_data=0,stddev=0.;i_data<n_data;i_data++){
        if(type==SID_DOUBLE)
          stddev+=pow((double)((double *)data)[i_data]-mean,2.);
        else if(type==SID_FLOAT)
          stddev+=pow((double)((float  *)data)[i_data]-mean,2.);
        else if(type==SID_INT)
          stddev+=pow((double)((int    *)data)[i_data]-mean,2.);
        else if(type==SID_SIZE_T)
          stddev+=pow((double)((size_t *)data)[i_data]-mean,2.);
        else
          SID_trap_error("Unknown variable type in calc_stddev",ERROR_LOGIC);
    }
    stddev=sqrt(stddev/(double)n_data);
  }
  return(stddev);
}
