#include <gbpLib.h>

double calc_sum(void   *data,
 		size_t  n_data,
                SID_Datatype type){
  int     i_data;
  double  sum;
  if(n_data<1)
    sum=0.;
  else{
    for(i_data=0,sum=0.;i_data<n_data;i_data++){
        if(type==SID_DOUBLE)
          sum+=(double)((double *)data)[i_data];
        else if(type==SID_FLOAT)
          sum+=(double)((float  *)data)[i_data];
        else if(type==SID_INT)
          sum+=(double)((int    *)data)[i_data];
        else if(type==SID_SIZE_T)
          sum+=(double)((size_t *)data)[i_data];
        else
          SID_trap_error("Unknown variable type in calc_sum",ERROR_LOGIC);
    }
  }
  return(sum);
}
