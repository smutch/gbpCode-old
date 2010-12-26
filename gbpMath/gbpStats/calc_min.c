#include <gbpLib.h>

double calc_min(void   *data,
 		size_t  n_data,
                SID_Datatype type){
  int     i_data;
  double  min;
  if(n_data<1)
    min=0.;
  else{
      if(type==SID_DOUBLE)
        min=(double)((double *)data)[0];
      else if(type==SID_FLOAT)
        min=(double)((float  *)data)[0];
      else if(type==SID_INT)
        min=(double)((int    *)data)[0];
      else if(type==SID_SIZE_T)
        min=(double)((size_t *)data)[0];
      else
        SID_trap_error("Unknown variable type in calc_min",ERROR_LOGIC);
    for(i_data=1;i_data<n_data;i_data++){
        if(type==SID_DOUBLE)
          min=MIN(min,(double)((double *)data)[i_data]);
        else if(type==SID_FLOAT)
          min=MIN(min,(double)((float  *)data)[i_data]);
        else if(type==SID_INT)
          min=MIN(min,(double)((int    *)data)[i_data]);
        else if(type==SID_SIZE_T)
          min=MIN(min,(double)((size_t *)data)[i_data]);
        else
          SID_trap_error("Unknown variable type in calc_min",ERROR_LOGIC);
    }
  }
  return(min);
}
