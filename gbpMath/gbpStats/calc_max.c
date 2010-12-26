#include <gbpLib.h>

double calc_max(void   *data,
 		size_t  n_data,
                SID_Datatype  type){
  int     i_data;
  double  max;
  if(n_data<1)
    max=0.;
  else{
      if(type==SID_DOUBLE)
        max=(double)((double *)data)[0];
      else if(type==SID_FLOAT)
        max=(double)((float  *)data)[0];
      else if(type==SID_INT)
        max=(double)((int    *)data)[0];
      else if(type==SID_SIZE_T)
        max=(double)((size_t *)data)[0];
      else
        SID_trap_error("Unknown variable type in calc_max\n",ERROR_LOGIC);
    for(i_data=1;i_data<n_data;i_data++){
        if(type==SID_DOUBLE)
          max=MAX(max,(double)((double *)data)[i_data]);
        else if(type==SID_FLOAT)
          max=MAX(max,(double)((float  *)data)[i_data]);
        else if(type==SID_INT)
          max=MAX(max,(double)((int    *)data)[i_data]);
        else if(type==SID_SIZE_T)
          max=MAX(max,(double)((size_t *)data)[i_data]);
        else
          SID_trap_error("Unknown datatype in calc_max\n",ERROR_LOGIC);
    }
  }
  return(max);
}
