#include <gbpLib.h>
#include <gbpStats.h>

double calc_median(void   *data,
		   size_t  n_data,
                   SID_Datatype     type){
  int     i_mid;
  double  median;
  size_t *index;
  int     flag_local=TRUE;

  if(n_data>0){

    sort(data,
         n_data,
         &index,
         type,
         flag_local,
         SORT_COMPUTE_INDEX,
         FALSE);

    i_mid=n_data/2;
    if(n_data%2){
        if(type==SID_DOUBLE)
          median=((double *)data)[index[i_mid]];
        else if(type==SID_FLOAT)
          median=((float *)data)[index[i_mid]];
        else if(type==SID_INT)
          median=((int *)data)[index[i_mid]];
        else if(type==SID_SIZE_T)
          median=((int *)data)[index[i_mid]];
        else
          SID_trap_error("type not supported in calc_median.",ERROR_LOGIC);
    }
    else{
        if(type==SID_DOUBLE)
          median=ONE_HALF*(((double *)data)[index[i_mid-1]]+
                           ((double *)data)[index[i_mid]]);
        else if(type==SID_FLOAT)
          median=ONE_HALF*(((float *)data)[index[i_mid-1]]+
                           ((float *)data)[index[i_mid]]);
        else if(type==SID_INT)
          median=ONE_HALF*(((int *)data)[index[i_mid-1]]+
                           ((int *)data)[index[i_mid]]);
        else if(type==SID_SIZE_T)
          median=ONE_HALF*(((size_t *)data)[index[i_mid-1]]+
                           ((size_t *)data)[index[i_mid]]);
        else
          SID_trap_error("type not supported in calc_median.",ERROR_LOGIC);
    }
    free(index);
  }
  else
    median=0.;
  return(median);
}
