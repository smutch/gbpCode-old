#include <gbpLib.h>
#include <gbpSort.h>
#include <gbpStats.h>

void calc_median(void   *data,
                 void   *result,
     	         size_t  n_data,
                 SID_Datatype  type,
                 int           mode){
  int     i_mid;
  double  median;
  size_t *index;

  if(n_data<1){
    if(type==SID_DOUBLE || check_mode_for_flag(mode,CALC_MODE_RETURN_DOUBLE))
      ((double *)result)[0]=0.;
    else if(type==SID_FLOAT)
      ((float  *)result)[0]=0.;
    else if(type==SID_INT)
      ((int    *)result)[0]=0;
    else if(type==SID_SIZE_T)
      ((size_t *)result)[0]=0;
    else
      SID_trap_error("Unknown variable type in calc_median",ERROR_LOGIC);
  }
  else{
    merge_sort(data,
               n_data,
               &index,
               type,
               SORT_COMPUTE_INDEX,
               FALSE);
    i_mid=n_data/2;
    if(n_data%2){
      if(type==SID_DOUBLE)
        median=(double)((double *)data)[index[i_mid]];
      else if(type==SID_FLOAT)
        median=(double)((float  *)data)[index[i_mid]];
      else if(type==SID_INT)
        median=(double)((int    *)data)[index[i_mid]];
      else if(type==SID_SIZE_T)
        median=(double)((size_t *)data)[index[i_mid]];
      else
        SID_trap_error("type not supported in calc_median.",ERROR_LOGIC);
    }
    else{
      if(type==SID_DOUBLE)
        median=ONE_HALF*(double)(((double *)data)[index[i_mid-1]]+
                                 ((double *)data)[index[i_mid]]);
      else if(type==SID_FLOAT)
        median=ONE_HALF*(double)(((float *)data)[index[i_mid-1]]+
                                 ((float *)data)[index[i_mid]]);
      else if(type==SID_INT)
        median=ONE_HALF*(double)(((int *)data)[index[i_mid-1]]+
                                 ((int *)data)[index[i_mid]]);
      else if(type==SID_SIZE_T)
        median=ONE_HALF*(double)(((size_t *)data)[index[i_mid-1]]+
                                 ((size_t *)data)[index[i_mid]]);
      else
        SID_trap_error("type not supported in calc_median.",ERROR_LOGIC);
    }
    SID_free(SID_FARG index);
    if(type==SID_DOUBLE || check_mode_for_flag(mode,CALC_MODE_RETURN_DOUBLE))
      ((double *)result)[0]=(double)median;
    else if(type==SID_FLOAT)
      ((float  *)result)[0]=(float)median;
    else if(type==SID_INT)
      ((int    *)result)[0]=(int)median;
    else if(type==SID_SIZE_T)
      ((size_t *)result)[0]=(size_t)median;
    else
      SID_trap_error("type not supported in calc_median.",ERROR_LOGIC);
  }
}
