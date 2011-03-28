#include <gbpLib.h>
#include <gbpStats.h>

double calc_median_global(void   *data_local,
		          size_t  n_data_local,
                          SID_Datatype type){
  int     i_mid;
  double  median;
  size_t *index;
  double  value_1,value_2;
  size_t  n_data;
  size_t  index_1,index_2;
  size_t *idx;
  size_t  i_lo,i_hi;
  size_t  index_lo,index_hi,index_mid;
  size_t  i_1,i_2;
  int     flag_1,flag_2;
  int     rank_1_local,rank_1,rank_2_local,rank_2;

  #if USE_MPI
  SID_trap_error("calc_median_global has not been debugged yet!",ERROR_LOGIC);
  MPI_Allreduce(&n_data_local,&n_data,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
  if(n_data>0){
    if(n_data>1){
      sort(data_local,
           n_data_local,
           &index,
           type,
           FALSE,
           TRUE,
           FALSE);
      index_1=n_data/2-1;
      index_2=n_data/2;
      sort(index,
           n_data_local,
           &idx,
           SID_SIZE_T,
           TRUE,
           TRUE,
           FALSE);
      i_lo    =0;
      i_hi    =n_data_local-1;
      index_lo=index[idx[i_lo]];
      index_hi=index[idx[i_hi]];
      do{     
        i_mid    =(i_lo+i_hi)/2;
        index_mid=index[idx[i_mid]];
        if(index_1<index_mid)
          i_hi=i_mid;
        else if(index_1>index_mid)
          i_lo=i_mid;
        else{ 
          i_lo=i_mid;
          i_hi=i_mid;
        }     
      } while(i_hi-i_lo>2);
      if(index_1==index[i_lo]){
        i_1   =i_lo;
        flag_1=TRUE;
      }
      else
        flag_1=FALSE;
      if(index_2==index[i_lo]){
        i_2   =i_lo;
        flag_2=TRUE;
      }
      else if(index_2==index[i_hi]){
        i_2   =i_hi;
        flag_2=TRUE;
      }
      else
        flag_2=FALSE;
      if(!(n_data%2)){
        if(flag_1){
          rank_1_local=SID.My_rank;
            if(type==SID_DOUBLE)
              value_1=((double *)data_local)[i_1];
            else if(type==SID_FLOAT)
              value_1=((float  *)data_local)[i_1];
            else if(type==SID_INT)
              value_1=((int    *)data_local)[i_1];
            else if(type==SID_SIZE_T)
              value_1=((size_t *)data_local)[i_1];
            else
              SID_trap_error("type not supported by calc_median_global.",ERROR_LOGIC);
          }
        }
        else
          rank_1_local=-1;
        MPI_Allreduce(&rank_1_local,&rank_1,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
        MPI_Bcast(&value_1,1,MPI_DOUBLE,rank_1,MPI_COMM_WORLD);
      }
      if(flag_2){
        rank_2_local=SID.My_rank;
          if(type==SID_DOUBLE)
            value_2=((double *)data_local)[i_2];
          else if(type==SID_FLOAT)
            value_2=((float  *)data_local)[i_2];
          else if(type==SID_INT)
            value_2=((int    *)data_local)[i_2];
          else if(type==SID_SIZE_T)
            value_2=((size_t *)data_local)[i_2];
          else
            SID_trap_error("type not supported by calc_median_global.",ERROR_LOGIC);
      }
      else
        rank_2_local=-1;
      MPI_Allreduce(&rank_2_local,&rank_2,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
      MPI_Bcast(&value_2,1,MPI_DOUBLE,rank_2,MPI_COMM_WORLD);          
      free(idx);
      free(index);

      if(n_data%2)
        median=value_2;
      else
        median=ONE_HALF*(value_1+value_2);
    }
    else
        if(type==SID_DOUBLE)
          median=((double *)data_local)[0];
        else if(type==SID_FLOAT)
          median=((float  *)data_local)[0];
        else if(type==SID_INT)
          median=((int    *)data_local)[0];
        else if(type==SID_SIZE_T)
          median=((size_t *)data_local)[0];
        else
          SID_trap_error("type not supported by calc_median_global.",ERROR_LOGIC);
  }
  else
    median=0.;
  #else
    median=calc_median(data_local,n_data_local,type);
  #endif
  return(median);
}
