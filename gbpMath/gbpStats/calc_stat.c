#include <common.h>
#include <ADaM.h>

size_t calc_stat(void     *data_local,
                 int      *mask_local,
                 size_t    n_data_local,
                 int       type,
                 mode_int  mode,
                 void     *return_value){
#ifdef USE_MPI
  size_t  *index;
  size_t  *idx;
  double   value_1_d,value_2_d;
  size_t   value_1_t,value_2_t;
  size_t   index_1,index_2;
  size_t   i_lo,i_hi;
  size_t   index_lo,index_hi,index_mid;
  size_t   i_1,i_2;
  int      flag_1,flag_2;
  int      rank_1_local,rank_1,rank_2_local,rank_2;
#endif
  size_t   i_mid;
  double   mean;
  double   result_d      =0.;
  double   result_local_d=0.;
  size_t   result_t      =0;
  size_t   result_local_t=0;
  size_t   i_data;
  size_t   n_data;
  size_t   n_data_used_local;
  size_t   n_data_used;
  size_t  *sort_index_local;
  mode_int operation;
  int      flag_local;
  int      flag_init;
  // Determine the desired operation
  if(check_mode_for_flag(mode,CALC_STAT_SUM))
    operation=CALC_STAT_SUM;
  else if(check_mode_for_flag(mode,CALC_STAT_MEAN))
    operation=CALC_STAT_MEAN;
  else if(check_mode_for_flag(mode,CALC_STAT_MIN))
    operation=CALC_STAT_MIN;
  else if(check_mode_for_flag(mode,CALC_STAT_MAX))
    operation=CALC_STAT_MAX;
  else if(check_mode_for_flag(mode,CALC_STAT_STDDEV))
    operation=CALC_STAT_STDDEV;
  else if(check_mode_for_flag(mode,CALC_STAT_MEDIAN))
    operation=CALC_STAT_MEDIAN;
  else
    SID_trap_error("Operation not specified correctly in calc_stat. (mode=%lld)",ERROR_LOGIC,mode);

  // Perform operation
  switch(operation){
    // --------------------------------
    // Compute the **MEAN** of the data
    // --------------------------------
  case CALC_STAT_MEAN:
    n_data_used=calc_stat(data_local,
                          mask_local,
                          n_data_local,
                          type,
                          CALC_STAT_SUM|CALC_STAT_RETURN_DOUBLE,
                          &result_d);
    result_d/=(double)n_data_used;
    result_t =(size_t)result_d;
    break;

    // -------------------------
    // Compute a sum of the data
    // -------------------------
  case CALC_STAT_SUM:
    switch(type){
    case ADaM_DOUBLE:
      result_local_d=0.;
      if(mask_local!=NULL){
        for(i_data=0,n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_d+=(double)((double *)data_local)[i_data];
            n_data_used_local++;
          }
        }
      }
      else{
        for(i_data=0;i_data<n_data_local;i_data++)
          result_local_d+=(double)((double *)data_local)[i_data];
        n_data_used_local=n_data_local;
      }
      result_d   =result_local_d;
      n_data_used=n_data_used_local;
#ifdef USE_MPI
      if(check_mode_for_flag(mode,CALC_STAT_GLOBAL)){
        MPI_Allreduce(&result_local_d,   &result_d,   1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&n_data_used_local,&n_data_used,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      }
#endif
      break;
    case ADaM_FLOAT:
      result_local_d=0.;
      if(mask_local!=NULL){
        for(i_data=0,n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_d+=(double)((float *)data_local)[i_data];
            n_data_used_local++;
          }
        }
      }
      else{
        for(i_data=0;i_data<n_data_local;i_data++)
          result_local_d+=(double)((float *)data_local)[i_data];
        n_data_used_local=n_data_local;
      }
      result_d   =result_local_d;
      n_data_used=n_data_used_local;
#ifdef USE_MPI
      if(check_mode_for_flag(mode,CALC_STAT_GLOBAL)){
        MPI_Allreduce(&result_local_d,   &result_d,   1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&n_data_used_local,&n_data_used,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      }
#endif
      break;
    case ADaM_INT:
      result_local_t=0;
      if(mask_local!=NULL){
        for(i_data=0,n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_t+=(size_t)((int *)data_local)[i_data];
            n_data_used_local++;
          }
        }
      }
      else{
        for(i_data=0;i_data<n_data_local;i_data++)
          result_local_t+=(size_t)((int *)data_local)[i_data];
        n_data_used_local=n_data_local;
      }
      result_t   =result_local_t;
      n_data_used=n_data_used_local;
#ifdef USE_MPI
      if(check_mode_for_flag(mode,CALC_STAT_GLOBAL)){
        MPI_Allreduce(&result_local_t,   &result_t,   1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&n_data_used_local,&n_data_used,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      }
#endif
      result_d=(double)result_t;
      break;
    case ADaM_SIZE_T:
      result_local_t=0;
      if(mask_local!=NULL){
        for(i_data=0,n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_t+=(size_t)((size_t *)data_local)[i_data];
            n_data_used_local++;
          }
        }
      }
      else{
        for(i_data=0;i_data<n_data_local;i_data++)
          result_local_t+=(size_t)((size_t *)data_local)[i_data];
        n_data_used_local=n_data_local;
      }
      result_t   =result_local_t;
      n_data_used=n_data_used_local;
#ifdef USE_MPI
      if(check_mode_for_flag(mode,CALC_STAT_GLOBAL)){
        MPI_Allreduce(&result_local_t,   &result_t,   1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&n_data_used_local,&n_data_used,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      }
#endif
      result_d=(double)result_t;
      break;      
    default:
      SID_trap_error("Unknown variable type {%d} in calc_stat\n",ERROR_LOGIC,type);
      break;
    }
    break;
    
    // -----------------------------------
    // Compute the **MINIMUM** of the data
    // -----------------------------------
  case CALC_STAT_MIN:
    switch(type){
    case ADaM_DOUBLE:
      if(mask_local!=NULL){
        for(i_data=0,flag_init=TRUE;flag_init && i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_d=((double *)data_local)[i_data];
            flag_init     =FALSE;
          }
        }
        for(n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_d=MIN(result_local_d,((double *)data_local)[i_data]);
            n_data_used_local++;
          }
        }
      }
      else{
        result_local_d=((double *)data_local)[0];
        for(i_data=1;i_data<n_data_local;i_data++)
          result_local_d=MIN(result_local_d,((double *)data_local)[i_data]);
        n_data_used_local=n_data_local;
      }
      result_d   =result_local_d;
      n_data_used=n_data_used_local;
#ifdef USE_MPI
      // STILL NEED TO CHECK HERE THAT ALL PROCS HAVE N_DATA_LOCAL>0 ... INCORRECT RESULTS OTHERWISE!
      if(check_mode_for_flag(mode,CALC_STAT_GLOBAL)){
        MPI_Allreduce(&result_local_d,   &result_d,   1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(&n_data_used_local,&n_data_used,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      }
#endif
      break;
    case ADaM_FLOAT:
      if(mask_local!=NULL){
        for(i_data=0,flag_init=TRUE;flag_init && i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_d=((float *)data_local)[i_data];
            flag_init     =FALSE;
          }
        }
        for(n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_d=MIN(result_local_d,((float *)data_local)[i_data]);
            n_data_used_local++;
          }
        }
      }
      else{
        result_local_d=((float *)data_local)[0];
        for(i_data=1;i_data<n_data_local;i_data++)
          result_local_d=MIN(result_local_d,((float *)data_local)[i_data]);
        n_data_used_local=n_data_local;
      }
      result_d   =result_local_d;
      n_data_used=n_data_used_local;
#ifdef USE_MPI
      // STILL NEED TO CHECK HERE THAT ALL PROCS HAVE N_DATA_LOCAL>0 ... INCORRECT RESULTS OTHERWISE!
      if(check_mode_for_flag(mode,CALC_STAT_GLOBAL)){
        MPI_Allreduce(&result_local_d,   &result_d,   1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(&n_data_used_local,&n_data_used,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      }
#endif
      break;
    case ADaM_INT:
      if(mask_local!=NULL){
        for(i_data=0,flag_init=TRUE;flag_init && i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_t=((int *)data_local)[i_data];
            flag_init     =FALSE;
          }
        }
        for(n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_t=MIN(result_local_t,((int *)data_local)[i_data]);
            n_data_used_local++;
          }
        }
      }
      else{
        result_local_t=((int *)data_local)[0];
        for(i_data=1;i_data<n_data_local;i_data++)
          result_local_t=MIN(result_local_t,((int *)data_local)[i_data]);
        n_data_used_local=n_data_local;
      }
      result_t   =result_local_t;
      n_data_used=n_data_used_local;
#ifdef USE_MPI
      // STILL NEED TO CHECK HERE THAT ALL PROCS HAVE N_DATA_LOCAL>0 ... INCORRECT RESULTS OTHERWISE!
      if(check_mode_for_flag(mode,CALC_STAT_GLOBAL)){
        MPI_Allreduce(&result_local_t,   &result_t,   1,MPI_SIZE_T,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(&n_data_used_local,&n_data_used,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      }
#endif
      result_d=(double)result_t;
      break;
    case ADaM_SIZE_T:
      if(mask_local!=NULL){
        for(i_data=0,flag_init=TRUE;flag_init && i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_t=((size_t *)data_local)[i_data];
            flag_init     =FALSE;
          }
        }
        for(n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_t=MIN(result_local_t,((size_t *)data_local)[i_data]);
            n_data_used_local++;
          }
        }
      }
      else{
        result_local_t=((size_t *)data_local)[0];
        for(i_data=1;i_data<n_data_local;i_data++)
          result_local_t=MIN(result_local_t,((size_t *)data_local)[i_data]);
        n_data_used_local=n_data_local;
      }
      result_t   =result_local_t;
      n_data_used=n_data_used_local;
#ifdef USE_MPI
      // STILL NEED TO CHECK HERE THAT ALL PROCS HAVE N_DATA_LOCAL>0 ... INCORRECT RESULTS OTHERWISE!
      if(check_mode_for_flag(mode,CALC_STAT_GLOBAL)){
        MPI_Allreduce(&result_local_t,   &result_t,   1,MPI_SIZE_T,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(&n_data_used_local,&n_data_used,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      }
#endif
      result_d=(double)result_t;
      break;      
    default:
      SID_trap_error("Unknown variable type {%d} in calc_stat\n",ERROR_LOGIC,type);
      break;
    }
    break;
    
    // -----------------------------------
    // Compute the **MAXIMUM** of the data
    // -----------------------------------
  case CALC_STAT_MAX:
    switch(type){
    case ADaM_DOUBLE:
      if(mask_local!=NULL){
        for(i_data=0,flag_init=TRUE;flag_init && i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_d=((double *)data_local)[i_data];
            flag_init     =FALSE;
          }
        }
        for(n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_d=MAX(result_local_d,((double *)data_local)[i_data]);
            n_data_used_local++;
          }
        }
      }
      else{
        result_local_d=((double *)data_local)[0];
        for(i_data=1;i_data<n_data_local;i_data++)
          result_local_d=MAX(result_local_d,((double *)data_local)[i_data]);
        n_data_used_local=n_data_local;
      }
      result_d   =result_local_d;
      n_data_used=n_data_used_local;
#ifdef USE_MPI
      // STILL NEED TO CHECK HERE THAT ALL PROCS HAVE N_DATA_LOCAL>0 ... INCORRECT RESULTS OTHERWISE!
      if(check_mode_for_flag(mode,CALC_STAT_GLOBAL)){
        MPI_Allreduce(&result_local_d,   &result_d,   1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(&n_data_used_local,&n_data_used,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      }
#endif
      break;
    case ADaM_FLOAT:
      if(mask_local!=NULL){
        for(i_data=0,flag_init=TRUE;flag_init && i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_d=((float *)data_local)[i_data];
            flag_init     =FALSE;
          }
        }
        for(n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_d=MAX(result_local_d,((float *)data_local)[i_data]);
            n_data_used_local++;
          }
        }
      }
      else{
        result_local_d=((float *)data_local)[0];
        for(i_data=1;i_data<n_data_local;i_data++)
          result_local_d=MAX(result_local_d,((float *)data_local)[i_data]);
        n_data_used_local=n_data_local;
      }
      result_d   =result_local_d;
      n_data_used=n_data_used_local;
#ifdef USE_MPI
      // STILL NEED TO CHECK HERE THAT ALL PROCS HAVE N_DATA_LOCAL>0 ... INCORRECT RESULTS OTHERWISE!
      if(check_mode_for_flag(mode,CALC_STAT_GLOBAL)){
        MPI_Allreduce(&result_local_d,   &result_d,   1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(&n_data_used_local,&n_data_used,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      }
#endif
      break;
    case ADaM_INT:
      if(mask_local!=NULL){
        for(i_data=0,flag_init=TRUE;flag_init && i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_t=((int *)data_local)[i_data];
            flag_init     =FALSE;
          }
        }
        for(n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_t=MAX(result_local_t,((int *)data_local)[i_data]);
            n_data_used_local++;
          }
        }
      }
      else{
        result_local_t=((int *)data_local)[0];
        for(i_data=1;i_data<n_data_local;i_data++)
          result_local_t=MAX(result_local_t,((int *)data_local)[i_data]);
        n_data_used_local=n_data_local;
      }
      result_t   =result_local_t;
      n_data_used=n_data_used_local;
#ifdef USE_MPI
      // STILL NEED TO CHECK HERE THAT ALL PROCS HAVE N_DATA_LOCAL>0 ... INCORRECT RESULTS OTHERWISE!
      if(check_mode_for_flag(mode,CALC_STAT_GLOBAL)){
        MPI_Allreduce(&result_local_t,   &result_t,   1,MPI_SIZE_T,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(&n_data_used_local,&n_data_used,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      }
#endif
      result_d=(double)result_t;
      break;
    case ADaM_SIZE_T:
      if(mask_local!=NULL){
        for(i_data=0,flag_init=TRUE;flag_init && i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_t=((size_t *)data_local)[i_data];
            flag_init     =FALSE;
          }
        }
        for(n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data]){
            result_local_t=MAX(result_local_t,((size_t *)data_local)[i_data]);
            n_data_used_local++;
          }
        }
      }
      else{
        result_local_t=((size_t *)data_local)[0];
        for(i_data=1;i_data<n_data_local;i_data++)
          result_local_t=MAX(result_local_t,((size_t *)data_local)[i_data]);
        n_data_used_local=n_data_local;
      }
      result_t   =result_local_t;
      n_data_used=n_data_used_local;
#ifdef USE_MPI
      // STILL NEED TO CHECK HERE THAT ALL PROCS HAVE N_DATA_LOCAL>0 ... INCORRECT RESULTS OTHERWISE!
      if(check_mode_for_flag(mode,CALC_STAT_GLOBAL)){
        MPI_Allreduce(&result_local_t,   &result_t,   1,MPI_SIZE_T,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(&n_data_used_local,&n_data_used,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      }
#endif
      result_d=(double)result_t;
      break;      
    default:
      SID_trap_error("Unknown variable type {%d} in calc_stat\n",ERROR_LOGIC,type);
      break;
    }
    break;
    
    // ----------------------------------------------
    // Compute the **STANDARD DEVIATION** of the data
    // ----------------------------------------------
  case CALC_STAT_STDDEV:
    if(check_mode_for_flag(mode,CALC_STAT_GLOBAL))
      n_data_used=calc_stat(data_local,
                            mask_local,
                            n_data_local,
                            type,
                            CALC_STAT_MEAN|CALC_STAT_GLOBAL|CALC_STAT_RETURN_DOUBLE,
                            &mean);
    else
      n_data_used=calc_stat(data_local,
                            mask_local,
                            n_data_local,
                            type,
                            CALC_STAT_MEAN|CALC_STAT_RETURN_DOUBLE,
                            &mean);
    switch(type){
      result_local_d=0.;
    case ADaM_DOUBLE:
      if(mask_local!=NULL){
        for(i_data=0,n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data])
            result_local_d+=pow((double)(((double *)data_local)[i_data])-mean,2.);
        }
      }
      else{
        for(i_data=0,n_data_used_local=0;i_data<n_data_local;i_data++)
          result_local_d+=pow((double)(((double *)data_local)[i_data])-mean,2.);
      }
      break;
    case ADaM_FLOAT:
      if(mask_local!=NULL){
        for(i_data=0,n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data])
            result_local_d+=pow((double)(((float *)data_local)[i_data])-mean,2.);
        }
      }
      else{
        for(i_data=0,n_data_used_local=0;i_data<n_data_local;i_data++)
          result_local_d+=pow((double)(((float *)data_local)[i_data])-mean,2.);
      }
      break;
    case ADaM_INT:
      if(mask_local!=NULL){
        for(i_data=0,n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data])
            result_local_d+=pow((double)(((int *)data_local)[i_data])-mean,2.);
        }
      }
      else{
        for(i_data=0,n_data_used_local=0;i_data<n_data_local;i_data++)
          result_local_d+=pow((double)(((int *)data_local)[i_data])-mean,2.);
      }
      break;
    case ADaM_SIZE_T:
      if(mask_local!=NULL){
        for(i_data=0,n_data_used_local=0;i_data<n_data_local;i_data++){
          if(mask_local[i_data])
            result_local_d+=pow((double)(((size_t *)data_local)[i_data])-mean,2.);
        }
      }
      else{
        for(i_data=0,n_data_used_local=0;i_data<n_data_local;i_data++)
          result_local_d+=pow((double)(((size_t *)data_local)[i_data])-mean,2.);
      }
      break;      
    default:
      SID_trap_error("Unknown variable type {%d} in calc_stat\n",ERROR_LOGIC,type);
      break;
    }
    result_d=result_local_d;
#ifdef USE_MPI
    if(check_mode_for_flag(mode,CALC_STAT_GLOBAL))
      MPI_Allreduce(&result_local_d,&result_d,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    result_d=sqrt(result_d/(double)n_data_used);
    result_t=(size_t)result_d;
    break;
    
    // ----------------------------------
    // Compute the **MEDIAN** of the data
    // ----------------------------------
  case CALC_STAT_MEDIAN:
    if(!check_mode_for_flag(mode,CALC_STAT_GLOBAL) || SID.n_proc==1){
      flag_local=SORT_LOCAL;
      sort(data_local,
           n_data_local,
           &sort_index_local,
           type,
           flag_local,
           SORT_COMPUTE_INDEX,
           FALSE);

      i_mid=n_data_local/2;
      if(n_data%2){
        switch(type){
        case ADaM_DOUBLE:
          result_d=(double)((double *)data_local)[sort_index_local[i_mid]];
          break;
        case ADaM_FLOAT:
          result_d=(double)((float *)data_local)[sort_index_local[i_mid]];
          break;
        case ADaM_INT:
          result_t=(size_t)((int *)data_local)[sort_index_local[i_mid]];
          break;
        case ADaM_SIZE_T:
          result_t=(size_t)((int *)data_local)[sort_index_local[i_mid]];
          break;
        default:
          SID_trap_error("type {%d} not supported in calc_stat.\n",ERROR_LOGIC,type);
          break;
        }
      }
      else{
        switch(type){
        case ADaM_DOUBLE:
          result_d=(double)(ONE_HALF*(((double *)data_local)[sort_index_local[i_mid-1]]+
                                      ((double *)data_local)[sort_index_local[i_mid]]));
          break;
        case ADaM_FLOAT:
          result_d=(double)(ONE_HALF*(((float *)data_local)[sort_index_local[i_mid-1]]+
                                      ((float *)data_local)[sort_index_local[i_mid]]));
          break;
        case ADaM_INT:
          result_t=(size_t)(ONE_HALF*(((int *)data_local)[sort_index_local[i_mid-1]]+
                                      ((int *)data_local)[sort_index_local[i_mid]]));
          break;
        case ADaM_SIZE_T:
          result_t=(size_t)(ONE_HALF*(((size_t *)data_local)[sort_index_local[i_mid-1]]+
                                      ((size_t *)data_local)[sort_index_local[i_mid]]));
          break;
        default:
          SID_trap_error("type {%d} not supported in calc_median.\n",ERROR_LOGIC,type);
          break;
        }
      }
      free(sort_index_local);
    }
#ifdef USE_MPI
    else{
      SID_trap_error("parallel global medians have not been debugged yet!",ERROR_LOGIC);
      MPI_Allreduce(&n_data_local,&n_data,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
      if(n_data>0){
        if(n_data>1){
          sort(data_local,
               n_data_local,
               &sort_index_local,
               type,
               FALSE,
               TRUE,
               FALSE);
          index_1=n_data/2-1;
          index_2=n_data/2;
          sort(sort_index_local,
               n_data_local,
               &idx,
               ADaM_SIZE_T,
               TRUE,
               TRUE,
               FALSE);
          i_lo    =0;
          i_hi    =n_data_local-1;
          index_lo=sort_index_local[idx[i_lo]];
          index_hi=sort_index_local[idx[i_hi]];
          do{     
            i_mid    =(i_lo+i_hi)/2;
            index_mid=sort_index_local[idx[i_mid]];
            if(index_1<index_mid)
              i_hi=i_mid;
            else if(index_1>index_mid)
              i_lo=i_mid;
            else{ 
              i_lo=i_mid;
              i_hi=i_mid;
            }     
          } while(i_hi-i_lo>2);
          if(index_1==sort_index_local[i_lo]){
            i_1   =i_lo;
            flag_1=TRUE;
          }
          else
            flag_1=FALSE;
          if(index_2==sort_index_local[i_lo]){
            i_2   =i_lo;
            flag_2=TRUE;
          }
          else if(index_2==sort_index_local[i_hi]){
            i_2   =i_hi;
            flag_2=TRUE;
          }
          else
            flag_2=FALSE;
          if(!(n_data%2)){
            if(flag_1){
              rank_1_local=SID.My_rank;
              switch(type){
              case ADaM_DOUBLE:
                value_1_d=(double)(((double *)data_local)[i_1]);
                break;
              case ADaM_FLOAT:
                value_1_d=(double)(((float  *)data_local)[i_1]);
                break;
              case ADaM_INT:
                value_1_t=(size_t)(((int    *)data_local)[i_1]);
                break;
              case ADaM_SIZE_T:
                value_1_t=(size_t)(((size_t *)data_local)[i_1]);
                break;
              default:
                SID_trap_error("type {%d} not supported by calc_median_global.\n",ERROR_LOGIC,type);
                break;
              }
            }
            else
              rank_1_local=-1;
            MPI_Allreduce(&rank_1_local,&rank_1,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
            MPI_Bcast(&value_1_t,1,MPI_SIZE_T,rank_1,MPI_COMM_WORLD);
            MPI_Bcast(&value_1_d,1,MPI_DOUBLE,rank_1,MPI_COMM_WORLD);
          }
          if(flag_2){
            rank_2_local=SID.My_rank;
            switch(type){
            case ADaM_DOUBLE:
              value_2_d=(double)(((double *)data_local)[i_2]);
              break;
            case ADaM_FLOAT:
              value_2_d=(double)(((float  *)data_local)[i_2]);
              break;
            case ADaM_INT:
              value_2_t=(size_t)(((int    *)data_local)[i_2]);
              break;
            case ADaM_SIZE_T:
              value_2_t=(size_t)(((size_t *)data_local)[i_2]);
              break;
            default:
              SID_trap_error("type {%d} not supported by calc_median_global.\n",ERROR_LOGIC,type);
              break;
            }
          }
          else
            rank_2_local=-1;
          MPI_Allreduce(&rank_2_local,&rank_2,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
          MPI_Bcast(&value_2_t,1,MPI_SIZE_T,rank_2,MPI_COMM_WORLD);          
          MPI_Bcast(&value_2_d,1,MPI_DOUBLE,rank_2,MPI_COMM_WORLD);          
          free(idx);
          free(sort_index_local);

          switch(type){
          case ADaM_DOUBLE:
          case ADaM_FLOAT:
            if(n_data%2)
              result_d=value_2_d;
            else
              result_d=ONE_HALF*(value_1_d+value_2_d);
            break;
          case ADaM_INT:
          case ADaM_SIZE_T:
            if(n_data%2)
              result_t=value_2_t;
            else
              result_t=ONE_HALF*(value_1_t+value_2_t);
            break;
          default:
            SID_trap_error("Unknown variable type {%d} in calc_stat\n",ERROR_LOGIC,type);
            break;
          }
        }
        else
          switch(type){
          case ADaM_DOUBLE:
            result_d=(double)(((double *)data_local)[0]);
            break;
          case ADaM_FLOAT:
            result_d=(double)(((float  *)data_local)[0]);
            break;
          case ADaM_INT:
            result_t=(size_t)(((int    *)data_local)[0]);
            break;
          case ADaM_SIZE_T:
            result_t=(size_t)(((size_t *)data_local)[0]);
            break;
          default:
            SID_trap_error("Unknown variable type {%d} in calc_stat\n",ERROR_LOGIC,type);
            break;
          }
      }
      else{
        result_t=0;
        result_d=0.;
      }
    }
#endif
    break;
  default:
    SID_trap_error("Operation not specified correctly in calc_stat. (mode=%lld)",ERROR_LOGIC,mode);
    break;
  }
  
  // At last, set the return values
  if(check_mode_for_flag(mode,CALC_STAT_RETURN_DOUBLE))
    ((double *)return_value)[0]=(double)result_d;
  else{
    switch(type){
    case ADaM_DOUBLE:
      ((double *)return_value)[0]=(double)result_d;
      break;
    case ADaM_FLOAT:
      ((float *)return_value)[0] =(float)result_d;
      break;
    case ADaM_INT:
      ((int *)return_value)[0]   =(int)result_t;
      break;
    case ADaM_SIZE_T:
      ((size_t *)return_value)[0]=(size_t)result_t;
      break;
    default:
      SID_trap_error("Unknown variable type {%d} in calc_stat\n",ERROR_LOGIC,type);
      break;
    }
  }
  return(n_data_used);
}
