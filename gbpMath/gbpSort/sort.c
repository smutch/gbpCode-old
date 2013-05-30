#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpSort.h>
#include <gbpDomain.h>

void sort(void          *sval,
	  size_t         nval,
	  size_t       **index,
	  SID_Datatype   data_type,
	  int            flag_local,
	  int            flag_compute_index,
	  int            flag_in_place){
  int     i_rank;
  int     j_rank;
  int     rank_send_to;
  int     rank_receive_from;
  size_t  i,ir,j,l,k;
  size_t  nval_tmp;
  size_t  nval_all;
  size_t *index_tmp;
  size_t *rank;
  int     data_type_size_i;
  size_t  data_type_size;
  void   *sval_tmp;

  //SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  SID_log("Performing sort...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Process passed arguments:
  //   ... check for nonsensical flag combinations
  if(flag_local   ==SORT_GLOBAL && 
     flag_in_place!=SORT_COMPUTE_NOT_INPLACE)
    SID_trap_error("Global in-place sorts are not possible.",ERROR_LOGIC);
  //   ... use a heap_sort for most situations (it uses less RAM!)
  else if(flag_local==SORT_LOCAL || SID.n_proc==1){
/*
    heap_sort(sval,
	      nval,
	      index,
	      data_type,
	      flag_compute_index,
	      flag_in_place);
*/
    merge_sort(sval,
              nval,
              index,
              data_type,
              flag_compute_index,
              flag_in_place);

  }
  //   ... perform a global search returning ranks
  else{
    #if USE_MPI
    // Rank local array in ascending order; must
    //   use merge_sort (despite extra RAM use) 
    //   because it is "stable"
    SID_log("Sorting local items...",SID_LOG_OPEN|SID_LOG_TIMER);
    merge_sort(sval,
	       nval,
	       index,
	       data_type,
	       SORT_COMPUTE_INDEX,
	       SORT_COMPUTE_NOT_INPLACE);

    // If this is being run in parallel, consider
    //   the arrays on the other ranks as well
    //   if flag_local is set to SORT_GLOBAL

    // ... create a local array to hold results for the local rank ...
    size_t *rank=NULL;

    // ... initialize it with the local rank values ...
    merge_sort(sval,
	       nval,
	       &rank,
	       data_type,
	       SORT_COMPUTE_RANK,
	       SORT_COMPUTE_NOT_INPLACE);
    SID_log("Done.",SID_LOG_CLOSE);

    // ... get the largest number of items on any rank ...
    SID_Allreduce(&nval,&nval_tmp,1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);

    // ... get the size of our data-type ...
    SID_Type_size(data_type,&data_type_size_i);
    data_type_size=(size_t)data_type_size_i;

    // Process the items on other ranks
    if(flag_local==SORT_GLOBAL && SID.n_proc>1){
      SID_log("Sorting distributed items...",SID_LOG_OPEN|SID_LOG_TIMER);

      // ... create two temporary arrays for exchanges ...
      sval_tmp =(void   *)SID_malloc((size_t)data_type_size*nval_tmp);
      index_tmp=(size_t *)SID_malloc(sizeof(size_t)*nval_tmp);

      // ... perform exchanges and sorts ...
      for(i_rank=1;i_rank<SID.n_proc;i_rank++){
        SID_log("Rank exchange %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_rank,SID.n_proc-1);

	// Perform a ring exchange of all data and their local sort indices
        set_exchange_ring_ranks(&rank_send_to,&rank_receive_from,i_rank);
        exchange_ring_buffer(sval,
                             data_type_size,
                             nval,
                             sval_tmp,
                             &nval_tmp,
                             i_rank);
        exchange_ring_buffer((*index),
                             sizeof(size_t),
                             nval,
                             index_tmp,
                             &nval_tmp,
                             i_rank);

        // Increment the local sort rank by dispplacements computed from other ranks
        if(data_type==SID_DOUBLE){
           double *sval_dtype    =(double *)sval;
           double *sval_tmp_dtype=(double *)sval_tmp;
	   if(rank_receive_from<SID.My_rank){
	     for(i=0,j=0;i<nval;i++){
               if(j<nval_tmp){
                  // Equal values from lower ranks come before (hence '>=')
                  while(sval_dtype[(*index)[i]]>=sval_tmp_dtype[index_tmp[j]]){
                    j++;
                    if(j==nval_tmp)
                      break;
                  }
               }
               rank[(*index)[i]]+=j;
             }
	   }
	   else{
	     for(i=0,j=0;i<nval;i++){
               if(j<nval_tmp){
                  // Equal values from larger ranks come after (hence '>')
                  while(sval_dtype[(*index)[i]]>sval_tmp_dtype[index_tmp[j]]){
                    j++;
                    if(j==nval_tmp)
                      break;
                  }
               }
               rank[(*index)[i]]+=j;
             }
	   }
        }
        else if(data_type==SID_FLOAT){
           float *sval_dtype    =(float *)sval;
           float *sval_tmp_dtype=(float *)sval_tmp;
	   if(rank_receive_from<SID.My_rank){
	     for(i=0,j=0;i<nval;i++){
               if(j<nval_tmp){
                  // Equal values from lower ranks come before (hence '>=')
                  while(sval_dtype[(*index)[i]]>=sval_tmp_dtype[index_tmp[j]]){
                    j++;
                    if(j==nval_tmp)
                      break;
                  }
               }
               rank[(*index)[i]]+=j;
             }
	   }
	   else{
	     for(i=0,j=0;i<nval;i++){
               if(j<nval_tmp){
                  // Equal values from larger ranks come after (hence '>')
                  while(sval_dtype[(*index)[i]]>sval_tmp_dtype[index_tmp[j]]){
                    j++;
                    if(j==nval_tmp)
                      break;
                  }
               }
               rank[(*index)[i]]+=j;
             }
	   }
        }
        else if(data_type==SID_INT){
           int *sval_dtype    =(int *)sval;
           int *sval_tmp_dtype=(int *)sval_tmp;
	   if(rank_receive_from<SID.My_rank){
	     for(i=0,j=0;i<nval;i++){
               if(j<nval_tmp){
                  // Equal values from lower ranks come before (hence '>=')
                  while(sval_dtype[(*index)[i]]>=sval_tmp_dtype[index_tmp[j]]){
                    j++;
                    if(j==nval_tmp)
                      break;
                  }
               }
               rank[(*index)[i]]+=j;
             }
	   }
	   else{
	     for(i=0,j=0;i<nval;i++){
               if(j<nval_tmp){
                  // Equal values from larger ranks come after (hence '>')
                  while(sval_dtype[(*index)[i]]>sval_tmp_dtype[index_tmp[j]]){
                    j++;
                    if(j==nval_tmp)
                      break;
                  }
               }
               rank[(*index)[i]]+=j;
             }
	   }
        }
        else if(data_type==SID_UNSIGNED){
           unsigned int *sval_dtype    =(unsigned int *)sval;
           unsigned int *sval_tmp_dtype=(unsigned int *)sval_tmp;
	   if(rank_receive_from<SID.My_rank){
	     for(i=0,j=0;i<nval;i++){
               if(j<nval_tmp){
                  // Equal values from lower ranks come before (hence '>=')
                  while(sval_dtype[(*index)[i]]>=sval_tmp_dtype[index_tmp[j]]){
                    j++;
                    if(j==nval_tmp)
                      break;
                  }
               }
               rank[(*index)[i]]+=j;
             }
	   }
	   else{
	     for(i=0,j=0;i<nval;i++){
               if(j<nval_tmp){
                  // Equal values from larger ranks come after (hence '>')
                  while(sval_dtype[(*index)[i]]>sval_tmp_dtype[index_tmp[j]]){
                    j++;
                    if(j==nval_tmp)
                      break;
                  }
               }
               rank[(*index)[i]]+=j;
             }
	   }
        }
        else if(data_type==SID_SIZE_T){
           size_t *sval_dtype    =(size_t *)sval;
           size_t *sval_tmp_dtype=(size_t *)sval_tmp;
	   if(rank_receive_from<SID.My_rank){
	     for(i=0,j=0;i<nval;i++){
               if(j<nval_tmp){
                  // Equal values from lower ranks come before (hence '>=')
                  while(sval_dtype[(*index)[i]]>=sval_tmp_dtype[index_tmp[j]]){
                    j++;
                    if(j==nval_tmp)
                      break;
                  }
               }
               rank[(*index)[i]]+=j;
             }
	   }
	   else{
	     for(i=0,j=0;i<nval;i++){
               if(j<nval_tmp){
                  // Equal values from larger ranks come after (hence '>')
                  while(sval_dtype[(*index)[i]]>sval_tmp_dtype[index_tmp[j]]){
                    j++;
                    if(j==nval_tmp)
                      break;
                  }
               }
               rank[(*index)[i]]+=j;
             }
	   }
        }
        else if(data_type==SID_LONG_LONG){
           long long *sval_dtype    =(long long *)sval;
           long long *sval_tmp_dtype=(long long *)sval_tmp;
	   if(rank_receive_from<SID.My_rank){
	     for(i=0,j=0;i<nval;i++){
               if(j<nval_tmp){
                  // Equal values from lower ranks come before (hence '>=')
                  while(sval_dtype[(*index)[i]]>=sval_tmp_dtype[index_tmp[j]]){
                    j++;
                    if(j==nval_tmp)
                      break;
                  }
               }
               rank[(*index)[i]]+=j;
             }
	   }
	   else{
	     for(i=0,j=0;i<nval;i++){
               if(j<nval_tmp){
                  // Equal values from larger ranks come after (hence '>')
                  while(sval_dtype[(*index)[i]]>sval_tmp_dtype[index_tmp[j]]){
                    j++;
                    if(j==nval_tmp)
                      break;
                  }
               }
               rank[(*index)[i]]+=j;
             }
	   }
        }
        else
          SID_trap_error("Unknown variable type in sort().",ERROR_LOGIC);

        SID_log("Done.",SID_LOG_CLOSE);
      } // loop over ranks

      SID_free(SID_FARG sval_tmp);
      SID_free(SID_FARG index_tmp);
      SID_log("Done.",SID_LOG_CLOSE);
    } // if rank exchange is needed
     
    // We want global ranks to be stored in 'index' at this point
    SID_free(SID_FARG (*index));
    (*index)=rank;

    // If we want a global index, then we need to sort the global ranks stored in 'index'
    if(flag_local        ==SORT_GLOBAL && 
       flag_compute_index!=SORT_COMPUTE_RANK){
      size_t *sort_ranks;
      SID_log("Generating sort indices from sort ranks...",SID_LOG_OPEN|SID_LOG_TIMER);
      // ... determine what the first array index is for each rank ..
      size_t first_index;
      for(i_rank=0,first_index=0;i_rank<SID.n_proc;i_rank++){
        nval_tmp=nval;
        SID_Bcast(&nval_tmp,sizeof(size_t),i_rank,SID.COMM_WORLD);
        if(i_rank<SID.My_rank)
          first_index+=nval_tmp;
      }
      // ... move the sort ranks (generated above) to a temp array and allocate
      //     a new array for the sort indices ...
      size_t *rank_rank;
      size_t  offset;
      size_t  i_val;
      sort_ranks=(*index);
      (*index)  =(size_t *)SID_calloc(sizeof(size_t)*nval);
      rank_rank =(size_t *)SID_malloc(sizeof(size_t)*nval_tmp);
      // ... loop over ranks, performing exchanges and setting indices ...
      for(i_rank=0,offset=0;i_rank<SID.n_proc;i_rank++){
        nval_tmp=nval;
        SID_Bcast(&nval_tmp,sizeof(size_t),i_rank,SID.COMM_WORLD);
        if(i_rank==SID.My_rank)
           memcpy(rank_rank,sort_ranks,nval*sizeof(size_t));
        SID_Bcast(rank_rank,nval_tmp*sizeof(size_t),i_rank,SID.COMM_WORLD);
        // ... scan the sort ranks we've just received and see if any of them
        //     are supposed to be pointed to by the indices stored locally ...
        for(i_val=0;i_val<nval_tmp;i_val++){
          if(rank_rank[i_val]>=first_index && rank_rank[i_val]<(nval+first_index))
            (*index)[rank_rank[i_val]-first_index]=offset+i_val;
        }
        offset+=nval_tmp;
      }
      SID_free(SID_FARG sort_ranks);
      SID_free(SID_FARG rank_rank);
      SID_log("Done.",SID_LOG_CLOSE);
    }
    #else
      SID_trap_error("Undefined behavior in sort().",ERROR_LOGIC);
    #endif
  }
  SID_log("Done.",SID_LOG_CLOSE);
}

