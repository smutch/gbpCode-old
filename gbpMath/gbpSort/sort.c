#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpDomain.h>
#include <gbpSort.h>

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
  size_t  nval_max;
  size_t  nval_tmp;
  size_t  nval_all;
  size_t *index_tmp;
  size_t *increment;
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
    // Rank local array in ascending order; must
    //   use merge_sort (despite extra RAM use) 
    //   because it is "stable"
    SID_log("Sorting local items...",SID_LOG_OPEN|SID_LOG_TIMER);
    merge_sort(sval,
	       nval,
	       index,
	       data_type,
	       SORT_COMPUTE_RANK,
	       SORT_COMPUTE_NOT_INPLACE);
    SID_log("Done.",SID_LOG_CLOSE);

    // If this is being run in parallel, consider
    //   the arrays on the other ranks as well
    //   if flag_local is set to SORT_GLOBAL
    #if USE_MPI
    if(flag_local==SORT_GLOBAL && SID.n_proc>1){
      SID_log("Sorting distributed items...",SID_LOG_OPEN|SID_LOG_TIMER);
      // ... get the largest number of items on any rank ...
      SID_Allreduce(&nval,&nval_max,1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
      // ... get the size of our data-type ...
      SID_Type_size(data_type,&data_type_size_i);
      data_type_size=(size_t)data_type_size_i;
      // ... create two temporary directories that are needed for the sort ...
      sval_tmp =(void   *)SID_malloc((size_t)data_type_size*(nval+nval_max));
      increment=(size_t *)SID_calloc(sizeof(size_t)*nval);

      // ... perform exchanges and sorts ...
      for(i_rank=1;i_rank<SID.n_proc;i_rank++){
        SID_log("Rank exchange %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_rank,SID.n_proc-1);

	// Ranks < My_rank are placed first in the new data list
	//  so that duplicate entries are placed in a unique order
        set_exchange_ring_ranks(&rank_send_to,&rank_receive_from,i_rank);
        if(rank_receive_from<SID.My_rank){
           exchange_ring_buffer(sval,
                                data_type_size,
                                nval,
                                sval_tmp,
                                &nval_tmp,
                                i_rank);
           memcpy(&(((char *)sval_tmp)[nval_tmp*data_type_size]),sval,nval*data_type_size);
        }
        else{
           exchange_ring_buffer(sval,
                                data_type_size,
                                nval,
                                &(((char *)sval_tmp)[nval*data_type_size]),
                                &nval_tmp,
                                i_rank);
           memcpy(sval_tmp,sval,nval*data_type_size);
        }
        nval_all=nval+nval_tmp;

        // A STABLE sort algorythm MUST be used here!
	merge_sort(sval_tmp,
		   nval_all,
		   &index_tmp,
		   data_type,
		   SORT_COMPUTE_RANK,
		   SORT_COMPUTE_INPLACE);

	// Figure-out where the local values lie in the new sort
	//   and increment rank for each if any values from
	//   receive_from_rank's array are less than it
	if(rank_receive_from<SID.My_rank){
	  for(i=0,j=0;i<nval;i++)
	    increment[i]+=index_tmp[i+nval_tmp]-(*index)[i];
	}
	else{
	  for(i=0,j=0;i<nval;i++)
	    increment[i]+=index_tmp[i]-(*index)[i];
	}

        // Free sort ranks
	SID_free(SID_FARG index_tmp);
        SID_log("Done.",SID_LOG_CLOSE);
      } // loop over ranks

      for(i=0;i<nval;i++)
	(*index)[i]+=increment[i];
    
      SID_free(SID_FARG increment);
      SID_free(SID_FARG sval_tmp);
      SID_log("Done.",SID_LOG_CLOSE);
    } // if rank exchange is needed

    if(flag_local        ==SORT_GLOBAL && 
       flag_compute_index!=SORT_COMPUTE_RANK){
      size_t *sort_ranks;
      SID_log("Generating sort indices from sort ranks...",SID_LOG_OPEN|SID_LOG_TIMER);
      // ... get the largest number of items on any rank ...
      SID_Allreduce(&nval,&nval_max,1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
      // ... get the size of our data-type ...
      SID_Type_size(data_type,&data_type_size_i);
      data_type_size=(size_t)data_type_size_i;
      // ... determine what the first element of the sort indices
      //     is for each rank ..
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
      (*index)  =(size_t *)SID_malloc(sizeof(size_t)*nval);
      rank_rank =(size_t *)SID_malloc(sizeof(size_t)*nval_max);
      // ... perform exchanges and set indices ...
      for(i_rank=0,offset=0;i_rank<SID.n_proc;i_rank++){
        nval_tmp=nval;
        SID_Bcast(&nval_tmp,sizeof(size_t),i_rank,SID.COMM_WORLD);
        if(i_rank==SID.My_rank)
           memcpy(rank_rank,sort_ranks,nval*sizeof(size_t));
        SID_Bcast(rank_rank,nval_tmp*sizeof(size_t),i_rank,SID.COMM_WORLD);
        // ... scan the sort ranks we've just received and see if any of them
        //     are supposed to be pointed to by the indices stored locally ...
        for(i_val=0;i_val<nval_tmp;i_val++){
          rank_rank[i_val]-=first_index;
          if(rank_rank[i_val]>=0 && rank_rank[i_val]<nval)
            (*index)[rank_rank[i_val]]=offset+i_val;
        }
        offset+=nval_tmp;
      }

      SID_free(SID_FARG sort_ranks);
      SID_free(SID_FARG rank_rank);
      SID_log("Done.",SID_LOG_CLOSE);
    }
    #endif
  }
  SID_log("Done.",SID_LOG_CLOSE);
}
