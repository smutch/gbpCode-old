#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpSort.h>

void heap_sort(void    *data_in,
               size_t   n_data,
               size_t **index,
               SID_Datatype data_type,
               int      flag_compute_index,
               int      flag_in_place){
  size_t  i,ir,j,l;
  void   *rra;
  size_t  rrb;
  void   *data;
  size_t  data_size;
  int     flag_indices_used;
  size_t *idx_sort;

  // Process passed arguments:
  //   Determine data element size
  if(data_type==SID_INT)
    data_size=sizeof(int);
  else if(data_type==SID_SIZE_T)
    data_size=sizeof(size_t);
  else if(data_type==SID_FLOAT)
    data_size=sizeof(float);
  else if(data_type==SID_DOUBLE)
    data_size=sizeof(double);
  else
    SID_trap_error("Unsupported data type {%d}",ERROR_LOGIC,data_type);

  // Allocate scratch space and sort indices (if needed)
  rra=SID_malloc(data_size);
  if(flag_in_place==SORT_COMPUTE_INPLACE)
    data=data_in;
  else if(flag_in_place==SORT_COMPUTE_NOT_INPLACE){
    data=SID_malloc(n_data*data_size);
    memcpy(data,data_in,n_data*data_size);
  }
  else{
    fprintf(stderr,"ERROR (heap_sort): flag_in_place {%d} must be SORT_COMPUTE_INPLACE||SORT_COMPUTE_NOT_INPLACE\n",
	    flag_in_place);
  }
  flag_indices_used=FALSE;
  if(flag_compute_index==SORT_COMPUTE_INDEX ||
     flag_compute_index==SORT_COMPUTE_RANK){
    flag_indices_used=TRUE;
    idx_sort =(size_t *)SID_malloc(sizeof(size_t)*n_data);
    for(i=0;i<n_data;i++)
      idx_sort[i]=i;
  }
  else if(flag_compute_index!=SORT_INPLACE_ONLY){
    fprintf(stderr,"ERROR (heap_sort): flag_compute_index {%d} must be SORT_COMPUTE_INDEX||SORT_COMPUTE_RANK||SORT_INPLACE_ONLY\n",
	    flag_compute_index);
  }

  // Check for nonsensical flag combinations
  if(flag_in_place==SORT_COMPUTE_NOT_INPLACE && flag_compute_index==SORT_INPLACE_ONLY){
    fprintf(stderr,"ERROR (heap_sort): flag combinations will generate no results.\n");
  }

  // Sort in ascending order -- modified version of NUMREC code hpsort.c
  l=(n_data >> 1)+1;
  ir=n_data;
  for (;;) {
    if (l > 1) {
      l--;
      if(data_type==SID_INT)
	((int *)rra)[0]=((int *)data)[l-1];
      else if(data_type==SID_SIZE_T)
	((size_t *)rra)[0]=((size_t *)data)[l-1];
      else if(data_type==SID_FLOAT)
	((float *)rra)[0]=((float *)data)[l-1];
      else if(data_type==SID_DOUBLE)
	((double *)rra)[0]=((double *)data)[l-1];
      if(flag_indices_used)
	rrb=idx_sort[l-1];
    } 
    else {
      if(data_type==SID_INT){
	((int *)rra)[0]    =((int *)data)[ir-1];
	((int *)data)[ir-1]=((int *)data)[0];
      }
      else if(data_type==SID_SIZE_T){
	((size_t *)rra)[0]    =((size_t *)data)[ir-1];
	((size_t *)data)[ir-1]=((size_t *)data)[0];
      }
      else if(data_type==SID_FLOAT){
	((float *)rra)[0]    =((float *)data)[ir-1];
	((float *)data)[ir-1]=((float *)data)[0];
      }
      else if(data_type==SID_DOUBLE){
	((double *)rra)[0]    =((double *)data)[ir-1];
	((double *)data)[ir-1]=((double *)data)[0];
      }
      if(flag_indices_used){
	rrb=idx_sort[ir-1];
	idx_sort[ir-1] =idx_sort[0];
      }
      if (--ir == 1) {
	if(data_type==SID_INT)
	  ((int *)data)[0]=((int *)rra)[0];
	else if(data_type==SID_SIZE_T)
	  ((size_t *)data)[0]=((size_t *)rra)[0];
	else if(data_type==SID_FLOAT)
	  ((float *)data)[0]=((float *)rra)[0];
	else if(data_type==SID_DOUBLE)
	  ((double *)data)[0]=((double *)rra)[0];
	if(flag_indices_used){
	  idx_sort[ir-1] =rrb;
	}
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
	if(data_type==SID_INT){
	  if (j < ir && ((int *)data)[j-1] < ((int *)data)[j]) j++;
	  if (((int *)rra)[0] < ((int *)data)[j-1]) {
	    ((int *)data)[i-1]=((int *)data)[j-1];
	    if(flag_indices_used)
	      idx_sort[i-1] =idx_sort[j-1];
	    i=j;
	    j <<= 1;
	  } 
	  else 
	    j=ir+1;
	}
	else if(data_type==SID_SIZE_T){
	  if (j < ir && ((size_t *)data)[j-1] < ((size_t *)data)[j]) j++;
	  if (((size_t *)rra)[0] < ((size_t *)data)[j-1]) {
	    ((size_t *)data)[i-1]=((size_t *)data)[j-1];
	    if(flag_indices_used)
	      idx_sort[i-1] =idx_sort[j-1];
	    i=j;
	    j <<= 1;
	  } 
	  else 
	    j=ir+1;
	}
	else if(data_type==SID_FLOAT){
	  if (j < ir && ((float *)data)[j-1] < ((float *)data)[j]) j++;
	  if (((float *)rra)[0] < ((float *)data)[j-1]) {
	    ((float *)data)[i-1]=((float *)data)[j-1];
	    if(flag_indices_used)
	      idx_sort[i-1] =idx_sort[j-1];
	    i=j;
	    j <<= 1;
	  } 
	  else 
	    j=ir+1;
	}
	else if(data_type==SID_DOUBLE){
	  if (j < ir && ((double *)data)[j-1] < ((double *)data)[j]) j++;
	  if (((double *)rra)[0] < ((double *)data)[j-1]) {
	    ((double *)data)[i-1]=((double *)data)[j-1];
	    if(flag_indices_used)
	      idx_sort[i-1] =idx_sort[j-1];
	    i=j;
	    j <<= 1;
	  } 
	  else 
	    j=ir+1;
	}
    }
    if(data_type==SID_INT){
      ((int *)data)[i-1]=((int *)rra)[0];
      if(flag_indices_used)
	idx_sort[i-1]=rrb;
    }
    else if(data_type==SID_SIZE_T){
      ((size_t *)data)[i-1]=((size_t *)rra)[0];
      if(flag_indices_used)
	idx_sort[i-1]=rrb;
    }
    else if(data_type==SID_FLOAT){
      ((float *)data)[i-1]=((float *)rra)[0];
      if(flag_indices_used)
	idx_sort[i-1]=rrb;
    }
    else if(data_type==SID_DOUBLE){
      ((double *)data)[i-1]=((double *)rra)[0];
      if(flag_indices_used)
	idx_sort[i-1]=rrb;
    }
  }
  if(!flag_in_place)
    SID_free((void **)&data);

  // Convert sort indices to ranks (if needed)
  if(flag_compute_index==SORT_COMPUTE_RANK){
    heap_sort(idx_sort,
              n_data,
              index,
	      SID_SIZE_T,
              SORT_COMPUTE_INDEX,
              SORT_COMPUTE_INPLACE);
    SID_free((void **)&idx_sort);
  }
  else if(flag_compute_index==SORT_COMPUTE_INDEX)
    (*index)=idx_sort;

  SID_free((void **)&rra);
}
