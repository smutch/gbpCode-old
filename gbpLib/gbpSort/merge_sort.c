#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpSort.h>

// This sort routine is RAM-inefficient, but it is STABLE ... used
//  when computing a distributed rank-list
void merge_sort(void    *data_in,
                size_t   n_data,
                size_t **index,
                SID_Datatype  data_type,
                int      flag_compute_index,
                int      flag_in_place){
  size_t  i,j;
  void   *data_sort;
  void   *scratch_d;
  size_t *scratch_i;
  size_t  data_size;
  size_t *idx_sort;
  int     flag_compute_index_sort;

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
  else if(data_type==SID_LONG_LONG)
    data_size=sizeof(long long);
  else
    SID_trap_error("Unsupported data type {%d}",ERROR_LOGIC,data_type);

  // Allocate scratch space and sort indices (if needed)
  scratch_d=SID_malloc(data_size*n_data);
  if(flag_in_place==SORT_COMPUTE_INPLACE)
    data_sort=data_in;
  else if(flag_in_place==SORT_COMPUTE_NOT_INPLACE){
    data_sort=SID_malloc(data_size*n_data);
    memcpy(data_sort,data_in,data_size*n_data);
  }
  else{
    fprintf(stderr,"ERROR (merge_sort): flag_in_place {%d} must be SORT_COMPUTE_INPLACE||SORT_COMPUTE_NOT_INPLACE\n",
	    flag_in_place);
  }
  flag_compute_index_sort=FALSE;
  if(flag_compute_index==SORT_COMPUTE_INDEX || 
     flag_compute_index==SORT_COMPUTE_RANK){
    flag_compute_index_sort=TRUE;
    scratch_i=(size_t *)SID_malloc(sizeof(size_t)*n_data);
    idx_sort =(size_t *)SID_malloc(sizeof(size_t)*n_data);
    for(i=0;i<n_data;i++)
      idx_sort[i]=i;
  }
  else if(flag_compute_index!=SORT_INPLACE_ONLY){
    fprintf(stderr,"ERROR (merge_sort): flag_compute_index {%d} must be SORT_COMPUTE_INDEX||SORT_COMPUTE_RANK||SORT_INPLACE_ONLY\n",
	    flag_compute_index);
  }

  // Check for nonsensical flag combinations
  if(flag_in_place==SORT_COMPUTE_NOT_INPLACE && flag_compute_index==SORT_INPLACE_ONLY)
    fprintf(stderr,"ERROR (mege_sort): flag combinations will generate no results.\n");

  // Perform the sort
  if(n_data>0)
    merge_helper(data_sort,
                 idx_sort,
                 0,
                 n_data,
                 scratch_d,
                 scratch_i,
                 data_type,
                 flag_compute_index_sort);
  SID_free((void **)&(scratch_d));
  if(flag_compute_index_sort)
    SID_free((void **)&scratch_i);
  if(flag_in_place==SORT_COMPUTE_NOT_INPLACE)
    SID_free((void **)&data_sort);

  // Invert sort indices to produce ranks (if needed)
  if(flag_compute_index==SORT_COMPUTE_RANK && n_data>0){
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
}

// This is the main recursive routine for the merge sort
void merge_helper(void   *data,
                  size_t *index,
                  size_t  left,
                  size_t  right,
                  void   *scratch_d,
                  size_t *scratch_i,
                  SID_Datatype    data_type,
                  int     flag_compute_index){
  // left is the index of the leftmost element of the subarray; right is one
  //   past the index of the rightmost element 
  size_t length=right-left;
  size_t midpoint_distance=length/2;
  size_t l = left, r = left + midpoint_distance;
  size_t i=0;
  // base case: one element
  if(right==left+1)
    return;
  else{
    // sort each (left and right) subarray
    merge_helper(data,
                 index,
                 left,
                 left+midpoint_distance,
                 scratch_d,
                 scratch_i,
                 data_type,
                 flag_compute_index);
    merge_helper(data,
                 index,
                 left+midpoint_distance,
                 right,
                 scratch_d,
                 scratch_i,
                 data_type,
                 flag_compute_index);

    // merge the arrays together using scratch for temporary storage
    for(i=0;i<length;i++){
      // Check to see if any elements remain in the left array; if so,
      //   we check if there are any elements left in the right array; if
      //   so, we compare them.  Otherwise, we know that the merge must
      //   use the element from the left array.  Then, we copy the
      //   sorted subarray back into the input data.
      if(data_type==SID_INT){
        if(l<left+midpoint_distance && (r==right||(((int *)data)[l]<=((int *)data)[r]))){
          ((int *)scratch_d)[i]=((int *)data)[l];
          if(flag_compute_index)
            scratch_i[i]=index[l];
          l++;
        }
        else{
          ((int *)scratch_d)[i]=((int *)data)[r];
          if(flag_compute_index)
            scratch_i[i]=index[r];
          r++;
        }
      }
      else if(data_type==SID_SIZE_T){
        if(l<left+midpoint_distance && (r==right||(((size_t *)data)[l]<=((size_t *)data)[r]))){
          ((size_t *)scratch_d)[i]=((size_t *)data)[l];
          if(flag_compute_index)
            scratch_i[i]=index[l];
          l++;
        }
        else{
          ((size_t *)scratch_d)[i]=((size_t *)data)[r];
          if(flag_compute_index)
            scratch_i[i]=index[r];
          r++;
        }
      }
      else if(data_type==SID_FLOAT){
        if(l<left+midpoint_distance && (r==right||(((float *)data)[l]<=((float *)data)[r]))){
          ((float *)scratch_d)[i]=((float *)data)[l];
          if(flag_compute_index)
            scratch_i[i]=index[l];
          l++;
        }
        else{
          ((float *)scratch_d)[i]=((float *)data)[r];
          if(flag_compute_index)
            scratch_i[i]=index[r];
          r++;
        }
      }
      else if(data_type==SID_DOUBLE){
        if(l<left+midpoint_distance && (r==right||(((double *)data)[l]<=((double *)data)[r]))){
          ((double *)scratch_d)[i]=((double *)data)[l];
          if(flag_compute_index)
            scratch_i[i]=index[l];
          l++;
        }
        else{
          ((double *)scratch_d)[i]=((double *)data)[r];
          if(flag_compute_index)
            scratch_i[i]=index[r];
          r++;
        }
      }
    }
  }
  if(data_type==SID_INT){
    for(i=left;i<right;i++){
      ((int *)data)[i]=((int *)scratch_d)[i-left];
      if(flag_compute_index)
	index[i]=scratch_i[i-left];
    }
  }
  else if(data_type==SID_SIZE_T){
    for(i=left;i<right;i++){
      ((size_t *)data)[i]=((size_t *)scratch_d)[i-left];
      if(flag_compute_index)
	index[i]=scratch_i[i-left];
    }
  }
  else if(data_type==SID_FLOAT){
    for(i=left;i<right;i++){
      ((float *)data)[i]=((float *)scratch_d)[i-left];
      if(flag_compute_index)
	index[i]=scratch_i[i-left];
    }
  }
  else if(data_type==SID_DOUBLE){
    for(i=left;i<right;i++){
      ((double *)data)[i]=((double *)scratch_d)[i-left];
      if(flag_compute_index)
	index[i]=scratch_i[i-left];
    }
  }
}
