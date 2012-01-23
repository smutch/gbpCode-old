#include <stdio.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpFITS.h>

void array_size_local(int n_d,int *d_i,size_t *n_data);
void array_size_local(int n_d,int *d_i,size_t *n_data){
  int i_d;
  if(n_d==0)
      (*n_data)=0;
  else{
    for(i_d=0,(*n_data)=1;i_d<n_d;i_d++)
        (*n_data)*=(size_t)(d_i[i_d]);
  }
}

void indices2index_local(int *indices,int n_d,int *d_i,size_t *index);
void indices2index_local(int *indices,int n_d,int *d_i,size_t *index){
  int i_d;
  for(i_d=0,(*index)=0;i_d<n_d;i_d++){
      (*index)*=(size_t)(d_i[i_d]);
      (*index)+=(size_t)(indices[i_d]);
  }
}

void index2indices_local(size_t index,int n_d,int *d_i,int *indices);
void index2indices_local(size_t index,int n_d,int *d_i,int *indices){
  int    i_d;
  size_t remainder;
  for(i_d=n_d-1,remainder=index;i_d>=0;i_d--){
    indices[i_d]=(int)(remainder%((size_t)(d_i[i_d])));
    remainder-=((size_t)indices[i_d]);
    remainder/=((size_t)d_i[i_d]);
  }
}

void transpose_array(void *data,SID_Datatype dtype,int n_d,int *d_i);
void transpose_array(void *data,SID_Datatype dtype,int n_d,int *d_i){
    void   *buffer;
    int     type_size;
    size_t  n_data;
    size_t  index;
    size_t  index_transpose;
    int    *indices;
    int    *indices_transpose;
    int    *d_transpose_i;
    int     i_d;

    // Get the element size for the data's given datatype
    SID_Type_size(dtype,&type_size);

    // Compute the number of elements in data
    array_size_local(n_d,d_i,&n_data);

    // Copy the data into a buffer
    buffer=(void *)SID_malloc((size_t)type_size*n_data);
    memcpy(buffer,data,(size_t)type_size*n_data);

    // Create the dimensions of the transposed array
    d_transpose_i=(int *)SID_malloc(sizeof(int)*n_d);
    for(i_d=0;i_d<n_d;i_d++)
        d_transpose_i[n_d-i_d-1]=d_i[i_d];

    // Perform the transpose
    indices          =(int *)SID_malloc(sizeof(int)*n_d);
    indices_transpose=(int *)SID_malloc(sizeof(int)*n_d);
    for(index=0;index<n_data;index++){
        // Convert the original element's index to its transposed index
        index2indices_local(index,n_d,d_i,indices);
        for(i_d=0;i_d<n_d;i_d++)
            indices_transpose[n_d-i_d-1]=indices[i_d];
        indices2index_local(indices_transpose,n_d,d_transpose_i,&index_transpose);
        // Copy over the data from the buffer
        if(dtype==SID_DOUBLE)
          ((double *)data)[index_transpose]=((double *)buffer)[index];        
        else if(dtype==SID_FLOAT){
          ((float *)data)[index_transpose]=((float *)buffer)[index];        
        }
        else if(dtype==SID_INT)
          ((int *)data)[index_transpose]=((int *)buffer)[index];
        else
          SID_trap_error("Unrecognized dtype in transpose_array",ERROR_LOGIC);
    }

    // Transpose the given array with the dimension sizes
    memcpy(d_i,d_transpose_i,sizeof(int)*n_d);
    
    // Clean-up
    SID_free(SID_FARG buffer);
    SID_free(SID_FARG d_transpose_i);
    SID_free(SID_FARG indices);
    SID_free(SID_FARG indices_transpose);
}

