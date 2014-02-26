#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpADaPS.h>

void ADaPS_store(ADaPS      **list,
                 void        *data,
                 const char  *name,
                 int          mode, ...){
  ADaPS   *new_item;
  size_t   n_subarray;
  size_t   data_size;
  va_list  vargs;
  va_start(vargs,mode);

  // Create the new item (and apply some defaults)
  new_item               =(ADaPS *)SID_malloc(sizeof(ADaPS));
  new_item->data         =NULL;
  new_item->mode         =mode;
  new_item->free_function=NULL;

  // Add data to the new item and interpret variable arguments
  if(check_mode_for_flag(mode,ADaPS_COPY)){
    data_size     =(size_t)va_arg(vargs,size_t);
    new_item->data=(void *)SID_malloc(data_size);
    memcpy(new_item->data,data,data_size);
  }
  else if(check_mode_for_flag(mode,ADaPS_COPY_SUBARRAY_DOUBLE)){
    n_subarray    =(size_t)va_arg(vargs,size_t);
    data_size     =sizeof(double)*n_subarray;
    new_item->data=(void *)SID_malloc(data_size);
    memcpy(new_item->data,data,data_size);
  }
  else if(check_mode_for_flag(mode,ADaPS_COPY_SUBARRAY_FLOAT)){
    n_subarray    =(size_t)va_arg(vargs,size_t);
    data_size     =sizeof(float)*n_subarray;
    new_item->data=(void *)SID_malloc(data_size);
    memcpy(new_item->data,data,data_size);
  }
  else if(check_mode_for_flag(mode,ADaPS_COPY_SUBARRAY_REAL)){
    n_subarray    =(size_t)va_arg(vargs,size_t);
    data_size     =sizeof(GBPREAL)*n_subarray;
    new_item->data=(void *)SID_malloc(data_size);
    memcpy(new_item->data,data,data_size);
  }
  else if(check_mode_for_flag(mode,ADaPS_COPY_SUBARRAY_SIZE_T)){
    n_subarray    =(size_t)va_arg(vargs,size_t);
    data_size     =sizeof(size_t)*n_subarray;
    new_item->data=(void *)SID_malloc(data_size);
    memcpy(new_item->data,data,data_size);
  }
  else if(check_mode_for_flag(mode,ADaPS_COPY_SUBARRAY_INT)){
    n_subarray    =(size_t)va_arg(vargs,size_t);
    data_size     =sizeof(int)*n_subarray;
    new_item->data=(void *)SID_malloc(data_size);
    memcpy(new_item->data,data,data_size);
  }
  else if(check_mode_for_flag(mode,ADaPS_SCALAR_DOUBLE)){
    data_size     =sizeof(double);
    new_item->data=(void *)SID_malloc(data_size);
    memcpy(new_item->data,data,data_size);
  }
  else if(check_mode_for_flag(mode,ADaPS_SCALAR_FLOAT)){
    data_size     =sizeof(float);
    new_item->data=(void *)SID_malloc(data_size);
    memcpy(new_item->data,data,data_size);
  }
  else if(check_mode_for_flag(mode,ADaPS_SCALAR_REAL)){
    data_size     =sizeof(GBPREAL);
    new_item->data=(void *)SID_malloc(data_size);
    memcpy(new_item->data,data,data_size);
  }
  else if(check_mode_for_flag(mode,ADaPS_SCALAR_SIZE_T)){
    data_size     =sizeof(size_t);
    new_item->data=(void *)SID_malloc(data_size);
    memcpy(new_item->data,data,data_size);
  }
  else if(check_mode_for_flag(mode,ADaPS_SCALAR_INT)){
    data_size     =sizeof(int);
    new_item->data=(void *)SID_malloc(data_size);
    memcpy(new_item->data,data,data_size);
  }
  else 
    new_item->data=data;

  // Give the new item its name
  vsprintf(new_item->name,name,vargs);

  // Remove any previous entries with this name
  ADaPS_remove(&(*list),new_item->name);

  // Place new item at the start of the list
  new_item->next=(*list);
  (*list)       =new_item;

  va_end(vargs);
}

