#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <gbpLib.h>
#include <gbpInterpolate.h>

void ADaPS_store_interp(ADaPS        **list,
                        void          *data,
                        const char    *name,
                        ...){
  ADaPS   *new_item;
  size_t   n_subarray;
  size_t   data_size;
  va_list  vargs;
  va_start(vargs,name);

  // Create the new item (and apply some defaults)
  new_item                      =(ADaPS *)SID_malloc(sizeof(ADaPS));
  new_item->data                =data;
  new_item->mode                =ADaPS_CUSTOM;
  new_item->free_function       =free_interpolate;
  new_item->free_function_params=NULL;
  data_size                     =sizeof(interp_info);

  // Give the new item its name
  vsprintf(new_item->name,name,vargs);

  // Remove any previous entries with this name
  ADaPS_remove(&(*list),new_item->name);

  // Place new item at the start of the list
  new_item->next=(*list);
  (*list)       =new_item;

  va_end(vargs);
}

