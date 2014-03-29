#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <gbpLib.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

void init_treenode_info_data(treenode_list_info  *list,
                             void               **rval,
                             SID_Datatype         data_type,
                             const char          *name,
                             ...){
  va_list  vargs;
  va_start(vargs,name);

  // Determine allocation size
  int dtype_size;
  SID_Type_size(data_type,&dtype_size);
  size_t alloc_size=dtype_size*list->n_list_alloc;

  // Create the new item (and apply some defaults)
  ADaPS                         *new_item;
  new_item                      =(ADaPS *)SID_malloc(sizeof(ADaPS));
  new_item->mode                =ADaPS_CUSTOM;
  new_item->free_function       =NULL;
  new_item->free_function_params=NULL;
  new_item->data                =SID_malloc(alloc_size);
  new_item->data_size           =alloc_size;

  // Give the new item its name
  vsprintf(new_item->name,name,vargs);

  // Remove any previous entries with this name
  ADaPS_remove(&(list->data),new_item->name);

  // Place new item at the start of the list
  new_item->next=list->data;
  list->data   =new_item;

  (*rval)=new_item->data;

  va_end(vargs);
}

