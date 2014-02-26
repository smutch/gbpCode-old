#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <gbpLib.h>
#include <gbpTrees.h>

void init_tree_data(tree_info    *trees,
                    void       ***rval,
                    size_t        data_size,
                    int           mode,
                    const char   *name,
                    ...){
  va_list  vargs;
  va_start(vargs,name);

  // Create the new item (and apply some defaults)
  ADaPS                           *new_item;
  store_tree_data_free_parms_info *params;
  params                        =(store_tree_data_free_parms_info *)SID_malloc(sizeof(store_tree_data_free_parms_info));
  params->n_snaps               =trees->n_snaps;
  new_item                      =(ADaPS *)SID_malloc(sizeof(ADaPS));
  new_item->mode                =ADaPS_CUSTOM;
  new_item->free_function       =free_tree_data;
  new_item->free_function_params=params;

  // Determine the size of the allocation
  new_item->data_size=trees->n_snaps*sizeof(void *);

  // Allocate arrays
  new_item->data=SID_malloc(sizeof(void *)*trees->n_snaps);
  void **data=(void **)new_item->data;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     size_t n_items;
     size_t alloc_size;
     if(mode==INIT_TREE_DATA_GROUPS)
        n_items=trees->n_groups_snap_local[i_snap];
     else if(mode==INIT_TREE_DATA_SUBGROUPS)
        n_items=trees->n_subgroups_snap_local[i_snap];
     else
        SID_trap_error("Invalid init_tree_data() mode (%d).",ERROR_LOGIC,mode);
     alloc_size          =data_size*n_items;
     data[i_snap]        =SID_malloc(alloc_size);
     new_item->data_size+=alloc_size;
  }

  // Give the new item its name
  vsprintf(new_item->name,name,vargs);

  // Remove any previous entries with this name
  ADaPS_remove(&(trees->data),new_item->name);

  // Place new item at the start of the list
  new_item->next=trees->data;
  trees->data   =new_item;

  (*rval)=data;

  va_end(vargs);
}

