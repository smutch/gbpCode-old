#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_core.h>

#define MAX_ITEM_SIZE_LOCAL 16

void swap_endian(char *string,int n_items,int item_byte_size){
   if(item_byte_size>MAX_ITEM_SIZE_LOCAL)
      SID_trap_error("An attempt has been made to swap endian of an item >%d bytes in size.  Increase the limit and recompile swap_endian.c",ERROR_LOGIC,MAX_ITEM_SIZE_LOCAL);
   char item_copy[MAX_ITEM_SIZE_LOCAL];
   for(int i_item=0;i_item<n_items;i_item++){
      memcpy(&item_copy[0],&string[i_item*item_byte_size],item_byte_size);
      for(int i_byte=0;i_byte<item_byte_size;i_byte++)
         string[i_item*item_byte_size+i_byte]=item_copy[item_byte_size-i_byte-1];
   }
}

