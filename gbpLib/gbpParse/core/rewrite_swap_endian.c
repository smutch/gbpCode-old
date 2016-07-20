#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_core.h>

void rewrite_swap_endian(FILE *fp_in,FILE *fp_out,int n_items,int item_byte_size,char *buffer_in){
   // Create the input buffer if one has not been passed to us
   char *buffer=NULL;
   if(buffer_in==NULL)
      buffer=(char *)SID_malloc(sizeof(char)*n_items*item_byte_size);
   else
      buffer=buffer_in;

   // Populate the input buffer
   size_t n_return;
   if((n_return=fread(buffer,item_byte_size,n_items,fp_in))!=(size_t)n_items)
      SID_trap_error("Failed to read %d %d-byte items (only got %d).",ERROR_IO_READ,n_items,item_byte_size,(int)n_return);

   // Swap endian
   swap_endian(buffer,n_items,item_byte_size);

   // Write to output file
   fwrite(buffer,item_byte_size,n_items,fp_out);
   
   // Free buffer (if necessary)
   if(buffer_in==NULL)
      SID_free(SID_FARG buffer);
}
   
