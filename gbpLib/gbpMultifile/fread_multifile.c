#include <string.h>
#include <gbpMultifile.h>

#define _FILE_OFFSET_BITS 64

int fread_multifile(fp_multifile_info *fp_in,void *data_out,int item_index){
  int n_skip;
  int r_val=0;

  if(item_index<0||item_index>=fp_in->n_items_total)
     SID_trap_error("item_index (%d) is out of range (0->%d) in fread_multifile_file().",ERROR_LOGIC,item_index,fp_in->n_items_total-1);

  // Skip to the right place (if need-be)
  if(item_index!=fp_in->i_item || item_index>fp_in->i_item_stop || item_index<fp_in->i_item_start){
     // We always have to scan forward, so if we're going backwards, we have to start from scratch
     if(item_index<fp_in->i_item)         fopen_multifile_nth_file(fp_in,0);
     while(item_index>fp_in->i_item_stop) fopen_multifile_nth_file(fp_in,fp_in->i_file+1);
     n_skip=item_index-fp_in->i_item;
     if(n_skip>0)
        fseeko(fp_in->fp_multifile,(off_t)(fp_in->data_size*n_skip),SEEK_CUR);
     else if(n_skip<0)
        SID_trap_error("Negative skips (%d) not supported in fread_multifile_file().",ERROR_LOGIC,n_skip);
     fp_in->i_item+=n_skip;
  }

  // Read data
  fread_verify(data_out,fp_in->data_size,1,fp_in->fp_multifile);

  // Set counter
  fp_in->i_item++;

  return(r_val);
}

int fread_multifile_raw(fp_multifile_info *fp_in,void *data_out,int item_index){
  int n_skip;
  int r_val=0;

  if(item_index<0||item_index>=fp_in->n_items_total)
     SID_trap_error("item_index (%d) is out of range (0->%d) in fread_multifile_file().",ERROR_LOGIC,item_index,fp_in->n_items_total-1);

  // Skip to the right place (if need-be)
  if(item_index!=fp_in->i_item || item_index>fp_in->i_item_stop || item_index<fp_in->i_item_start){
     // We always have to scan forward, so if we're going backwards, we have to start from scratch
     if(item_index<fp_in->i_item_start)   fopen_multifile_nth_file(fp_in,0);
     while(item_index>fp_in->i_item_stop) fopen_multifile_nth_file(fp_in,fp_in->i_file+1);
     n_skip=item_index-fp_in->i_item;
     if(n_skip>0)
        fseeko(fp_in->fp_multifile,(off_t)(fp_in->data_size*n_skip),SEEK_CUR);
     fp_in->i_item+=n_skip;
  }

  // Read data
  fread_verify(data_out,fp_in->data_size,1,fp_in->fp_multifile);

  // Set counter
  fp_in->i_item++;

  return(r_val);
}

