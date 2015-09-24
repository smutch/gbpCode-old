#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

void write_treenode_list_data(tree_info *trees,const char *filename_out_root,treenode_list_info *list){
  tree_node_info **list_in         =list->list;
  int              n_list_local       =list->n_list_local;
  int              flag_groups_list=list->flag_groups_list;
  char            *catalog_name    =list->catalog_name;

  // Open file
  char  filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_out,"%s_%s_data.txt",filename_out_root,list->catalog_name);
  FILE *fp_props_out=fopen(filename_out,"w");

  SID_log("Writing treenode_list data to {%s}...",SID_LOG_OPEN,filename_out);

  // Write the header
  write_treenode_list_data_header(trees,list,fp_props_out);

  // Count the number of data arrays
  int    n_data=0;
  ADaPS *current=list->data;
  while(current!=NULL){
     n_data++;
     current=current->next;
  }

  // Create an array of pointers for each data array
  SID_Datatype  *typ_data=(SID_Datatype  *)SID_malloc(sizeof(SID_Datatype)*n_data);
  void         **ptr_data=(void         **)SID_malloc(sizeof(void       *)*n_data);
  int            i_data  =0;
  current=list->data;
  while(current!=NULL){
     if(current->data_type!=SID_INT && current->data_type!=SID_DOUBLE)
        SID_trap_error("Unsupported data type in write_treenode_list_data() (1).",ERROR_LOGIC);
     typ_data[i_data]=current->data_type;
     ptr_data[i_data]=current->data;
     i_data++;
     current=current->next;
  }

  // Perform rank-ordered write
  int j_list=0;
  for(int i_rank=0;i_rank<SID.n_proc;i_rank++){
     // Master Rank does all the writing
     if(SID.My_rank==i_rank || SID.I_am_Master){
        int   n_list_i;
        void *data_i;
        // Generate properties
        n_list_i=n_list_local;
        if(i_rank!=MASTER_RANK)
           SID_Sendrecv(&n_list_local,
                        1,
                        SID_INT,
                        MASTER_RANK,
                        1918270,
                        &n_list_i,
                        1,
                        SID_INT,
                        i_rank,
                        1918270,
                        SID.COMM_WORLD);
        // Allocate buffers
        void   **buffer      =(void  **)SID_malloc(sizeof(void *)*n_data);
        size_t  *data_written=(size_t *)SID_calloc(sizeof(size_t)*n_data);
        int     *dtype_size  =(int    *)SID_malloc(sizeof(int)   *n_data);
        int      n_buffer_max=4*1024;
        int      n_remaining =n_list_i;
        for(int i_data=0;i_data<n_data;i_data++)
           SID_Type_size(typ_data[i_data],&(dtype_size[i_data]));
        if(i_rank!=MASTER_RANK){
           for(int i_data=0;i_data<n_data;i_data++)
              buffer[i_data]=SID_malloc(n_buffer_max*dtype_size[i_data]);
        }
        else{
           for(int i_data=0;i_data<n_data;i_data++)
             buffer[i_data]=ptr_data[i_data];
           n_buffer_max=n_remaining;
        }
        // Fill buffers and perform write in chunks
        while(n_remaining>0){
           int n_buffer_i=MAX(n_buffer_max,n_remaining);
           // Buffer exchange between ranks
           if(i_rank!=MASTER_RANK){
              for(int i_data=0;i_data<n_data;i_data++){
                 int buffer_size=n_buffer_i*dtype_size[i_data];
                 SID_Sendrecv(&(((char *)(ptr_data[i_data]))[data_written[i_data]]),buffer_size,SID_CHAR,MASTER_RANK,1978271,buffer[i_data],buffer_size,SID_CHAR,i_rank,1978271,SID.COMM_WORLD);
              }
           }
           // Write buffer
           if(SID.I_am_Master){
              for(int i_buffer=0;i_buffer<n_buffer_i;i_buffer++,j_list++){
                 fprintf(fp_props_out,"%4d",j_list);
                 for(int i_data=0;i_data<n_data;i_data++){
                    if     (typ_data[i_data]==SID_INT)    fprintf(fp_props_out," %4d",    ((int    *)(buffer[i_data]))[i_buffer]);
                    else if(typ_data[i_data]==SID_DOUBLE) fprintf(fp_props_out," %11.4le",((double *)(buffer[i_data]))[i_buffer]);
                    else SID_trap_error("Unsupported data type in write_treenode_list_data() (2).",ERROR_LOGIC);
                 }
                 fprintf(fp_props_out,"\n");
              }
           }
           // Update counters
           n_remaining-=n_buffer_i;
           for(int i_data=0;i_data<n_data;i_data++)
              data_written[i_data]+=(size_t)n_buffer_i*(size_t)dtype_size[i_data];
        }
        // Free buffers
        if(i_rank!=MASTER_RANK){
           for(int i_data=0;i_data<n_data;i_data++)
              SID_free(SID_FARG buffer[i_data]);
        }
        SID_free(SID_FARG buffer);
        SID_free(SID_FARG data_written);
        SID_free(SID_FARG dtype_size);
     } // if i_rank
     SID_Barrier(SID.COMM_WORLD);
  } // for i_rank
  fclose(fp_props_out);
  SID_free(SID_FARG ptr_data);
  SID_free(SID_FARG typ_data);
  SID_log("Done.",ERROR_LOGIC);
}

