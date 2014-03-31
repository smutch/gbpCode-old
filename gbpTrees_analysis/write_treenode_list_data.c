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

void write_treenode_list_data(tree_info *trees,treenode_list_info *list){
  tree_node_info **list_in         =list->list;
  int              n_list_in       =list->n_list;
  int              flag_groups_list=list->flag_groups_list;
  char            *catalog_name    =list->catalog_name;

  // Open file
  char  filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_out,"%s_data.txt",list->catalog_name);
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
        if(SID.I_am_Master)
           n_list_i=n_list_in;
        else
           SID_Sendrecv(&n_list_in,
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
        // Write properties
        for(int i_list=0;i_list<n_list_i;i_list++,j_list++){
           int n_write;
           if(SID.I_am_Master)
              fprintf(fp_props_out,"%4d",j_list);
           for(int i_write=0;i_write<n_data;i_write++){
              // Perform data exchange (if needed)
              if(typ_data[i_write]==SID_INT){
                 int data_i;
                 if(i_rank!=MASTER_RANK)
                    SID_Sendrecv(&data_i,1,SID_INT,MASTER_RANK,1978271,&(((int *)(ptr_data[i_write]))[i_list]),1,SID_INT,i_rank,1978271,SID.COMM_WORLD);
                 else
                    data_i=((int *)(ptr_data[i_write]))[i_list];
                 if(SID.I_am_Master)
                    fprintf(fp_props_out," %4d",data_i);
              }
              else if(typ_data[i_write]==SID_DOUBLE){
                 double data_d;
                 if(i_rank!=MASTER_RANK)
                    SID_Sendrecv(&data_d,1,SID_DOUBLE,MASTER_RANK,1978272,&(((double *)(ptr_data[i_write]))[i_list]),1,SID_DOUBLE,i_rank,1978272,SID.COMM_WORLD);
                 else
                    data_d=((double *)(ptr_data[i_write]))[i_list];
                 if(SID.I_am_Master)
                    fprintf(fp_props_out," %le",data_d);
              }
              else
                 SID_trap_error("Unsupported data type in write_treenode_list_data() (2).",ERROR_LOGIC);
           } // i_write
           if(SID.I_am_Master) fprintf(fp_props_out,"\n");
        }
     } // if i_rank
     SID_Barrier(SID.COMM_WORLD);
  } // for i_rank
  fclose(fp_props_out);
  SID_free(SID_FARG ptr_data);
  SID_free(SID_FARG typ_data);
  SID_log("Done.",ERROR_LOGIC);
}

