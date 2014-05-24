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

void write_treenode_list_properties(tree_info *trees,const char *filename_out_root,treenode_list_info *list){
  tree_node_info **list_in         =list->list;
  int              n_list_in       =list->n_list;
  int              flag_groups_list=list->flag_groups_list;
  char            *catalog_name    =list->catalog_name;

  // Open file
  char  filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_out,"%s_%s_properties.txt",filename_out_root,list->catalog_name);
  FILE *fp_props_out=fopen(filename_out,"w");

  SID_log("Writing treenode_list properties to {%s}...",SID_LOG_OPEN,filename_out);

  // Write the header
  write_treenode_list_properties_header(trees,list,fp_props_out);

  // Fetch the pointer to the properties
  halo_properties_info **properties;
  if(flag_groups_list){
     if(trees->group_properties!=NULL)
        properties=trees->group_properties;
     else
        SID_trap_error("Group properties are not defined.  They probably have not been read.",ERROR_LOGIC);
  }
  else{
     if(trees->subgroup_properties!=NULL)
        properties=trees->subgroup_properties;
     else
        SID_trap_error("Subgroup properties are not defined.  They probably have not been read.",ERROR_LOGIC);
  }


  // Perform rank-ordered write
  int j_list=0;
  for(int i_rank=0;i_rank<SID.n_proc;i_rank++){
     // Master Rank does all the writing
     if(SID.My_rank==i_rank || SID.I_am_Master){
        int   n_list_i;
        void *data_i;
        // Generate properties
        if(i_rank==0)
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
           // Point to the halo to be processed 
           tree_node_info *current_halo;
           if(i_rank==SID.My_rank)
              current_halo=list_in[i_list];

           int n_write;
           if(i_rank==0)
              fprintf(fp_props_out,"%4d",j_list);
           for(int i_write=0;i_write<WRITE_TREENODE_LIST_PROPERTIES_N;i_write++){
              // Set the i_write'th data
              SID_Datatype  data_type;
              int           data_i;
              double        data_d;
              if(i_rank==SID.My_rank)
                 write_treenode_list_properties_set_ith(trees,
                                                        current_halo,
                                                        i_write,
                                                        NULL,
                                                        &data_type,
                                                        &data_i,
                                                        &data_d);
              if(i_rank!=MASTER_RANK){
                 if(data_type==SID_INT)
                    SID_Sendrecv(&data_i,1,data_type,MASTER_RANK,2178271,&data_i,1,data_type,i_rank,2178271,SID.COMM_WORLD);
                 else if(data_type==SID_DOUBLE) 
                    SID_Sendrecv(&data_d,1,data_type,MASTER_RANK,2178272,&data_d,1,data_type,i_rank,2178272,SID.COMM_WORLD);
                 else  
                    SID_trap_error("Unsupported data type in write_treenode_list_data() (2).",ERROR_LOGIC);
              }

              // Write the data with the master rank
              if(SID.I_am_Master){
                 if(data_type==SID_INT)         fprintf(fp_props_out," %4d",data_i);
                 else if(data_type==SID_DOUBLE) fprintf(fp_props_out," %le",data_d);
                 else                           SID_trap_error("Unsupported data type in write_treenode_list_data() (2).",ERROR_LOGIC);
              }
           } // i_write
           if(SID.I_am_Master) fprintf(fp_props_out,"\n");
        }
     } // if i_rank
     SID_Barrier(SID.COMM_WORLD);
  } // for i_rank
  fclose(fp_props_out);
  SID_log("Done.",ERROR_LOGIC);
}

