#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees.h>

void write_trees_vertical(tree_info **trees,
                          int        *n_halos_tree_local,
                          int         n_trees_local,
                          int        *tree_lo_file,
                          int        *tree_hi_file,
                          int        *n_halos_file,
                          int         n_files,
                          char       *filename_root_out,
                          char       *group_text_prefix){
  int              i_rank;
  int              i_write;
  int              i_tree,j_tree;
  int              flag_write_init;
  int              n_trees_file;
  int              n_halos;
  int              n_trees;
  int              n_halos_written;
  int              n_halos_tree_written;
  int              n_trees_written;
  tree_node_info  *current;
  tree_node_info  *next;
  char             filename_out[MAX_FILENAME_LENGTH];
  SID_fp           fp_out;

  // Write headers first
  calc_sum_global(n_halos_tree_local,&n_halos,n_trees_local,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
  calc_sum_global(&n_trees_local,    &n_trees,1,            SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
  SID_log("Writing %d %sgroup trees (structure size=%lld bytes)...",SID_LOG_OPEN|SID_LOG_TIMER,n_trees,group_text_prefix,sizeof(halo_info));
  SID_log("Writing headers...",SID_LOG_OPEN);
  for(i_rank=0,i_write=0,i_tree=0,n_halos_written=0,n_trees_written=0,flag_write_init=TRUE;i_rank<SID.n_proc;i_rank++){
    if(SID.My_rank==i_rank){
      sprintf(filename_out,"%s.%sgroup_trees.%d",filename_root_out,group_text_prefix,i_write);
      if(flag_write_init){
        SID_fopen(filename_out,"w",&fp_out);
        n_trees_file=tree_hi_file[i_write]-tree_lo_file[i_write]+1;
        SID_fwrite(&n_trees_file,           sizeof(int),1,&fp_out);
        SID_fwrite(&(n_halos_file[i_write]),sizeof(int),1,&fp_out);
      }
      else
        SID_fopen(filename_out,"a",&fp_out);
      for(j_tree=0;j_tree<n_trees_local;i_tree++,j_tree++){
        SID_fwrite(&(n_halos_tree_local[i_tree]),sizeof(int),1,&fp_out);
        n_halos_written+=n_halos_tree_local[i_tree];
        n_trees_written++;
        if(i_tree==tree_hi_file[i_write]){
          i_write++;
          if(i_write<n_files){
            flag_write_init=TRUE;
            SID_fclose(&fp_out);
            sprintf(filename_out,  "%s.%sgroup_trees.%d",filename_root_out,group_text_prefix,i_write);
            SID_fopen(filename_out,"w",&fp_out);
            n_trees_file=tree_hi_file[i_write]-tree_lo_file[i_write]+1;
            SID_fwrite(&n_trees_file,           sizeof(int),1,&fp_out);
            SID_fwrite(&(n_halos_file[i_write]),sizeof(int),1,&fp_out);
          }
          else if(n_halos_written!=n_halos || n_trees_written!=n_trees)
            SID_trap_error("Not all halos (%d=?%d) and/or trees (%d=?%d) have been written.",ERROR_LOGIC,n_halos_written,n_halos,n_trees_written,n_trees);
        }
      }
      SID_fclose(&fp_out);
    }
    // Update the other ranks 
    SID_Bcast(&i_write,        sizeof(int),i_rank,SID.COMM_WORLD);
    SID_Bcast(&i_tree,         sizeof(int),i_rank,SID.COMM_WORLD);
    SID_Bcast(&n_halos_written,sizeof(int),i_rank,SID.COMM_WORLD);
    SID_Bcast(&n_trees_written,sizeof(int),i_rank,SID.COMM_WORLD);
    SID_Bcast(&flag_write_init,sizeof(int),i_rank,SID.COMM_WORLD);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Writing %d halos...",SID_LOG_OPEN,n_halos);
  for(i_rank=0,i_write=0,i_tree=0,n_halos_written=0,n_trees_written=0,flag_write_init=TRUE;i_rank<SID.n_proc;i_rank++){
    if(SID.My_rank==i_rank){
      sprintf(filename_out,"%s.%sgroup_trees.%d",filename_root_out,group_text_prefix,i_write);
      SID_fopen(filename_out,"a",&fp_out);
      for(j_tree=0;j_tree<n_trees_local;i_tree++,j_tree++){
        n_halos_tree_written=0;
        current=trees[i_tree]->root;
        while(current!=NULL){
          if(current->descendant==NULL)
            n_halos_tree_written+=write_tree_vertical_halos_recursive(current,&fp_out,NULL);
          current=current->next;
        }
        if(n_halos_tree_written!=n_halos_tree_local[i_tree])
          SID_trap_error("Number of halos written is not right (i.e. %d!=%d) for tree #%d",ERROR_LOGIC,n_halos_written,n_halos_tree_local[i_tree],i_tree);
        n_halos_written+=n_halos_tree_written;
        n_trees_written++;
        if(i_tree==tree_hi_file[i_write]){
          i_write++;
          if(i_write<n_files){
            flag_write_init=TRUE;
            SID_fclose(&fp_out);
            sprintf(filename_out,  "%s.%sgroup_trees.%d",filename_root_out,group_text_prefix,i_write);
            SID_fopen(filename_out,"a",&fp_out);
            n_trees_file=tree_hi_file[i_write]-tree_lo_file[i_write]+1;
          }
          else if(n_halos_written!=n_halos || n_trees_written!=n_trees)
            SID_trap_error("Not all halos (%d=?%d) and/or trees (%d=?%d) have been written.",ERROR_LOGIC,n_halos_written,n_halos,n_trees_written,n_trees);
        }
      }
      SID_fclose(&fp_out);
    }
    // Update the other ranks 
    SID_Bcast(&i_write,        sizeof(int),i_rank,SID.COMM_WORLD);
    SID_Bcast(&i_tree,         sizeof(int),i_rank,SID.COMM_WORLD);
    SID_Bcast(&n_halos_written,sizeof(int),i_rank,SID.COMM_WORLD);
    SID_Bcast(&n_trees_written,sizeof(int),i_rank,SID.COMM_WORLD);
  }
  SID_log("Done.",SID_LOG_CLOSE);
  SID_log("Halos written: %d",SID_LOG_COMMENT,n_halos_written);
  SID_log("Trees written: %d",SID_LOG_COMMENT,n_trees_written);
  SID_log("Done.",SID_LOG_CLOSE);
}

