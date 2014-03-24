#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpTrees_build.h>

void read_a_list(const char  *filename_root_in,
                 double     **a_list_in,
                 int         *n_a_list){
   // Create snapshot expansion factor list
   FILE   *fp_in;
   char    filename_snap_list_in[MAX_FILENAME_LENGTH];
   sprintf(filename_snap_list_in,"%s/a_list.txt",filename_root_in);
   SID_log("Reading snapshot list...",SID_LOG_OPEN);
   if(SID.I_am_Master){
      // Read the list
      fp_in      =fopen(filename_snap_list_in,"r");
      (*n_a_list)=count_lines_data(fp_in);
   }
   SID_Bcast(n_a_list,sizeof(int),MASTER_RANK,SID.COMM_WORLD);
   (*a_list_in)=(double *)SID_malloc(sizeof(double)*(*n_a_list));
   if(SID.I_am_Master){
      char   *line=NULL;
      size_t  line_length=0;
      for(int i_read=0;i_read<(*n_a_list);i_read++){
         grab_next_line_data(fp_in,&line,&line_length);
         grab_double(line,1,&((*a_list_in)[i_read]));
      }
      fclose(fp_in);
   }
   SID_Bcast(a_list_in,sizeof(double)*(*n_a_list),MASTER_RANK,SID.COMM_WORLD);
   SID_log("Done.",SID_LOG_CLOSE);
}

