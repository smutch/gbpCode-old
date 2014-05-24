#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpTrees_build.h>

void write_a_list(const char *filename_snap_list_in,
                  const char *filename_root_out,
                  int         i_read_start,
                  int         i_read_stop,
                  int         i_read_step){

   // Create snapshot expansion factor list
   if(SID.I_am_Master){
      char filename_snap_list_out[MAX_FILENAME_LENGTH];
      sprintf(filename_snap_list_out,"%s/a_list.txt",filename_root_out);
      SID_log("Creating snapshot list...",SID_LOG_OPEN);
  
      // Read the original list
      char   *line=NULL;
      size_t  line_length=0;
      FILE   *fp_in;
      int     i_read;
      int     n_snaps;
      double *a_list_in;
      double *a_list_out;
      fp_in     =fopen(filename_snap_list_in, "r");
      n_snaps   =count_lines_data(fp_in);
      a_list_in =(double *)SID_malloc(sizeof(double)*n_snaps);
      a_list_out=(double *)SID_malloc(sizeof(double)*n_snaps);
      for(i_read=0;i_read<n_snaps;i_read++){
         grab_next_line_data(fp_in,&line,&line_length);
         grab_double(line,1,&(a_list_in[i_read]));
      }
      fclose(fp_in);
  
      // Select the snapshots we've used
      int i_next;
      int n_keep;
      for(i_read=i_read_stop,i_next=i_read_stop,n_keep=0;i_read>=i_read_start;i_read--){
         if(i_read==i_next){
            a_list_out[n_keep++]=a_list_in[i_read];
            i_next-=i_read_step;
         }
      }
  
      // Write them to a file
      int   i_write;
      FILE *fp_out;
      fp_out=fopen(filename_snap_list_out,"w");
      for(i_write=n_keep-1;i_write>=0;i_write--)
         fprintf(fp_out,"%le\n",a_list_out[i_write]);
      fclose(fp_out);
 
      SID_free(SID_FARG line); 
      SID_free(SID_FARG a_list_in);
      SID_free(SID_FARG a_list_out);
  
      SID_log("Done.",SID_LOG_CLOSE);
   }
   SID_Barrier(SID.COMM_WORLD);
}

