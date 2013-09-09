#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void write_tree_run_parameters(char *filename_root_out,
                               int   i_read_start,
                               int   i_read_stop,
                               int   i_read_step,
                               int   n_search,
                               int   flag_fix_bridges,
                               int   flag_compute_fragmented,
                               int   flag_compute_ghosts){
   // Open file
   FILE *fp_out;
   char  filename_out[MAX_FILENAME_LENGTH];
   mkdir(filename_root_out,02755);
   sprintf(filename_out,"%s/run_parameters.txt",filename_root_out);
   if((fp_out=fopen(filename_out,"w"))==NULL)
      SID_trap_error("Could not open file {%s}",ERROR_IO_OPEN,filename_root_out);

   // Write file
   fprintf(fp_out,"#snapshot_start          %d\n",i_read_start);
   fprintf(fp_out,"#snapshot_stop           %d\n",i_read_stop);
   fprintf(fp_out,"#snapshot_step           %d\n",i_read_step);
   fprintf(fp_out,"#tree_scan_range         %d\n",n_search);
   fprintf(fp_out,"#flag_fix_bridges        %d\n",flag_fix_bridges);
   fprintf(fp_out,"#flag_compute_fragmented %d\n",flag_compute_fragmented);
   fprintf(fp_out,"#flag_compute_ghosts     %d\n",flag_compute_ghosts);

   // Close file
   fclose(fp_out);
}

