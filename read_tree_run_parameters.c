#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void read_tree_run_parameters(char *filename_root_in,
                              int  *i_read_start,
                              int  *i_read_stop,
                              int  *i_read_step,
                              int  *n_search,
                              int  *flag_fix_bridges,
                              int  *flag_compute_fragmented,
                              int  *flag_compute_ghosts){
   // Open file
   FILE *fp_in;
   char  filename_in[MAX_FILENAME_LENGTH];
   sprintf(filename_in,"%s/run_parameters.txt",filename_root_in);
   if((fp_in=fopen(filename_in,"r"))==NULL)
      SID_trap_error("Could not open file {%s}",ERROR_IO_OPEN,filename_in);

   // Read file
   char   *line=NULL;
   size_t  line_length=0;
   char    variable[64];
   int     value;
   int     i_variable;
   for(i_variable=0;i_variable<7;i_variable++){
      grab_next_line(fp_in,&line,&line_length);
      grab_word(line,1,variable);
      grab_int( line,2,&value);
      switch(i_variable){
         case 0:
            if(!strcmp(variable,"#snapshot_start"))
               (*i_read_start)=value;
            else
               SID_trap_error("Invalid run_parameters.txt file.",ERROR_LOGIC);
            break;
         case 1:
            if(!strcmp(variable,"#snapshot_stop"))
               (*i_read_stop)=value;
            else
               SID_trap_error("Invalid run_parameters.txt file.",ERROR_LOGIC);
            break;
         case 2:
            if(!strcmp(variable,"#snapshot_step"))
               (*i_read_step)=value;
            else
               SID_trap_error("Invalid run_parameters.txt file.",ERROR_LOGIC);
            break;
         case 3:
            if(!strcmp(variable,"#tree_scan_range"))
               (*n_search)=value;
            else
               SID_trap_error("Invalid run_parameters.txt file.",ERROR_LOGIC);
            break;
         case 4:
            if(!strcmp(variable,"#flag_fix_bridges"))
               (*flag_fix_bridges)=value;
            else
               SID_trap_error("Invalid run_parameters.txt file.",ERROR_LOGIC);
            break;
         case 5:
            if(!strcmp(variable,"#flag_compute_fragmented"))
               (*flag_compute_fragmented)=value;
            else
               SID_trap_error("Invalid run_parameters.txt file.",ERROR_LOGIC);
            break;
         case 6:
            if(!strcmp(variable,"#flag_compute_ghosts"))
               (*flag_compute_ghosts)=value;
            else
               SID_trap_error("Invalid run_parameters.txt file.",ERROR_LOGIC);
            break;
      }
   }

   // Close file
   fclose(fp_in);
}

