#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpSims_Swap_Endian.h>

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  if(argc!=5) SID_trap_error("Invalid syntax.",ERROR_SYNTAX);
  char filename_SSimPL_in[MAX_FILENAME_LENGTH];
  char filename_halo_type[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_in, argv[1]);
  strcpy(filename_halo_type, argv[2]);
  int i_snap_lo=atoi(argv[3]);
  int i_snap_hi=atoi(argv[4]);

  // Set which files will be processed
  int flag_process_halos    =TRUE;
  int flag_process_catalogs =TRUE;
  int flag_process_grids    =TRUE;
  int flag_process_snapshots=TRUE;

  // Loop over the given snapshot range
  SID_log("Processing group/subgroup statistics for files #%d->#%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_snap_lo,i_snap_hi);
  for(int i_snap=i_snap_lo;i_snap<=i_snap_hi;i_snap++){
    char filename_in[MAX_FILENAME_LENGTH];
    SID_log("Processing snapshot #%03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_snap);
    // Process halos
    if(flag_process_halos){
       sprintf(filename_in, "%s/halos/",filename_SSimPL_in);
       check_integrity_halos(filename_in,filename_halo_type,i_snap);
    }
    // Process catalogs
    if(flag_process_catalogs){
       sprintf(filename_in, "%s/catalogs/",filename_SSimPL_in, filename_halo_type);
       check_integrity_catalogs(filename_in,filename_halo_type,i_snap);
    }
    // Process grids
    if(flag_process_grids){
       sprintf(filename_in, "%s/grids/snapshot_%03d_dark_grid.dat",filename_SSimPL_in, i_snap);
       check_integrity_grids(filename_in);
    }
    // Process snapshots and their smooth files
    //if(flag_process_snapshots){
    //   int IDs_byte_size;
    //   int i_region=-1;
    //   if(check_integrity_snapshot(filename_SSimPL_in,i_region,i_snap,&IDs_byte_size))
    //      check_integrity_smooth(filename_SSimPL_in,i_region,i_snap,IDs_byte_size);
    //   i_region++;
    //   while(check_integrity_snapshot(filename_SSimPL_in,i_region,i_snap,&IDs_byte_size)){
    //      check_integrity_smooth(filename_SSimPL_in,i_region,i_snap,IDs_byte_size);
    //      i_region++;
    //   }
    //}
    SID_log("Done.",SID_LOG_CLOSE);
  }

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

