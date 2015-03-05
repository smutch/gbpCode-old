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
  if(argc!=7) SID_trap_error("Invalid syntax.",ERROR_SYNTAX);
  char filename_SSimPL_in[MAX_FILENAME_LENGTH];
  char filename_SSimPL_out[MAX_FILENAME_LENGTH];
  char filename_halo_type[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_in, argv[1]);
  strcpy(filename_SSimPL_out,argv[2]);
  strcpy(filename_halo_type, argv[3]);
  int i_snap_lo=atoi(argv[4]);
  int i_snap_hi=atoi(argv[5]);
  int mode_in  =atoi(argv[6]);

  // Parse the mode flag
  int mode=SWAP_SSIMPL_ENDIAN_TO_NATIVE;
  if(mode_in!=0)
     mode=SWAP_SSIMPL_ENDIAN_FROM_NATIVE;

  // Create output directories
  char dir_snapshots[MAX_FILENAME_LENGTH];
  char dir_smooth[MAX_FILENAME_LENGTH];
  char dir_grids[MAX_FILENAME_LENGTH];
  char dir_halos[MAX_FILENAME_LENGTH];
  char dir_catalogs[MAX_FILENAME_LENGTH];
  sprintf(dir_snapshots,"%s/snapshots",filename_SSimPL_out);
  sprintf(dir_smooth,   "%s/smooth",   filename_SSimPL_out);
  sprintf(dir_grids,    "%s/grids",    filename_SSimPL_out);
  sprintf(dir_halos,    "%s/halos",    filename_SSimPL_out);
  sprintf(dir_catalogs, "%s/catalogs", filename_SSimPL_out);
  mkdir(filename_SSimPL_out,02755);
  mkdir(dir_snapshots,      02755);
  mkdir(dir_smooth,         02755);
  mkdir(dir_grids,          02755);
  mkdir(dir_halos,          02755);
  mkdir(dir_catalogs,       02755);

  // Set which files will be processed
  int flag_process_halos    =TRUE;
  int flag_process_catalogs =TRUE;
  int flag_process_grids    =TRUE;
  int flag_process_snapshots=TRUE;

  // Loop over the given snapshot range
  SID_log("Processing group/subgroup statistics for files #%d->#%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_snap_lo,i_snap_hi);
  for(int i_snap=i_snap_lo;i_snap<=i_snap_hi;i_snap++){
    char filename_in[MAX_FILENAME_LENGTH];
    char filename_out[MAX_FILENAME_LENGTH];
    SID_log("Processing snapshot #%03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_snap);
    // Process halos
    if(flag_process_halos){
       sprintf(filename_in, "%s/halos/",filename_SSimPL_in);
       sprintf(filename_out,"%s/halos/",filename_SSimPL_out);
       swap_endian_halos(filename_in,filename_out,filename_halo_type,i_snap,mode);
    }
    // Process catalogs
    if(flag_process_catalogs){
       sprintf(filename_in, "%s/catalogs/",filename_SSimPL_in, filename_halo_type);
       sprintf(filename_out,"%s/catalogs/",filename_SSimPL_out,filename_halo_type);
       swap_endian_catalogs(filename_in,filename_out,filename_halo_type,i_snap,mode);
    }
    // Process grids
    if(flag_process_grids){
       sprintf(filename_in, "%s/grids/snapshot_%03d_dark_grid.dat",filename_SSimPL_in, i_snap);
       sprintf(filename_out,"%s/grids/snapshot_%03d_dark_grid.dat",filename_SSimPL_out,i_snap);
       swap_endian_grids(filename_in,filename_out,mode);
    }
    // Process snapshots and their smooth files
    if(flag_process_snapshots){
       int IDs_byte_size;
       int i_region=-1;
       if(swap_endian_snapshot(filename_SSimPL_in,filename_SSimPL_out,i_region,i_snap,mode,&IDs_byte_size))
          swap_endian_smooth(filename_SSimPL_in,filename_SSimPL_out,i_region,i_snap,mode,IDs_byte_size);
       i_region++;
       while(swap_endian_snapshot(filename_SSimPL_in,filename_SSimPL_out,i_region,i_snap,mode,&IDs_byte_size)){
          swap_endian_smooth(filename_SSimPL_in,filename_SSimPL_out,i_region,i_snap,mode,IDs_byte_size);
          i_region++;
       }
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

