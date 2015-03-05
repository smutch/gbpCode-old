#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpClustering.h>
#include <gbpSims_Swap_Endian.h>

void swap_endian_grids(const char *filename_in,const char *filename_out,int mode){

  SID_log("Swapping endian of grids...",SID_LOG_OPEN);

  // Sanity check
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE) && check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE))
     SID_trap_error("Invalid mode flag (%d) in swap_endian_grids().",ERROR_LOGIC,mode);

  // Open input and output files
  FILE *fp_in =NULL;
  FILE *fp_out=NULL;
  if((fp_in=fopen(filename_in,"r"))==NULL)
     SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);
  if((fp_out=fopen(filename_out,"w"))==NULL)
     SID_trap_error("Could not open {%s} for writing.",ERROR_IO_OPEN,filename_out);

  // Read the needed header information and rewind
  int    n[3];
  double L[3];
  int    n_grids;
  fread(&(n[0]), sizeof(int),   1,fp_in);
  fread(&(n[1]), sizeof(int),   1,fp_in);
  fread(&(n[2]), sizeof(int),   1,fp_in);
  fread(&(L[0]), sizeof(double),1,fp_in);
  fread(&(L[1]), sizeof(double),1,fp_in);
  fread(&(L[2]), sizeof(double),1,fp_in);
  fread(&n_grids,sizeof(int),   1,fp_in);
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE)){
     swap_endian((char *)(n),       3,sizeof(int));
     swap_endian((char *)(&n_grids),1,sizeof(int));
  }
  int grid_size=n[0]*n[1]*n[2];
  rewind(fp_in);

  // Create a read buffer
  char *buffer=(char *)SID_malloc(sizeof(char)*grid_size*sizeof(fftw_real));

  // Process the file
  rewrite_swap_endian(fp_in,fp_out,3,sizeof(int),   buffer);
  rewrite_swap_endian(fp_in,fp_out,3,sizeof(double),buffer);
  rewrite_swap_endian(fp_in,fp_out,2,sizeof(int),   buffer);
  for(int i_grid=0;i_grid<n_grids;i_grid++){
     rewrite_swap_endian(fp_in,fp_out,GRID_IDENTIFIER_SIZE,sizeof(char),     buffer);
     rewrite_swap_endian(fp_in,fp_out,grid_size,           sizeof(fftw_real),buffer);
  }

  // Free the read buffer
  SID_free(SID_FARG buffer);

  // Close files
  fclose(fp_in);
  fclose(fp_out);
  
  SID_log("Done.",SID_LOG_CLOSE);
}

