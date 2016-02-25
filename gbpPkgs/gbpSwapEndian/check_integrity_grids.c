#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpClustering.h>
#include <gbpSwapEndian.h>

void check_integrity_grids(const char *filename_in){

  SID_log("Swapping endian of grids...",SID_LOG_OPEN);

  // Open input files
  FILE *fp_in =NULL;
  if((fp_in=fopen(filename_in,"r"))==NULL)
     SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);

  // Read the needed header information and rewind
  int    n[3];
  double L[3];
  int    n_grids;
  int    scheme;
  size_t n_bytes=0;
  fread_verify(&(n[0]), sizeof(int),   1,fp_in);n_bytes+=sizeof(int);
  fread_verify(&(n[1]), sizeof(int),   1,fp_in);n_bytes+=sizeof(int);
  fread_verify(&(n[2]), sizeof(int),   1,fp_in);n_bytes+=sizeof(int);
  fread_verify(&(L[0]), sizeof(double),1,fp_in);n_bytes+=sizeof(double);
  fread_verify(&(L[1]), sizeof(double),1,fp_in);n_bytes+=sizeof(double);
  fread_verify(&(L[2]), sizeof(double),1,fp_in);n_bytes+=sizeof(double);
  fread_verify(&n_grids,sizeof(int),   1,fp_in);n_bytes+=sizeof(int);
  fread_verify(&scheme, sizeof(int),   1,fp_in);n_bytes+=sizeof(int);
  int grid_size=n[0]*n[1]*n[2];

  // Create a read buffer
  char *buffer=(char *)SID_malloc(sizeof(char)*grid_size*sizeof(fftw_real));

  // Process the file
  for(int i_grid=0;i_grid<n_grids;i_grid++){
     fread_verify(buffer,GRID_IDENTIFIER_SIZE,sizeof(char),     fp_in);n_bytes+=GRID_IDENTIFIER_SIZE*sizeof(char);
     fread_verify(buffer,grid_size,           sizeof(fftw_real),fp_in);n_bytes+=grid_size*sizeof(fftw_real);
  }

  // Check that we are at the end of the file
  char test;fread(&test,1,1,fp_in);
  if(!feof(fp_in))
     SID_trap_error("There are stray bytes at the end of {%s} (%lld bytes read).",ERROR_LOGIC,filename_in,n_bytes);

  // Free the read buffer
  SID_free(SID_FARG buffer);

  // Close files
  fclose(fp_in);
  
  SID_log("Done.",SID_LOG_CLOSE);
}

