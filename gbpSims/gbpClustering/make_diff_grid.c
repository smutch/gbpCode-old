#define  _MAIN
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpCosmo.h>
#include <gbpClustering.h>

#define _FILE_OFFSET_BITS 64

int main(int argc, char *argv[]){

  // Initialization -- MPI etc.
  SID_init(&argc,&argv,NULL);

  SID_log("Writing the difference between the grids in {%s} and {%s} to stdout...",SID_LOG_OPEN,argv[1],argv[2]);

  FILE   *fp1;
  int     nx1,ny1,nz1,ng1,as1;
  double  Lx1,Ly1,Lz1;
  fp1=fopen(argv[1],"r");
  fread(&nx1,sizeof(int),   1,fp1);
  fread(&ny1,sizeof(int),   1,fp1);
  fread(&nz1,sizeof(int),   1,fp1);
  fread(&Lx1,sizeof(double),1,fp1);
  fread(&Ly1,sizeof(double),1,fp1);
  fread(&Lz1,sizeof(double),1,fp1);
  fread(&ng1,sizeof(int),   1,fp1);
  fread(&as1,sizeof(int),   1,fp1);

  FILE   *fp2;
  int     nx2,ny2,nz2,ng2,as2;
  double  Lx2,Ly2,Lz2;
  fp2=fopen(argv[2],"r");
  fread(&nx2,sizeof(int),   1,fp2);
  fread(&ny2,sizeof(int),   1,fp2);
  fread(&nz2,sizeof(int),   1,fp2);
  fread(&Lx2,sizeof(double),1,fp2);
  fread(&Ly2,sizeof(double),1,fp2);
  fread(&Lz2,sizeof(double),1,fp2);
  fread(&ng2,sizeof(int),   1,fp2);
  fread(&as2,sizeof(int),   1,fp2);

  if(nx1!=nx2) SID_trap_error("nx's don't match (ie. %d!=%d)",ERROR_LOGIC,nx1,nx2);
  if(ny1!=ny2) SID_trap_error("ny's don't match (ie. %d!=%d)",ERROR_LOGIC,ny1,ny2);
  if(nz1!=nz2) SID_trap_error("nz's don't match (ie. %d!=%d)",ERROR_LOGIC,nz1,nz2);
  if(Lx1!=Lx2) SID_trap_error("Lx's don't match (ie. %le!=%le)",ERROR_LOGIC,Lx1,Lx2);
  if(Ly1!=Ly2) SID_trap_error("Ly's don't match (ie. %le!=%le)",ERROR_LOGIC,Ly1,Ly2);
  if(Lz1!=Lz2) SID_trap_error("Lz's don't match (ie. %le!=%le)",ERROR_LOGIC,Lz1,Lz2);
  if(ng1!=ng2) SID_trap_error("The number of grids don't match (ie. %d!=%d)",ERROR_LOGIC,ng1,ng2);

  int i_x,i_y,i_z;
  fftw_real d1,d2;
  int i_grid;
  char *grid_identifier_1;
  char *grid_identifier_2;
  grid_identifier_1=(char *)SID_malloc(GRID_IDENTIFIER_SIZE*sizeof(char));
  grid_identifier_2=(char *)SID_malloc(GRID_IDENTIFIER_SIZE*sizeof(char));
  for(i_grid=0;i_grid<ng1;i_grid++){
     fread(grid_identifier_1,sizeof(char),GRID_IDENTIFIER_SIZE,fp1);
     fread(grid_identifier_2,sizeof(char),GRID_IDENTIFIER_SIZE,fp2);
     if(strcmp(grid_identifier_1,grid_identifier_2))
        SID_log_warning("grid identifiers don't match (ie. {%s}!={%s})",ERROR_LOGIC,grid_identifier_1,grid_identifier_2);
     SID_log("Processing {%s} ...",SID_LOG_OPEN,grid_identifier_1);
     for(i_x=0;i_x<nx1;i_x++){
        for(i_y=0;i_y<ny1;i_y++){
           for(i_z=0;i_z<nz1;i_z++){
                fread(&d1,sizeof(fftw_real),1,fp1);
                fread(&d2,sizeof(fftw_real),1,fp2);
                if(d2!=d1)
                   fprintf(stdout,"%3d %4d %4d %4d %le %le %le\n",i_grid,i_x,i_y,i_z,(double)d1,(double)d2,(double)(d2-d1));
           }
        }
     }
     SID_log("Done.",SID_LOG_CLOSE);
  }
  SID_free(SID_FARG grid_identifier_1);
  SID_free(SID_FARG grid_identifier_2);
  fclose(fp1); 
  fclose(fp2); 

  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}

