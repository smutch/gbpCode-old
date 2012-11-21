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

  FILE   *fp1;
  int     nx1,ny1,nz1,nb1;
  double  Lx1,Ly1,Lz1;
  fp1=fopen(argv[1],"r");
  fread(&nx1,sizeof(int),   1,fp1);
  fread(&ny1,sizeof(int),   1,fp1);
  fread(&nz1,sizeof(int),   1,fp1);
  fread(&Lx1,sizeof(double),1,fp1);
  fread(&Ly1,sizeof(double),1,fp1);
  fread(&Lz1,sizeof(double),1,fp1);
  fread(&nb1,sizeof(int),   1,fp1);

  FILE   *fp2;
  int     nx2,ny2,nz2,nb2;
  double  Lx2,Ly2,Lz2;
  fp2=fopen(argv[2],"r");
  fread(&nx2,sizeof(int),   1,fp2);
  fread(&ny2,sizeof(int),   1,fp2);
  fread(&nz2,sizeof(int),   1,fp2);
  fread(&Lx2,sizeof(double),1,fp2);
  fread(&Ly2,sizeof(double),1,fp2);
  fread(&Lz2,sizeof(double),1,fp2);
  fread(&nb2,sizeof(int),   1,fp2);

  if(nx1!=nx2) SID_trap_error("nx's don't match (ie. %d!=%d)",ERROR_LOGIC,nx1,nx2);
  if(ny1!=ny2) SID_trap_error("ny's don't match (ie. %d!=%d)",ERROR_LOGIC,ny1,ny2);
  if(nz1!=nz2) SID_trap_error("nz's don't match (ie. %d!=%d)",ERROR_LOGIC,nz1,nz2);
  if(Lx1!=Lx2) SID_trap_error("Lx's don't match (ie. %le!=%le)",ERROR_LOGIC,Lx1,Lx2);
  if(Ly1!=Ly2) SID_trap_error("Ly's don't match (ie. %le!=%le)",ERROR_LOGIC,Ly1,Ly2);
  if(Lz1!=Lz2) SID_trap_error("Lz's don't match (ie. %le!=%le)",ERROR_LOGIC,Lz1,Lz2);
  if(nb1!=nb2) SID_trap_error("The number of boxes don't match (ie. %d!=%d)",ERROR_LOGIC,nb1,nb2);

  int i_x,i_y,i_z;
  fftw_real d1,d2;
  for(i_x=0;i_x<nx1;i_x++){
     for(i_y=0;i_y<ny1;i_y++){
        for(i_z=0;i_z<nz1;i_z++){
             fread(&d1,sizeof(fftw_real),1,fp1);
             fread(&d2,sizeof(fftw_real),1,fp2);
             fprintf(stdout,"%4d %4d %4d %le %le %le\n",i_x,i_y,i_z,(double)d1,(double)d2,(double)(d2-d1));
        }
     }
  }

  fclose(fp1); 
  fclose(fp2); 

  SID_exit(ERROR_NONE);
}


