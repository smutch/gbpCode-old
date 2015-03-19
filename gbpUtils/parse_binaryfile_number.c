#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
int main(int argc, char *argv[]){

    SID_init(&argc,&argv,NULL,NULL);

    // Parse arguments
    if(argc<4)
       SID_trap_error("Invlaid syntax.",ERROR_SYNTAX);
    char filename[MAX_FILENAME_LENGTH];
    char dtype[32];
    int  offset;
    int  i_arg=1;
    strcpy(filename,argv[i_arg++]);
    offset=atoi(    argv[i_arg++]);

    // Open file
    FILE *fp;
    if((fp=fopen(filename,"r"))==NULL)
       SID_trap_error("Could not open file {%s}.",ERROR_IO_OPEN,filename);

    while(i_arg<argc){
       // Move to offset
       fseeko(fp,(off_t)(offset),SEEK_SET);

       // Parse entry
       strcpy(dtype,argv[i_arg++]);
       if(!strcmp(dtype,"int")){
          int buffer;
          fread(&buffer,sizeof(int),1,fp);
          printf("%d\n",buffer);
          offset+=sizeof(int);
       }
       else if(!strcmp(dtype,"long")){
          long long buffer;
          fread(&buffer,sizeof(long long),1,fp);
          printf("%lld\n",buffer);
          offset+=sizeof(long long);
       }
       else if(!strcmp(dtype,"float")){
          float buffer;
          fread(&buffer,sizeof(float),1,fp);
          printf("%f\n",buffer);
          offset+=sizeof(float);
       }
       else if(!strcmp(dtype,"double")){
          double buffer;
          fread(&buffer,sizeof(double),1,fp);
          printf("%llf\n",buffer);
          offset+=sizeof(double);
       }
       else
          SID_trap_error("Unsupported dtype {%s}.",ERROR_LOGIC,dtype);
    }

    // Close file
    fclose(fp);

    SID_exit(ERROR_NONE);
}

