/********************************************************/
/*                                                      */
/*                    random_number.c                   */
/*                    ---------------                   */
/*                                                      */
/********************************************************/
#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
int main(int argc, char *argv[]){
    int     column_1,column_2;
    char    filename[MAX_FILENAME_LENGTH];
    char   *line=NULL;
    size_t  line_length=0;
    FILE   *fp;

    // Check syntax and read-in input parameters 
    if(argc!=4){
        printf("\n  Syntax: %s filename column_1 column_2\n",
            argv[0]);
        return(ERROR_SYNTAX);
    }

    SID_init(&argc,&argv,NULL);
    strcpy(filename,   argv[1]);
    column_1=(int)atoi(argv[2]);
    column_2=(int)atoi(argv[3]);

    // Open and read file
    SID_log("Reading file...",SID_LOG_OPEN|SID_LOG_TIMER);
    size_t  n_data,i_data;
    float  *data_1;
    float  *data_2;
    fp    =fopen(filename,"r");
    n_data=(size_t)count_lines_data(fp);
    data_1=(float *)SID_malloc(sizeof(float)*n_data);
    data_2=(float *)SID_malloc(sizeof(float)*n_data);
    for(i_data=0;i_data<n_data;i_data++){
      grab_next_line_data(fp,&line,&line_length);
      grab_float(line,column_1,&(data_1[i_data]));
      grab_float(line,column_2,&(data_2[i_data]));
    }
    SID_log("Done.",SID_LOG_CLOSE);

    float *result;
    result=(float *)SID_malloc(sizeof(float)*n_data);
    SID_log("Performing multiplication...",SID_LOG_OPEN|SID_LOG_TIMER);
    for(i_data=0;i_data<n_data;i_data++) result[i_data]=data_1[i_data]*data_2[i_data];
    SID_log("Done.",SID_LOG_CLOSE);

    // Perform multiplication
    SID_log("Performing multiplication...",SID_LOG_OPEN|SID_LOG_TIMER);
    calc_array_multiply(data_1,data_2,NULL,n_data,SID_FLOAT,CALC_MODE_DEFAULT);
    SID_log("Done.",SID_LOG_CLOSE);

    // Write Result
    SID_log("Writing results...",SID_LOG_OPEN|SID_LOG_TIMER);
    for(i_data=0;i_data<n_data;i_data++) printf("%lf %lf\n",data_1[i_data],result[i_data]);
    SID_free(SID_FARG result);
    SID_log("Done.",SID_LOG_CLOSE);

    // Clean-up
    SID_free(SID_FARG data_1);
    SID_free(SID_FARG data_2);

    SID_exit(ERROR_NONE);
}

