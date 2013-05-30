#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>

int main(int argc, char *argv[]){
    char    filename_in[256];
    char    filename_out[256];
    char    temp_char[2],temp_char_old[2];
    int     data_column;
    int     n_bins;
    char   *line_in;
    size_t  line_length=0;
    double *data;
    int    *data_bin;
    size_t *data_bin_index;
    int     i_data,j_data,i_bin;
    int     cumulator;
    double  temp_value;
    size_t  n_data;
    int     i,j,k,l;
    int     j_bin;
    int     flag;
    int     i_median;
    double  min,max,mean,median,std_dev,sum;
    double  min_bin,max_bin,bin_size;
    double *bin;
    int    *hist;
    FILE   *unit_in;
    FILE   *unit_out;
    int     ir;
    double  rra;
    double *bin_median;
    int     bin_start;

    SID_init(&argc,&argv,NULL);

    // Process inputs
    strcpy(filename_in,argv[1]);
    data_column  =atoi(argv[2]);
    SID_log("Sorting column #%d of {%s}...",SID_LOG_OPEN,data_column,filename_in);

    // Open data file
    if((unit_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Error opening input file {%s}.",ERROR_IO_OPEN,filename_in);

    // Count number of lines in the input file
    n_data=count_lines_data(unit_in);
    SID_log("  n_data   = %11d",    SID_LOG_COMMENT,n_data);

    // Decide how many items go to each rank
    int n_left;
    int n_local;
    int n_rank;
    int i_rank;
    int n_max;
    for(i_rank=0,n_left=n_data,n_max=0;i_rank<SID.n_proc;i_rank++){
       if(i_rank==SID.n_proc-1)
         n_rank=n_left;
       else
         n_rank=(int)((float)n_left/(float)(SID.n_proc-i_rank));
       if(i_rank==SID.My_rank){
          n_local=n_rank;
       }
       n_left-=n_rank;
       n_max  =MAX(n_max,n_rank);
    }
 
    // Allocate memory for the data
    data=(double *)SID_malloc(sizeof(double)*n_local);
double *data_all      =NULL;
size_t *data_all_index=NULL;
size_t *data_all_rank =NULL;
data_all=(double *)SID_malloc(sizeof(double)*n_data);

    // Read the file
    double *buffer;
    buffer=(double *)SID_malloc(sizeof(double)*n_max);
    for(i_rank=0,n_left=n_data,j=0;i_rank<SID.n_proc;i_rank++){
       if(i_rank==SID.n_proc-1)
         n_rank=n_left;
       else
         n_rank=(int)((float)n_left/(float)(SID.n_proc-i_rank));
       // Perform read
       if(SID.I_am_Master){
         for(i=0;i<n_rank;i++){
           grab_next_line_data(unit_in,&line_in,&line_length);
           grab_double(line_in,data_column,&(buffer[i]));
         }
       }
       SID_Bcast(buffer,sizeof(double)*n_rank,MASTER_RANK,SID.COMM_WORLD);
       if(i_rank==SID.My_rank){
         n_local=n_rank;
         memcpy(data,buffer,sizeof(double)*n_local);
       }
if(SID.I_am_Master){
memcpy(&(data_all[j]),buffer,sizeof(double)*n_rank);
}
       j     +=n_rank;
       n_left-=n_rank;
    }
    fclose(unit_in);

    // Perform sort
    size_t *data_index=NULL;
    size_t *data_rank =NULL;
    sort(data,(size_t)n_local,&data_index,SID_DOUBLE,SORT_GLOBAL,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
    sort(data,(size_t)n_local,&data_rank, SID_DOUBLE,SORT_GLOBAL,SORT_COMPUTE_RANK, SORT_COMPUTE_NOT_INPLACE);

if(SID.I_am_Master){
merge_sort(data_all,(size_t)n_data,&data_all_index,SID_DOUBLE,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
merge_sort(data_all,(size_t)n_data,&data_all_rank, SID_DOUBLE,SORT_COMPUTE_RANK, SORT_COMPUTE_NOT_INPLACE);
}

    // Print results
    size_t *buffer_2;
    size_t *buffer_3;
    buffer_2=(size_t *)SID_malloc(sizeof(size_t)*n_max);
    buffer_3=(size_t *)SID_malloc(sizeof(size_t)*n_max);
    for(i_rank=0,j=0;i_rank<SID.n_proc;i_rank++){
       if(SID.My_rank==i_rank){
          n_rank=n_local;
          memcpy(buffer,  data,      n_rank*sizeof(double));
          memcpy(buffer_2,data_index,n_rank*sizeof(size_t));
          memcpy(buffer_3,data_rank, n_rank*sizeof(size_t));
       }
       SID_Bcast(&n_rank,  sizeof(int),          i_rank,SID.COMM_WORLD);
       SID_Bcast(buffer,   sizeof(double)*n_rank,i_rank,SID.COMM_WORLD);
       SID_Bcast(buffer_2, sizeof(size_t)*n_rank,i_rank,SID.COMM_WORLD);
       SID_Bcast(buffer_3, sizeof(size_t)*n_rank,i_rank,SID.COMM_WORLD);
       if(SID.I_am_Master){
         for(i=0;i<n_rank;i++,j++){
           printf("%3d %6d %le %le %6zd %6zd %6zd %6zd\n",i_rank,j,buffer[i],data_all[j],buffer_2[i],data_all_index[j],buffer_3[i],data_all_rank[j]);
         }
       }
    }
    SID_free(SID_FARG buffer);
    SID_free(SID_FARG buffer_2);
    SID_free(SID_FARG buffer_3);
    SID_free(SID_FARG data);
    SID_free(SID_FARG data_index);
    SID_free(SID_FARG data_rank);
    SID_free(SID_FARG data_all_index);
    SID_free(SID_FARG data_all_rank);

    SID_exit(ERROR_NONE);
}
