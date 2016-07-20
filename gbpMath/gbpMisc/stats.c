#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <common.h>

int main(int argc, char *argv[]){
  char    filename_in[256];
  int     data_column;
  char   *line_in;
  int     line_length=0;
  double *data;
  int     i_data;
  double  min,max,mean,median,std_dev,sum;
  size_t  n_data_local,n_data;
  FILE   *unit_in;

  // Initialize gbpLib
  SID_init(&argc,&argv,NULL,NULL);

  // Check syntax and read-in input parameters
  if(argc!=3)
    SID_trap_error("Syntax: %s file_in column",ERROR_SYNTAX,argv[0]);
  strcpy(filename_in, argv[1]);
  data_column=atoi(argv[2]);

  // Open input file
  if((unit_in=fopen(filename_in,"r"))==NULL)
    SID_trap_error("Error opening input file {%s}.\n",ERROR_IO_READ,filename_in);

  // Count number of lines in the input file
  n_data_local=count_lines_data(unit_in);
 
  // Allocate memory for the data. Read it and sort it in ascending order 
  data=(double *)malloc(sizeof(double)*n_data_local);
  for(i_data=0;i_data<n_data_local;i_data++){
    grab_next_line_data(unit_in,&line_in,&line_length);
    grab_double(line_in,data_column,&(data[i_data]));
  }
  fclose(unit_in);

  // Build statistics of data
  n_data=calc_stat(data,NULL,n_data_local,ADaM_DOUBLE,CALC_STAT_GLOBAL|CALC_STAT_MIN,   &min);
  n_data=calc_stat(data,NULL,n_data_local,ADaM_DOUBLE,CALC_STAT_GLOBAL|CALC_STAT_MAX,   &max);
  n_data=calc_stat(data,NULL,n_data_local,ADaM_DOUBLE,CALC_STAT_GLOBAL|CALC_STAT_SUM,   &sum);
  n_data=calc_stat(data,NULL,n_data_local,ADaM_DOUBLE,CALC_STAT_GLOBAL|CALC_STAT_MEAN,  &mean);
  n_data=calc_stat(data,NULL,n_data_local,ADaM_DOUBLE,CALC_STAT_GLOBAL|CALC_STAT_STDDEV,&std_dev);
  n_data=calc_stat(data,NULL,n_data_local,ADaM_DOUBLE,CALC_STAT_GLOBAL|CALC_STAT_MEDIAN,&median);

  // Print stats
  printf("Stats for column #%d of %s:\n",data_column,filename_in);
  printf("  min      = %11.4e\n",min);
  printf("  max      = %11.4e\n",max);
  printf("  sum      = %11.4e\n",sum);
  printf("  mean     = %11.4e\n",mean);
  printf("  median   = %11.4e\n",median);
  printf("  std_dev  = %11.4e\n",std_dev);
  printf("  n_data   =  %lld\n", n_data);

  // Clean-up
  free(data);
  SID_exit(ERROR_NONE);
}
