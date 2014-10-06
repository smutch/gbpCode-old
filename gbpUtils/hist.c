/********************************************************/
/*                                                      */
/*                       hist.c                         */
/*                       ------                         */
/*                                                      */
/*  This program takes as input a data file with an     */
/*  unspecified number of columns and generates a       */
/*  histogram (and various supporting statistics) for   */
/*  a chosen column.  The histogram is written to an    */
/*  output file while the statistics are written to     */
/*  stdout.                                             */
/*                                                      */
/*  SYNTAX: hist filein fileout column n_bins [min max] */
/*  ------                                              */
/*                                                      */
/*  Created by:   Greg Poole                            */
/*  Last updated: Oct 25/2000                           */
/*                                                      */
/********************************************************/
#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>

#define CLEAN_EXIT        0
#define SYNTAX_ERROR      1
#define OPEN_INPUT_ERROR  2
#define OPEN_OUTPUT_ERROR 3

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
    int     i_data,j_data,i_bin;
    int     cumulator_up;
    int     cumulator_dn;
    double  temp_value;
    size_t  n_data;
    size_t  n_data_used;
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

    /*********************************************/
    /* Check syntax and read-in input parameters */
    /*********************************************/
    if(argc<4 || argc>6 || argc==5){
        printf("\n  Syntax: %s file_in column n_bins [min max]\n\n",
            argv[0]);
        return(SYNTAX_ERROR);
    }
    else if(argc==6){
        min_bin=(double)atof(argv[4]);
        max_bin=(double)atof(argv[5]);
    }
    strcpy(filename_in,argv[1]);
    data_column  =atoi(argv[2]);
    n_bins       =atoi(argv[3]);

    // If the column is negative, take the log of it
    int flag_take_log=FALSE;
    if(data_column<0){
       data_column*=-1;
       flag_take_log=TRUE;
    }

    /*******************/
    /* Open input file */
    /*******************/
    if((unit_in=fopen(filename_in,"r"))==NULL){
        fprintf(stderr,"Error opening input file {%s}.\n",filename_in);
        return(OPEN_INPUT_ERROR);
    }

    SID_log("Compiling stats and histogram of column #%d of {%s}...",SID_LOG_OPEN,data_column,filename_in);
    if(flag_take_log) SID_log("Taking the log of the data.",SID_LOG_COMMENT);
    SID_log("Stats:",SID_LOG_COMMENT,data_column);

    /*******************************************/
    /* Count number of lines in the input file */
    /*******************************************/
    n_data=count_lines_data(unit_in);
    SID_log("  n_data   = %11d",    SID_LOG_COMMENT,n_data);
 
    /************************************************************************/
    /* Allocate memory for the data. Read it and sort it in ascending order */
    /************************************************************************/
    data=(double *)SID_malloc(sizeof(double)*n_data);
    for(i=0;i<n_data;i++){
      grab_next_line_data(unit_in,&line_in,&line_length);
      grab_double(line_in,data_column,&(data[i]));
    }
    fclose(unit_in);

    // Take the log of the data if we have been asked to
    if(flag_take_log) for(i=0;i<n_data;i++) data[i]=take_log10(data[i]);

    /******************************************/
    /* Build simple statistics of data column */
    /******************************************/
    calc_min(data,&min,n_data,   SID_DOUBLE,CALC_MODE_DEFAULT);
    SID_log("  min      = %11.4e",SID_LOG_COMMENT,min,CALC_MODE_DEFAULT);
    calc_max(data,&max,n_data,   SID_DOUBLE,CALC_MODE_DEFAULT);
    SID_log("  max      = %11.4e",SID_LOG_COMMENT,max,CALC_MODE_DEFAULT);
    calc_sum(data,&sum,n_data,   SID_DOUBLE,CALC_MODE_DEFAULT);
    SID_log("  sum      = %11.4e",SID_LOG_COMMENT,sum,CALC_MODE_DEFAULT);
    calc_mean(data,&mean,n_data,  SID_DOUBLE,CALC_MODE_DEFAULT);
    SID_log("  mean     = %11.4e",SID_LOG_COMMENT,mean,CALC_MODE_DEFAULT);
    calc_median(data,&median,n_data,SID_DOUBLE,CALC_MODE_DEFAULT);
    SID_log("  median   = %11.4e",SID_LOG_COMMENT,median,CALC_MODE_DEFAULT);
    calc_stddev(data,&std_dev,n_data,SID_DOUBLE,CALC_MODE_DEFAULT);
    SID_log("  std_dev  = %11.4e",SID_LOG_COMMENT,std_dev);

    if(argc!=6){
        min_bin=min;
        max_bin=max;
    }
    bin_size  =(max_bin-min_bin)/(double)(n_bins);
    bin       =(double *)SID_malloc(sizeof(double)*(n_bins+1));
    bin_median=(double *)SID_malloc(sizeof(double)* n_bins);
    hist      =(int    *)SID_malloc(sizeof(int)*    n_bins);
    bin[0] =min_bin;
    hist[0]=0;
    for(i=1;i<n_bins;i++){
        bin[i] =((double)i)*bin_size+min_bin;
        hist[i]=0;
    }
    bin[n_bins]=max_bin;
    for(i=0;i<n_bins;i++)
       bin_median[i]=0.5*(bin[i]+bin[i+1]); // default value for empty bins

    SID_log("Properties of histogram:",SID_LOG_COMMENT);
    SID_log("  n_bins   = %11d",   SID_LOG_COMMENT,n_bins);
    SID_log("  min_bin  = %11.4le",SID_LOG_COMMENT,min_bin);
    SID_log("  max_bin  = %11.4le",SID_LOG_COMMENT,max_bin);
    SID_log("  bin_size = %11.4le",SID_LOG_COMMENT,bin_size);

    /*********************/
    /* Compile histogram */
    /*********************/
    SID_log("Compiling histogram...",SID_LOG_OPEN|SID_LOG_TIMER);
    SID_log("Sorting...",SID_LOG_OPEN|SID_LOG_TIMER);
    size_t *data_index=NULL;
    merge_sort(data,(size_t)n_data,&data_index,SID_DOUBLE,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
    SID_log("Done.",SID_LOG_CLOSE);
    SID_log("Counting...",SID_LOG_OPEN|SID_LOG_TIMER);
    int data_bin_last=-1;
    int i_data_start =-1;
    int i_data_stop  =-1;
    int flag_finalize=FALSE;
    for(i_data=0,n_data_used=0;i_data<n_data;i_data++){
      size_t idx_i     =data_index[i_data];
      double data_bin_d=(data[idx_i]-min_bin)/bin_size;
      int    data_bin;
      if(data_bin_d<0)
         data_bin=-1;
      else
         data_bin=(int)floor(data_bin_d); 
      if(data_bin>=0 && data_bin<n_bins){
         if(data_bin>data_bin_last){
            if(data_bin_last>=0){
               switch(hist[data_bin_last]%2){
                 case 0:
                   bin_median[data_bin_last]=0.5*(data[data_index[i_data_start+hist[data_bin_last]/2-1]]+
                                                  data[data_index[i_data_start+hist[data_bin_last]/2+1]]);
                   break;
                 case 1:
                   bin_median[data_bin_last]=data[data_index[i_data_start+hist[data_bin_last]/2]];
                   break;
               }
            }
            i_data_start =i_data;
            flag_finalize=TRUE;
         }
         data_bin_last=data_bin;
         i_data_stop  =i_data;
         hist[data_bin_last]++;
         n_data_used++;
      }
    }
    if(flag_finalize){
       if(data_bin_last>=0){
          switch(hist[data_bin_last]%2){
            case 0:
              bin_median[data_bin_last]=0.5*(data[data_index[i_data_start+hist[data_bin_last]/2-1]]+
                                             data[data_index[i_data_start+hist[data_bin_last]/2+1]]);
              break;
            case 1:
              bin_median[data_bin_last]=data[data_index[i_data_start+hist[data_bin_last]/2]];
              break;
          }
       }
       flag_finalize=FALSE;
    }
    SID_log("Done.",SID_LOG_CLOSE);
    SID_free(SID_FARG data_index);
    SID_log("Done.",SID_LOG_CLOSE);
       
    /********************/
    /* Open output file */
    /********************/
    sprintf(filename_out,"%s.hist.%05d",filename_in,data_column);
    if((unit_out=fopen(filename_out,"w"))==NULL){
        fprintf(stderr,"Error opening output file {%s}.\n",filename_out);
        SID_free(SID_FARG data);
        return(OPEN_OUTPUT_ERROR);
    }

    /****************/
    /* Write header */
    /****************/
    fprintf(unit_out,"# Analysis of column #%d of {%s}:\n",data_column,filename_in);
    fprintf(unit_out,"#\n");
    fprintf(unit_out,"# Stats:\n");
    fprintf(unit_out,"#   n_data   = %zd\n",   n_data);
    fprintf(unit_out,"#   min      = %11.4e\n",min);
    fprintf(unit_out,"#   max      = %11.4e\n",max);
    fprintf(unit_out,"#   sum      = %11.4e\n",sum);
    fprintf(unit_out,"#   mean     = %11.4e\n",mean);
    fprintf(unit_out,"#   median   = %11.4e\n",median);
    fprintf(unit_out,"#   std_dev  = %11.4e\n",std_dev);
    fprintf(unit_out,"#\n");
    fprintf(unit_out,"# Properties of histogram:\n");
    fprintf(unit_out,"#   n_bins   = %11d\n",  n_bins);
    fprintf(unit_out,"#   min_bin  = %11.4e\n",min_bin);
    fprintf(unit_out,"#   max_bin  = %11.4e\n",max_bin);
    fprintf(unit_out,"#   bin_size = %11.4e\n",bin_size);
    fprintf(unit_out,"#\n");
    fprintf(unit_out,"# Column: (1) Histogram bin lo\n");
    fprintf(unit_out,"#         (2) Histogram bin median\n");
    fprintf(unit_out,"#         (3) Histogram bin hi\n");
    fprintf(unit_out,"#         (4) Bin counts\n");
    fprintf(unit_out,"#         (5) Cumulative bin counts (>lo bin)\n");
    fprintf(unit_out,"#         (6) Cumulative bin counts (<hi bin)\n");
    fprintf(unit_out,"#\n");

    /*******************/
    /* Write histogram */
    /*******************/
    SID_log("Writing to {%s}...",SID_LOG_OPEN,filename_out);
    bin_start=0;
    for(i=0,cumulator_up=0,cumulator_dn=n_data_used;i<n_bins;i++){
      cumulator_up+=hist[i];
      fprintf(unit_out,"%11.4f %11.4f %11.4f %d %d %d\n",
              bin[i],bin_median[i],bin[i+1],hist[i],cumulator_dn,cumulator_up);
      cumulator_dn-=hist[i];
    }
    fclose(unit_out);
    SID_log("Done.",SID_LOG_CLOSE);

    /*************************/
    /* Free allocated memory */
    /*************************/
    SID_free(SID_FARG data);
    SID_free(SID_FARG bin);
    SID_free(SID_FARG bin_median);
    SID_free(SID_FARG hist);
    SID_log("Done.",SID_LOG_CLOSE);

    SID_exit(ERROR_NONE);
}
