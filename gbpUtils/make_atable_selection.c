#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>

int main(int argc, char *argv[]){

    SID_init(&argc,&argv,NULL);

    if(argc!=5){
        printf("\n  Syntax: %s file_in column check_option value\n\n",
            argv[0]);
        return(ERROR_SYNTAX);
    }
    char filename_in[MAX_FILENAME_LENGTH];
    char filename_out[MAX_FILENAME_LENGTH];
    char check_option_in[8];
    strcpy(filename_in,argv[1]);
    int data_column=atoi(argv[2]);
    strcpy(check_option_in,argv[3]);
    double value=(double)atof(argv[4]);
    // Make sure the check_option is valid
    int check_option=-1;
    if(!strcmp(check_option_in,"gt"))      check_option=0;
    else if(!strcmp(check_option_in,"ge")) check_option=1;
    else if(!strcmp(check_option_in,"lt")) check_option=2;
    else if(!strcmp(check_option_in,"le")) check_option=3;
    else SID_trap_error("Invalid check_option {%s}.",ERROR_SYNTAX,check_option_in);

    SID_log("Selecting lines from {%s}...",SID_LOG_OPEN,filename_in);

    // Open input file
    FILE *fp_in;
    if((fp_in=fopen(filename_in,"r"))==NULL){
        fprintf(stderr,"Error opening input file {%s}.\n",filename_in);
        return(ERROR_IO_OPEN);
    }

    // Open output file
    FILE *fp_out;
    sprintf(filename_out,"%s.col_%d_%s_%le",filename_in,data_column,check_option_in,value);
    if((fp_out=fopen(filename_out,"w"))==NULL){
        fprintf(stderr,"Error opening output file {%s}.\n",filename_out);
        return(ERROR_IO_OPEN);
    }

    // Count number of lines in the input file 
    int n_lines=count_lines(fp_in);

    // Perform trim
    int flag_keep;
    int n_trim=0;
    char   *line_in=NULL;
    size_t  line_length=0;
    SID_log("%d lines being processed...",SID_LOG_OPEN,n_lines);
    for(int i=0;i<n_lines;i++){
      grab_next_line(fp_in,&line_in,&line_length);
      flag_keep=TRUE;
      
      if(!check_comment(line_in)){
         double data;
         grab_double(line_in,data_column,&data);
         switch(check_option){
            case 0:
              flag_keep=(data>value);
              break;
            case 1:
              flag_keep=(data>=value);
              break;
            case 2:
              flag_keep=(data<value);
              break;
            case 3:
              flag_keep=(data<=value);
              break;
         }
      }
      if(flag_keep) fprintf(fp_out,"%s",line_in);
      else          n_trim++;
    }
    SID_log("%d removed...",SID_LOG_CONTINUE,n_trim);
    fclose(fp_in);
    fclose(fp_out);
    SID_free(SID_FARG line_in);
    SID_log("Done.",SID_LOG_CLOSE);

    SID_log("Output written to {%s}.",SID_LOG_COMMENT,filename_out);

    SID_log("Done.",SID_LOG_CLOSE);
    SID_exit(ERROR_NONE);
}

