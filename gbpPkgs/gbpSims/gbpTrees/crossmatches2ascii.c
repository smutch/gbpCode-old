#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){
  int     n_search;
  int     i_halo;
  char    group_text_prefix[4];
  int     n_files;
  int     k_read;
  int     max_n_groups;
  int     l_read;
  int     n_groups;
  int    *n_particles_i;
  int    *n_particles_j;
  int     j_read;
  int     mode;
  int     n_groups_i;
  int     n_groups_j;
  int     j_halo;
  int     i_read;
  int     i_read_start;
  int     i_read_stop;
  SID_fp  fp_in;

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  char filename_in[2][MAX_FILENAME_LENGTH];
  char filename_out[MAX_FILENAME_LENGTH];
  strcpy(filename_in[0],argv[1]);
  strcpy(filename_in[1],argv[2]);
  strcpy(filename_out,  argv[3]);
  SID_log("Creating ascii version of {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in);

  // Allocate arrays
  int    **match    =(int    **)SID_malloc(2*sizeof(int *));
  size_t **indices  =(size_t **)SID_malloc(2*sizeof(size_t *));
  float  **score    =(float  **)SID_malloc(2*sizeof(float *));

  for(int i_read=0;i_read<2;i_read++){
     // Read header information
     SID_log("Reading {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in[i_read]);
     int i_read_in;
     int j_read_in;
     int n_groups_1;
     int n_groups_2;
     SID_fopen(filename_in[i_read],"r",&fp_in);
     SID_fread(&i_read_in, sizeof(int),1,&fp_in);SID_log("i_read    =%d",SID_LOG_COMMENT,i_read_in);
     SID_fread(&j_read_in, sizeof(int),1,&fp_in);SID_log("j_read    =%d",SID_LOG_COMMENT,j_read_in);
     SID_fread(&n_groups_1,sizeof(int),1,&fp_in);SID_log("n_groups_i=%d",SID_LOG_COMMENT,n_groups_1);
     SID_fread(&n_groups_2,sizeof(int),1,&fp_in);SID_log("n_groups_j=%d",SID_LOG_COMMENT,n_groups_2);
     if(i_read==0){
        n_groups_i=n_groups_1;
        n_groups_j=n_groups_2;
     }
     else if(n_groups_i!=n_groups_2 || n_groups_j!=n_groups_1)
        SID_trap_error("Input match file halo counts are not symetric (ie %d!=%d or %d!=%d).",ERROR_LOGIC,n_groups_i,n_groups_1,n_groups_j,n_groups_2);

     // Allocate RAM
     match[i_read]    =(int    *)SID_malloc(sizeof(int)   *n_groups_1);
     indices[i_read]  =(size_t *)SID_malloc(sizeof(size_t)*n_groups_1);
     score[i_read]    =(float  *)SID_malloc(sizeof(float) *n_groups_1);

     // Read matches
     for(k_read=0;k_read<n_groups_1;k_read++)
        SID_fread(&(match[i_read][k_read]),sizeof(int),1,&fp_in);

     // Read indices
     for(k_read=0;k_read<n_groups_1;k_read++)
        SID_fread(&(indices[i_read][k_read]),sizeof(size_t),1,&fp_in);

     // Read scores
     for(k_read=0;k_read<n_groups_1;k_read++)
        SID_fread(&(score[i_read][k_read]),sizeof(float),1,&fp_in);

     // Close file
     SID_fclose(&fp_in);
     SID_log("Done.",SID_LOG_CLOSE);
  }

  // Compute 2way match flags
  int *flag_2way=(int *)SID_malloc(sizeof(int)*n_groups_i);
  for(int i_group=0;i_group<n_groups_i;i_group++){
     int match_i=match[0][i_group];
     if(match_i>=0){
        int match_j=match[1][match_i];
        if(match_j==i_group) flag_2way[i_group]=TRUE;
        else                 flag_2way[i_group]=FALSE;
     }
     else
        flag_2way[i_group]=FALSE;
  }

  // Print results
  SID_log("Writing to {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_out);
  FILE *fp_out=fopen(filename_out,"w");
  fprintf(fp_out,"# Ascii dump of {%s}\n",filename_in);
  int i_column=1;
  fprintf(fp_out,"#\n");
  fprintf(fp_out,"# Column (%02d): Halo index in catalog A\n",i_column++);
  fprintf(fp_out,"#        (%02d): Halo match in catalog B\n",i_column++);
  fprintf(fp_out,"#        (%02d): Halo match sort index\n",  i_column++);
  fprintf(fp_out,"#        (%02d): Halo match score\n",       i_column++);
  fprintf(fp_out,"#        (%02d): 2-way match?\n",           i_column++);
  for(k_read=0;k_read<n_groups_i;k_read++)
     fprintf(fp_out,"%7d %7d %7lld %10.3le %d\n",k_read,match[0][k_read],indices[0][k_read],score[0][k_read],flag_2way[k_read]);
  fclose(fp_out);
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  for(int i_read=0;i_read<2;i_read++){
     SID_free(SID_FARG match[i_read]);
     SID_free(SID_FARG indices[i_read]);
     SID_free(SID_FARG score[i_read]); 
  } 
  SID_free(SID_FARG match);
  SID_free(SID_FARG indices);
  SID_free(SID_FARG score); 
  SID_free(SID_FARG flag_2way);  
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

