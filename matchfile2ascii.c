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
  char    filename_in[MAX_FILENAME_LENGTH];
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

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  strcpy(filename_in,argv[1]);

  // Read header information
  int i_read_in;
  int j_read_in;
  int n_groups_1;
  int n_groups_2;
  SID_log("Reading header information...",SID_LOG_OPEN);
  SID_fopen(filename_in,"r",&fp_in);
  SID_fread(&i_read_in, sizeof(int),1,&fp_in);SID_log("i_read    =%d",SID_LOG_COMMENT,i_read_in);
  SID_fread(&j_read_in, sizeof(int),1,&fp_in);SID_log("j_read    =%d",SID_LOG_COMMENT,j_read_in);
  SID_fread(&n_groups_i,sizeof(int),1,&fp_in);SID_log("n_groups_i=%d",SID_LOG_COMMENT,n_groups_i);
  SID_fread(&n_groups_j,sizeof(int),1,&fp_in);SID_log("n_groups_j=%d",SID_LOG_COMMENT,n_groups_j);

  // Read matches
  int *match;
  match=(int *)SID_malloc(sizeof(int)*n_groups_i);
  for(k_read=0;k_read<n_groups_i;k_read++)
     SID_fread(&(match[k_read]),sizeof(int),1,&fp_in);

  // Read indices
  size_t *indices;
  indices=(size_t *)SID_malloc(sizeof(size_t)*n_groups_i);
  for(k_read=0;k_read<n_groups_i;k_read++)
     SID_fread(&(indices[k_read]),sizeof(size_t),1,&fp_in);

  // Read scores
  float *score;
  score=(float *)SID_malloc(sizeof(float)*n_groups_i);
  for(k_read=0;k_read<n_groups_i;k_read++)
     SID_fread(&(score[k_read]),sizeof(float),1,&fp_in);

  // Close file
  SID_fclose(&fp_in);
  
  // Print results
  for(k_read=0;k_read<n_groups_i;k_read++)
     printf("%7d %7d %7d %10.3le\n",k_read,match[k_read],indices[k_read],score[k_read]);

  // Clean-up
  SID_free(SID_FARG match);
  SID_free(SID_FARG indices);
  SID_free(SID_FARG score);  
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

