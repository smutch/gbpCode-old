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

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  char filename_SSimPL_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_root,argv[1]);
  SID_log("Checking the integrity of the match files for {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_SSimPL_root);

  // Convert filename_root to filename
  switch(mode){
     case MATCH_SUBGROUPS:
     sprintf(group_text_prefix,"sub");
     break;
     case MATCH_GROUPS:
     sprintf(group_text_prefix,"");
     break;
  }

  // Set the standard SSiMPL match file path
  char filename_root_in[MAX_FILENAME_LENGTH];
  sprintf(filename_root_in,"%s/trees/matches/",filename_SSimPL_root);

  // Read header information
  int i_read_in;
  int j_read_in;
  int n_groups_1;
  int n_groups_2;
  sprintf(filename_in,"%s/%03d/%sgroup_matches_%03d_%03d.dat",filename_root_in,i_read,group_text_prefix,i_read,j_read);
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

  // Read halo sizes from header file
  SID_log("Reading halo sizes...",SID_LOG_OPEN);
  sprintf(filename_in,"%s/%sgroup_matches_header.dat",filename_root_in,group_text_prefix);
  SID_fopen(filename_in,"r",&fp_in);
  SID_fread(&i_read_start,sizeof(int),1,&fp_in);SID_log("snap start  =%d",SID_LOG_COMMENT,i_read_start);
  SID_fread(&i_read_stop, sizeof(int),1,&fp_in);SID_log("snap stop   =%d",SID_LOG_COMMENT,i_read_stop);
  SID_fread(&n_search,    sizeof(int),1,&fp_in);SID_log("search range=%d",SID_LOG_COMMENT,n_search);
  SID_fread(&n_files,     sizeof(int),1,&fp_in);SID_log("# of files  =%d",SID_LOG_COMMENT,n_files);
  int *n_p=(int *)SID_malloc(sizeof(int)*n_groups_i);
  for(k_read=0,max_n_groups=0;k_read<MAX(n_files,i_read);k_read++){
     SID_fread(&l_read,  sizeof(int),1,&fp_in);
     SID_fread(&n_groups,sizeof(int),1,&fp_in);
     if(l_read==i_read){
        if(n_groups!=n_groups_i)
           SID_trap_error("Incorrect number of halos to be read (ie %d!=%d).",ERROR_LOGIC,n_groups,n_groups_i);
        SID_fread(n_p,sizeof(int),n_groups,&fp_in);
     }
     else
        SID_fseek(&fp_in,sizeof(int),n_groups,SID_SEEK_CUR);
     if(mode==MATCH_GROUPS)
        SID_fseek(&fp_in,sizeof(int),n_groups,SID_SEEK_CUR);
  }
  SID_fclose(&fp_in);
  SID_log("Done.",SID_LOG_CLOSE);
  
  // Print results
  for(k_read=0;k_read<n_groups_i;k_read++)
     printf("%7d %7d %7d %7lld %10.3le\n",k_read,n_p[k_read],match[k_read],indices[k_read],score[k_read]);

  // Clean-up
  SID_free(SID_FARG match);
  SID_free(SID_FARG indices);
  SID_free(SID_FARG score);  
  SID_free(SID_FARG n_p);  
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

