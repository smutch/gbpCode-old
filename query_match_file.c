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
  char    filename_root_in[MAX_FILENAME_LENGTH];
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
  int    *match_ids;
  float  *match_score;
  size_t *match_index;
  int     j_halo;
  int     match;
  int     i_read;
  int     i_read_start;
  int     i_read_stop;
  SID_fp  fp_in;

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  strcpy(filename_root_in,argv[1]);
  if(!strcmp(argv[2],"groups") || !strcmp(argv[2],"group"))
     mode=MATCH_GROUPS;
  else if(!strcmp(argv[2],"subgroups") || !strcmp(argv[2],"subgroup"))
     mode=MATCH_SUBGROUPS;
  else{
     SID_log("Invalid mode selection {%s}.  Should be 'group' or 'subgroup'.",SID_LOG_COMMENT,argv[2]);
     SID_exit(ERROR_SYNTAX);
  }
  i_read=atoi(argv[3]);
  i_halo=atoi(argv[4]);
  SID_log("Searching match information for halo #%d in file #%d from {%s}...",SID_LOG_OPEN,i_halo,i_read,filename_root_in);

  // Convert filename_root to filename
  switch(mode){
     case MATCH_SUBGROUPS:
     sprintf(group_text_prefix,"sub");
     break;
     case MATCH_GROUPS:
     sprintf(group_text_prefix,"");
     break;
  }
  char filename_base[MAX_FILENAME_LENGTH];
  strcpy(filename_base,filename_root_in);
  strip_path(filename_base);

  // Read header information
  SID_log("Reading header information...",SID_LOG_OPEN);
  sprintf(filename_in,"%s/%sgroup_matches_header.dat",filename_root_in,group_text_prefix);
  SID_fopen(filename_in,"r",&fp_in);
  SID_fread(&i_read_start,sizeof(int),1,&fp_in);SID_log("snap start  =%d",SID_LOG_COMMENT,i_read_start);
  SID_fread(&i_read_stop, sizeof(int),1,&fp_in);SID_log("snap stop   =%d",SID_LOG_COMMENT,i_read_stop);
  SID_fread(&n_search,    sizeof(int),1,&fp_in);SID_log("search range=%d",SID_LOG_COMMENT,n_search);
  SID_fread(&n_files,     sizeof(int),1,&fp_in);SID_log("# of files  =%d",SID_LOG_COMMENT,n_files);
  for(k_read=0,max_n_groups=0;k_read<n_files;k_read++){
     SID_fread(&l_read,  sizeof(int),1,&fp_in);
     SID_fread(&n_groups,sizeof(int),1,&fp_in);
     SID_fseek(&fp_in,   sizeof(int),n_groups,SID_SEEK_CUR);
     if(mode==MATCH_GROUPS)
        SID_fseek(&fp_in,   sizeof(int),n_groups,SID_SEEK_CUR);
     max_n_groups=MAX(max_n_groups,n_groups);
  }
  SID_log("Max # groups=%d",SID_LOG_COMMENT,max_n_groups);
  SID_fclose(&fp_in);
  SID_log("Done.",SID_LOG_CLOSE);

  // Initialize some arrays
  n_particles_i=(int    *)SID_malloc(sizeof(int)   *max_n_groups);
  n_particles_j=(int    *)SID_malloc(sizeof(int)   *max_n_groups);
  match_ids    =(int    *)SID_malloc(sizeof(int)   *max_n_groups);
  match_index  =(size_t *)SID_malloc(sizeof(size_t)*max_n_groups);
  match_score  =(float  *)SID_malloc(sizeof(float) *max_n_groups);

  // Loop over all matching combinations
  SID_log("Processing forward matches...",SID_LOG_COMMENT);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  for(j_read=MIN(i_read_stop,i_read+n_search);j_read>=MAX(0,i_read-n_search);j_read--){
     if(i_read!=j_read){
        read_matches(filename_root_in,
                     i_read,
                     j_read,
                     mode,
                     &n_groups_i,
                     &n_groups_j,
                     n_particles_i,
                     n_particles_j,
                     NULL,
                     NULL,
                     match_ids,
                     match_score,
                     match_index);

        // Write desired information
        match=match_ids[i_halo];
        if(match>=0)
          printf("(%3d,%6d,%4d)->(%3d,%6d,%4d) score=%10.3f index=%zu\n",
                 i_read,i_halo,
                 n_particles_i[i_halo],
                 j_read,match,
                 n_particles_j[match],
                 match_score[i_halo],
                 match_index[i_halo]);
        else
          printf("(%3d,%6d,%4d)->(%3d,%6d,%4d) score=%10.3f index=%zu\n",
                 i_read,i_halo,
                 n_particles_i[i_halo],
                 j_read,match,-1,
                 match_score[i_halo],
                 match_index[i_halo]);
     }    
  }
  SID_log("Done.",SID_LOG_CLOSE);
  SID_log("Processing backwards matches...",SID_LOG_COMMENT);
  j_read=i_read;
  j_halo=i_halo;
  for(i_read=MIN(i_read_stop,j_read+n_search);i_read>=MAX(0,j_read-n_search);i_read--){
     if(i_read!=j_read){
        read_matches(filename_root_in,
                     i_read,
                     j_read,
                     mode,
                     &n_groups_i,
                     &n_groups_j,
                     n_particles_i,
                     n_particles_j,
                     NULL,
                     NULL,
                     match_ids,
                     match_score,
                     match_index);

        // Write desired information
        for(i_halo=0,match=-1;i_halo<n_groups_i && match<0;i_halo++){
          if(match_ids[i_halo]==j_halo){
            printf("(%3d,%6d,%4d)->(%3d,%6d,%4d) score=%10.3f index=%zu\n",
                    i_read,i_halo,
                    n_particles_i[i_halo],
                    j_read,j_halo,
                    n_particles_j[j_halo],
                    match_score[i_halo],
                    match_index[i_halo]);
          }
        }
     }    
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  SID_free(SID_FARG n_particles_i);
  SID_free(SID_FARG n_particles_j);
  SID_free(SID_FARG match_ids);
  SID_free(SID_FARG match_index);
  SID_free(SID_FARG match_score);
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

