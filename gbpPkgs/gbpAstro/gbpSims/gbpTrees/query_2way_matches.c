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
  char    filename_SSimPL_root[MAX_FILENAME_LENGTH];
  char    filename_in[MAX_FILENAME_LENGTH];
  char    group_text_prefix[4];
  int     n_files;
  int     k_read;
  int     max_n_groups;
  int     l_read;
  int     n_groups;
  int     j_read;
  int     mode;
  int     n_groups_i;
  int     n_groups_j;
  int     j_halo;
  int     match;
  int     i_read;
  int     i_read_start;
  int     i_read_stop;
  SID_fp  fp_in;

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  strcpy(filename_SSimPL_root,argv[1]);
  if(!strcmp(argv[2],"groups") || !strcmp(argv[2],"group"))
     mode=MATCH_GROUPS;
  else if(!strcmp(argv[2],"subgroups") || !strcmp(argv[2],"subgroup"))
     mode=MATCH_SUBGROUPS;
  else{
     SID_log("Invalid mode selection {%s}.  Should be 'group' or 'subgroup'.",SID_LOG_COMMENT,argv[2]);
     SID_exit(ERROR_SYNTAX);
  }
  i_read=atoi(argv[3]);
  j_read=atoi(argv[4]);
  SID_log("Searching match information for halo #%d in file #%d of {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,i_halo,i_read,filename_SSimPL_root);

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

  // Set the output file
  char filename_base[MAX_FILENAME_LENGTH];
  char filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_base,filename_SSimPL_root);
  if(!strcmp(&(filename_base[strlen(filename_base)-1]),"/"))
     strcpy(&(filename_base[strlen(filename_base)-1]),"\0");
  strip_path(filename_base);
  sprintf(filename_out,"%s_%d_%d_2way_matches.txt",filename_base,i_read,j_read);

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
  int    *n_particles_i       =(int    *)SID_malloc(sizeof(int)   *max_n_groups);
  int    *n_particles_j       =(int    *)SID_malloc(sizeof(int)   *max_n_groups);
  int    *match_forward_ids   =(int    *)SID_malloc(sizeof(int)   *max_n_groups);
  size_t *match_forward_index =(size_t *)SID_malloc(sizeof(size_t)*max_n_groups);
  float  *match_forward_score =(float  *)SID_malloc(sizeof(float) *max_n_groups);
  char   *match_forward_2way  =(char   *)SID_malloc(sizeof(char)  *max_n_groups);
  int    *match_backward_ids  =(int    *)SID_malloc(sizeof(int)   *max_n_groups);
  size_t *match_backward_index=(size_t *)SID_malloc(sizeof(size_t)*max_n_groups);
  float  *match_backward_score=(float  *)SID_malloc(sizeof(float) *max_n_groups);
  char   *match_backward_2way =(char   *)SID_malloc(sizeof(char)  *max_n_groups);

  // Loop over all matching combinations
  SID_log("Reading forward matches...",SID_LOG_OPEN|SID_LOG_TIMER);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  read_matches(filename_root_in,
               i_read,
               j_read,
               max_n_groups,
               mode,
               &n_groups_i,
               &n_groups_j,
               n_particles_i,
               n_particles_j,
               NULL,
               NULL,
               match_forward_ids,
               match_forward_score,
               match_forward_index,
               match_forward_2way,
               0.);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_log("Processing backwards matches...",SID_LOG_OPEN|SID_LOG_TIMER);
  read_matches(filename_root_in,
               j_read,
               i_read,
               max_n_groups,
               mode,
               &n_groups_j,
               &n_groups_i,
               n_particles_j,
               n_particles_i,
               NULL,
               NULL,
               match_backward_ids,
               match_backward_score,
               match_backward_index,
               match_backward_2way,
               0.);
  SID_log("Done.",SID_LOG_CLOSE);

  // Open output file
  FILE *fp_out;
  fp_out=fopen(filename_out,"w");
  int i_column=1;
  fprintf(fp_out,"# Column (%02d): Halo index for snapshot %d\n",    i_column++,i_read);
  fprintf(fp_out,"#        (%02d): Halo index for snapshot %d\n",    i_column++,j_read);
  fprintf(fp_out,"#        (%02d): No. particles in snapshot %d\n",  i_column++,i_read);
  fprintf(fp_out,"#        (%02d): No. particles in snapshot %d\n",  i_column++,j_read);
  fprintf(fp_out,"#        (%02d): Forward  match score\n",          i_column++);
  fprintf(fp_out,"#        (%02d): Forward  match score/max score\n",i_column++);
  fprintf(fp_out,"#        (%02d): Backward match score\n",          i_column++);
  fprintf(fp_out,"#        (%02d): Backward match score/max score\n",i_column++);
  for(int i_halo=0;i_halo<n_groups_i;i_halo++){
     int j_halo=match_forward_ids[i_halo];
     if(match_forward_2way[i_halo]){
        if(j_halo<0 || j_halo>n_groups_j)
           SID_trap_error("There's an invalid match id (ie. %d<0 || %d>%d)  attached to a 2-way match!",ERROR_LOGIC,j_halo,j_halo,n_groups_j);
        fprintf(fp_out,"%7d %7d %6d %6d %10.3le %10.3le %10.3le %10.3le\n",
                i_halo,
                j_halo,
                n_particles_i[i_halo],
                n_particles_j[j_halo],
                match_forward_score[i_halo],
                match_forward_score[i_halo]/maximum_match_score((double)n_particles_i[i_halo]),
                match_backward_score[j_halo],
                match_backward_score[j_halo]/maximum_match_score((double)n_particles_j[j_halo]));
     }                
  }
  fclose(fp_out);

  // Clean-up
  SID_free(SID_FARG n_particles_i);
  SID_free(SID_FARG n_particles_j);
  SID_free(SID_FARG match_forward_ids);
  SID_free(SID_FARG match_forward_index);
  SID_free(SID_FARG match_forward_score);
  SID_free(SID_FARG match_forward_2way);
  SID_free(SID_FARG match_backward_ids);
  SID_free(SID_FARG match_backward_index);
  SID_free(SID_FARG match_backward_score);
  SID_free(SID_FARG match_backward_2way);
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

