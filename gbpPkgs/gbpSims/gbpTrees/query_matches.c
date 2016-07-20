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
  int    *n_particles_i;
  int    *n_particles_j;
  int     j_read;
  int     mode;
  int     n_groups_i;
  int     n_groups_j;
  int    *match_ids;
  float  *match_score;
  size_t *match_index;
  char   *match_2way;
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
  i_halo=atoi(argv[4]);
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
  sprintf(filename_out,"%s_%d_%d_matches.txt",filename_base,i_read,i_halo);

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
  match_2way   =(char   *)SID_malloc(sizeof(char)  *max_n_groups);

  // Open output file
  FILE *fp_out;
  fp_out=fopen(filename_out,"w");
  int i_column=1;

  fprintf(fp_out,"# Parameters dicating goodness of fit:\n");
  fprintf(fp_out,"#    F_GOODNESS_OF_MATCH  =%le\n",F_GOODNESS_OF_MATCH);
  fprintf(fp_out,"#    F_MAX_MATCH_SCORE_MIN=%le\n",F_MAX_MATCH_SCORE_MIN);
  fprintf(fp_out,"#    MIN_MATCH_SCORE      =%le\n",MIN_MATCH_SCORE);
  fprintf(fp_out,"# Column (%02d): Snapshot\n",                i_column++);
  fprintf(fp_out,"#        (%02d): Halo index\n",              i_column++);
  fprintf(fp_out,"#        (%02d): No. particles\n",           i_column++);
  fprintf(fp_out,"#        (%02d): Matched to snapshot\n",     i_column++);
  fprintf(fp_out,"#        (%02d): Matched to index\n",        i_column++);
  fprintf(fp_out,"#        (%02d): Matched to No. particles\n",i_column++);
  fprintf(fp_out,"#        (%02d): Match score\n",             i_column++);
  fprintf(fp_out,"#        (%02d): Match score f_goodness\n",  i_column++);
  fprintf(fp_out,"#        (%02d): Match sort index\n",        i_column++);
  fprintf(fp_out,"#        (%02d): 2-way or 1-way match?\n",   i_column++);
  fprintf(fp_out,"#        (%02d): Good or bad match?\n",      i_column++);

  // Loop over all matching combinations
  SID_log("Processing forward matches...",SID_LOG_OPEN|SID_LOG_TIMER);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  for(j_read=MIN(i_read_stop,i_read+n_search);j_read>=MAX(0,i_read-n_search);j_read--){
     if(i_read!=j_read){
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
                     match_ids,
                     match_score,
                     match_index,
                     match_2way,
                     0.);

        // Check for goodness of match
        char goodness_of_match_text[5];
        int n_particles_j_i=0;
        if(match_ids[i_halo]>=0)
           n_particles_j_i=n_particles_j[match_ids[i_halo]];
        if(check_validity_of_match(n_particles_i[i_halo],match_score[i_halo],F_GOODNESS_OF_MATCH))
           sprintf(goodness_of_match_text,"good");
        else
           sprintf(goodness_of_match_text,"bad");

        // Check for 2-way match
        char twoway_match_text[5];
        if(match_2way[i_halo])
           sprintf(twoway_match_text,"2-way");
        else
           sprintf(twoway_match_text,"1-way");

        // Write desired information
        match=match_ids[i_halo];
        if(match>=0)
          fprintf(fp_out,"%3d %7d %6d %3d %7d %6d %10.3le %10.3le %7zu %s %s\n",
                         i_read,i_halo,
                         n_particles_i[i_halo],
                         j_read,match,
                         n_particles_j[match],
                         match_score[i_halo],
                         match_score_f_goodness(match_score[i_halo],n_particles_i[i_halo]),
                         match_index[i_halo],
                         twoway_match_text,
                         goodness_of_match_text);
        else
          fprintf(fp_out,"%3d %7d %6d %3d %7d %6d %10.3le %10.3le %7zu %s %s\n",
                         i_read,i_halo,
                         n_particles_i[i_halo],
                         j_read,match,-1,
                         match_score[i_halo],
                         match_score_f_goodness(match_score[i_halo],n_particles_i[i_halo]),
                         match_index[i_halo],
                         twoway_match_text,
                         goodness_of_match_text);
     }    
  }
  SID_log("Done.",SID_LOG_CLOSE);
  SID_log("Processing backwards matches...",SID_LOG_OPEN|SID_LOG_TIMER);
  j_read=i_read;
  j_halo=i_halo;
  for(i_read=MIN(i_read_stop,j_read+n_search);i_read>=MAX(0,j_read-n_search);i_read--){
     if(i_read!=j_read){
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
                     match_ids,
                     match_score,
                     match_index,
                     match_2way,
                     0.);

        // Write desired information
        for(i_halo=0,match=-1;i_halo<n_groups_i && match<0;i_halo++){
          if(match_ids[i_halo]==j_halo){

            // Check for goodness of match
            char goodness_of_match_text[5];
            int n_particles_j_i=0;
            if(match_ids[i_halo]>=0)
               n_particles_j_i=n_particles_j[j_halo];
            if(check_validity_of_match(n_particles_i[i_halo],match_score[i_halo],F_GOODNESS_OF_MATCH))
               sprintf(goodness_of_match_text,"good");
            else
               sprintf(goodness_of_match_text,"bad");

            // Check for 2-way match
            char twoway_match_text[5];
            if(match_2way[i_halo])
               sprintf(twoway_match_text,"2-way");
            else
               sprintf(twoway_match_text,"1-way");

            fprintf(fp_out,"%3d %7d %6d %3d %7d %6d %10.3le %10.3le %7zu %s %s\n",
                           i_read,i_halo,
                           n_particles_i[i_halo],
                           j_read,j_halo,
                           n_particles_j[j_halo],
                           match_score[i_halo],
                           match_score_f_goodness(match_score[i_halo],n_particles_i[i_halo]),
                           match_index[i_halo],
                           twoway_match_text,
                           goodness_of_match_text);
          }
        }
     }    
  }
  SID_log("Done.",SID_LOG_CLOSE);
  SID_log("Output written to {%s}",SID_LOG_COMMENT,filename_out);

  fclose(fp_out);

  // Clean-up
  SID_free(SID_FARG n_particles_i);
  SID_free(SID_FARG n_particles_j);
  SID_free(SID_FARG match_ids);
  SID_free(SID_FARG match_index);
  SID_free(SID_FARG match_score);
  SID_free(SID_FARG match_2way);
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

