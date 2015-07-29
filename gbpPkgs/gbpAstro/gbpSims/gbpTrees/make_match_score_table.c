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
  SID_log("Writing match score table...",SID_LOG_OPEN);

  FILE *fp=stdout;
  fprintf(fp,"# Table of match score values.\n");
  fprintf(fp,"# Match Rank Index=\n",MATCH_SCORE_RANK_INDEX);
  fprintf(fp,"# Column (01): n_particles\n");
  fprintf(fp,"#        (02): max score       (n_particles)\n");
  fprintf(fp,"#        (03): min match score (n_particles)\n");
  fprintf(fp,"#        (04): max score       (0.90*n_particles)\n");
  fprintf(fp,"#        (05): max score       (0.80*n_particles)\n");
  fprintf(fp,"#        (06): max score       (0.70*n_particles)\n");
  fprintf(fp,"#        (07): max score       (0.60*n_particles)\n");
  fprintf(fp,"#        (08): max score       (0.50*n_particles)\n");
  fprintf(fp,"#        (09): max score       (0.40*n_particles)\n");
  fprintf(fp,"#        (10): max score       (0.30*n_particles)\n");
  fprintf(fp,"#        (11): max score       (0.20*n_particles)\n");
  fprintf(fp,"#        (12): max score       (0.10*n_particles)\n");
  fprintf(fp,"#        (13): max score       (0.50*n_particles)\n");
  fprintf(fp,"#        (14): max score       (0.01*n_particles)\n");
  for(double n_particles=1.;n_particles<1e9;n_particles*=1.2){
     fprintf(fp,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
             n_particles,
             maximum_match_score(n_particles),
             minimum_match_score(n_particles),
             maximum_match_score(0.9*n_particles),
             maximum_match_score(0.8*n_particles),
             maximum_match_score(0.7*n_particles),
             maximum_match_score(0.6*n_particles),
             maximum_match_score(0.5*n_particles),
             maximum_match_score(0.4*n_particles),
             maximum_match_score(0.3*n_particles),
             maximum_match_score(0.2*n_particles),
             maximum_match_score(0.1*n_particles),
             maximum_match_score(0.05*n_particles),
             maximum_match_score(0.01*n_particles));
  }
  if(fp!=stdout && fp!=stderr)
     fclose(fp);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

