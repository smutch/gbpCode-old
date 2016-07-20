#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  float match_weight_rank_index=MATCH_SCORE_RANK_INDEX;
  char filename_halo_root_in_1[256];
  char filename_halo_root_in_2[256];
  char filename_root_matches[256];
  int  i_read;
  strcpy(filename_halo_root_in_1,argv[1]);
  strcpy(filename_halo_root_in_2,argv[2]);
  strcpy(filename_root_matches,  argv[3]);
  i_read                   =atoi(argv[4]);
  if(argc==6)
     match_weight_rank_index=atof(argv[5]);


  SID_log("Constructing matches between catalogs {%s} and {%s} for snapshot #%d (weight index=%f)...",SID_LOG_OPEN|SID_LOG_TIMER,
          filename_halo_root_in_1,
          filename_halo_root_in_2,
          i_read,
          match_weight_rank_index);

  // Perform matching
  compute_cross_catalog_matches(filename_halo_root_in_1,
                                filename_halo_root_in_2,
                                filename_root_matches,
                                i_read,
                                i_read,
                                match_weight_rank_index);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

