#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){
  SID_init(&argc,&argv,NULL,NULL);

  // Parse command line
  char filename_from[MAX_FILENAME_LENGTH];
  char filename_to[MAX_FILENAME_LENGTH];
  strcpy(filename_from,     argv[1]);
  int column_from   =  atoi(argv[2]);
  strcpy(filename_to,       argv[3]);
  int column_to     =  atoi(argv[4]);

  SID_log("Computing match score from column #%d in {%s} to column #%d in {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,column_from,filename_from,column_to,filename_to);

  // Read 'from' IDs
  SID_log("Reading 'from' catalog...",SID_LOG_OPEN|SID_LOG_TIMER);
  FILE   *fp_from    =fopen(filename_from,"r");
  size_t  line_length=0;
  char   *line       =NULL;
  int     n_ids_from =count_lines_data(fp_from);
  SID_log("(%d items)...",SID_LOG_CONTINUE,n_ids_from);
  size_t *ids_from   =(size_t *)SID_malloc(sizeof(size_t)*n_ids_from);
  for (int i_from=0;i_from<n_ids_from;i_from++){
     grab_next_line_data(fp_from,&line,&line_length);
     grab_size_t(line,column_from,&(ids_from[i_from]));
  }
  fclose(fp_from);
  SID_log("Done.",SID_LOG_CLOSE);

  // Read 'to' IDs
  SID_log("Reading 'to' catalog...",SID_LOG_OPEN|SID_LOG_TIMER);
  FILE   *fp_to      =fopen(filename_to,"r");
  int     n_ids_to   =count_lines_data(fp_to);
  SID_log("(%d items)...",SID_LOG_CONTINUE,n_ids_to);
  size_t *ids_to     =(size_t *)SID_malloc(sizeof(size_t)*n_ids_to);
  for (int i_to=0;i_to<n_ids_to;i_to++){
     grab_next_line_data(fp_to,&line,&line_length);
     grab_size_t(line,column_to,&(ids_to[i_to]));
  }
  fclose(fp_to);
  SID_log("Done.",SID_LOG_CLOSE);

  // Perform sorting
  SID_log("Sorting...",SID_LOG_OPEN|SID_LOG_TIMER);
  size_t *ids_from_index=NULL;
  size_t *ids_to_index  =NULL;
  merge_sort(ids_from,n_ids_from,&ids_from_index,SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  merge_sort(ids_to,  n_ids_to,  &ids_to_index,  SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  SID_log("Done.",SID_LOG_CLOSE);

  // Compute and report score
  double score=0;
  for(int i_from=0,i_to=0;i_from<n_ids_from;i_from++){
     size_t idx_from=ids_from_index[i_from];
     while(ids_to[ids_to_index[i_to]]<ids_from[idx_from] && i_to<(n_ids_to-1)) i_to++;
     if(ids_to[ids_to_index[i_to]]==ids_from[idx_from])
        score+=pow((double)(idx_from+1),MATCH_SCORE_RANK_INDEX);
  }
  SID_log("Score = %le (%.3lf x goodness criterion)",SID_LOG_COMMENT,score,score/maximum_match_score(F_GOODNESS_OF_MATCH*(double)n_ids_from));

  // Clean-up
  SID_free(SID_FARG ids_from);
  SID_free(SID_FARG ids_to);
  SID_free(SID_FARG ids_from_index);
  SID_free(SID_FARG ids_to_index);
  SID_free(SID_FARG line);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

