#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

void free_treenode_list(treenode_list_info **list){

  SID_log("Freeing treenode list...",SID_LOG_OPEN);

  SID_free  (SID_FARG (*list)->list);
  ADaPS_free(SID_FARG (*list)->data);
  SID_free  (SID_FARG (*list));
 
  SID_log("Done.",SID_LOG_CLOSE);
}

