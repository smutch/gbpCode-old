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

void free_treenode_hist(treenode_hist_info **hist){
  SID_free  (SID_FARG (*hist)->array);
  SID_free  (SID_FARG (*hist));
}

