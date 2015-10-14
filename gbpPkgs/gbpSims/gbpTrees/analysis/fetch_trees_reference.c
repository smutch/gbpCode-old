#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

tree_info *fetch_trees_reference(tree_info *trees){
   if(trees->trees_reference!=NULL)
      return(trees->trees_reference);
   else
      return(trees);
}

