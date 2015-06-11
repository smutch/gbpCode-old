#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

int check_goodness_of_match(int n_particles,float match_score,double f_goodness_of_match){
   if(match_score>=F_MAX_MATCH_SCORE_MIN*maximum_match_score((double)n_particles))
      return(TRUE);
   if(match_score<maximum_match_score(f_goodness_of_match*(double)n_particles))
      return(FALSE);
   else
      return(TRUE);
}

