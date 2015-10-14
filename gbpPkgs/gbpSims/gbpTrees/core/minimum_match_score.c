#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

float minimum_match_score(double n_particles){
   float check_1=maximum_match_score(F_GOODNESS_OF_MATCH*(double)n_particles);
   float check_2=F_MAX_MATCH_SCORE_MIN*maximum_match_score((double)n_particles);
   float check_3=MIN_MATCH_SCORE;
   return(MAX(MAX(check_1,check_2),check_3));
}

