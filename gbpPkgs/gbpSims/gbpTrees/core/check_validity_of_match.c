#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

int check_validity_of_match(int n_particles_use,float match_score,double f_goodness_of_match){
   return(match_score>=minimum_match_score((double)n_particles_use,f_goodness_of_match));
}
