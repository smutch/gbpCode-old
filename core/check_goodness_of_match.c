#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

#define READ_MATCHES_GOODNESS_A   1.4421650782960167E+01
#define READ_MATCHES_GOODNESS_B   2.7010890377990435E+01
#define READ_MATCHES_GOODNESS_C   ONE_THIRD
#define READ_MATCHES_GOODNESS_D  -2.4490239433470165E+00

int check_goodness_of_match(int n_particles,float match_score){
   float  min_score=(float)pow(READ_MATCHES_GOODNESS_A+READ_MATCHES_GOODNESS_B*READ_MATCHES_GOODNESS_FS*(double)n_particles,READ_MATCHES_GOODNESS_C)+READ_MATCHES_GOODNESS_D;
   if(match_score<min_score)
      return(FALSE);
   else
      return(TRUE);
}

