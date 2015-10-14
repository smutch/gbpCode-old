#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gsl/gsl_sf_zeta.h>

#define READ_MATCHES_GOODNESS_A   1.4421650782960167E+01
#define READ_MATCHES_GOODNESS_B   2.7010890377990435E+01
#define READ_MATCHES_GOODNESS_C   ONE_THIRD
#define READ_MATCHES_GOODNESS_D  -2.4490239433470165E+00

#define MAXIMUM_SCORE_COEFF    2.61   // Only if MATCH_SCORE_RANK_INDEX=-1.5.  This comes from Wolfram Alpha.  Query "sum k^(-1.5) from k = 1 to N"
#define EULER_MASCHERONI_CONST 0.5772 // Used if MATCH_SCORE_RANK_INDEX=-1

float maximum_match_score(double n_particles){
   if(n_particles<1)
      return(0);
   if(MATCH_SCORE_RANK_INDEX==(-1.5))
      return(MAXIMUM_SCORE_COEFF-gsl_sf_hzeta(-MATCH_SCORE_RANK_INDEX,n_particles));
   else if(MATCH_SCORE_RANK_INDEX==(-1.))
      return(EULER_MASCHERONI_CONST+take_ln(n_particles));
   else if(MATCH_SCORE_RANK_INDEX==(-TWO_THIRDS))
      return((float)pow(READ_MATCHES_GOODNESS_A+READ_MATCHES_GOODNESS_B*n_particles,READ_MATCHES_GOODNESS_C)+READ_MATCHES_GOODNESS_D);
   else
      SID_trap_error("Invalid MATCH_SCORE_RANK_INDEX passed to maximum_match_score().  Please update the code to accomodate MATCH_SCORE_RANK_INDEX=%lf.",
                     ERROR_LOGIC,MATCH_SCORE_RANK_INDEX);
}

