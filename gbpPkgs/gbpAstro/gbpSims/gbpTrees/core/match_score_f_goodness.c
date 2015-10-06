#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

float match_score_f_goodness(float match_score,int n_particles_in){
   // If score is <=0, return zero ...
   if(match_score<=0.)
      return(0.);
   double n_particles  =(double)n_particles_in;
   float  max_score    =maximum_match_score(n_particles);
   // ... else if score is >= the maximum score, return one ...
   if(match_score>=max_score)
      return(1.);
   // ... else, use a bisection algorythm to find the f_goodness implied by
   //    the given match score and particle count
   float f_goodness_hi=1.;
   float f_goodness_lo=0.;
   float f_goodness;
   float threshold =0.005; 
   float difference=threshold;
   int   iterations=0;
   int   n_max     =1000;
   do{ 
      f_goodness=0.5*(f_goodness_lo+f_goodness_hi);
      float test_score=maximum_match_score(f_goodness*n_particles);
      difference=(match_score-test_score)/match_score;
      if(difference<0.){
         f_goodness_hi=f_goodness;
         difference*=-1;
      }
      else
         f_goodness_lo=f_goodness;
      iterations++;
   } while (difference>=threshold && (++iterations)<n_max);
   return(f_goodness);
}

