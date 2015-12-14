#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

int check_validity_of_match(int n_particles_i,int n_particles_j,float match_score,double f_goodness_of_match){
   // To deal with cases where particle counts suddenly change (during 
   //    halo fragmentation or when the identity of the central subhalo
   //    changes), we use the smaller of the two particle counts.
   int n_particles_use=MIN(n_particles_i,n_particles_j);
   return(match_score>=minimum_match_score((double)n_particles_use));
}
