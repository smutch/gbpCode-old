#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double W_k_tophat(double kR){
  return(3.*(sin(kR)-kR*cos(kR))/pow(kR,3.));
}

