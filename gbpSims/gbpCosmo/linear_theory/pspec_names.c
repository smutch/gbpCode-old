#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

void pspec_names(int   mode,
       int   component,
       char *mode_name,
       char *component_name){
  switch(mode){
  case PSPEC_LINEAR_TF:
    sprintf(mode_name,"TF");
    break;
  case PSPEC_LINEAR_BBKS:
    sprintf(mode_name,"BBKS");
    break;
  case PSPEC_NONLINEAR_JAIN:
    sprintf(mode_name,"JAIN");
    break;
  case PSPEC_NONLINEAR_PD:
    sprintf(mode_name,"NL_PD");
    break;
  case PSPEC_NONLINEAR_SMITH:
    sprintf(mode_name,"NL_Smith");
    break;
  default:
    fprintf(stderr,"Invalid PSPEC mode in pspec_names!\n");
    break;
  }
  switch(component){
  case PSPEC_ALL_MATTER:
    sprintf(component_name,"all");
    break;
  case PSPEC_DARK_MATTER:
    sprintf(component_name,"dark");
    break;
  case PSPEC_BARYON:
    sprintf(component_name,"gas");
    break;
  default:
    fprintf(stderr,"Invalid PSPEC component in pspec_names!\n");
    break;
  }
}

