#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <gbpLib.h>
#include <gbpSPH.h>

void free_types(char ***pname,int n_species){
  int     i;
  for(i=0;i<n_species;i++)
    SID_free(SID_FARG (*pname)[i]);
  SID_free(SID_FARG (*pname));
}
