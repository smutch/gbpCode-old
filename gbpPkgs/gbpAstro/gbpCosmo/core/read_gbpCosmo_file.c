#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

void read_gbpCosmo_file(cosmo_info **cosmo,const char *filename_in){
  char   Name[PARAMETER_STRING_LENGTH];
  double Omega_M;
  double Omega_k;
  double Omega_Lambda;
  double Omega_b;
  double f_gas;
  double h_Hubble;
  double sigma_8;
  double n_spectral;

  // Define the parameter file
  parameter_list_info *parameter_list=NULL;
  init_parameter_list(&parameter_list);
  add_parameter_to_list(parameter_list,"Name",        SID_CHAR,  PARAMETER_MODE_OPTIONAL);
  add_parameter_to_list(parameter_list,"Omega_Lambda",SID_DOUBLE,PARAMETER_MODE_DEFAULT);
  add_parameter_to_list(parameter_list,"Omega_M",     SID_DOUBLE,PARAMETER_MODE_DEFAULT);
  add_parameter_to_list(parameter_list,"Omega_k",     SID_DOUBLE,PARAMETER_MODE_DEFAULT);
  add_parameter_to_list(parameter_list,"Omega_b",     SID_DOUBLE,PARAMETER_MODE_DEFAULT);
  add_parameter_to_list(parameter_list,"h_Hubble",    SID_DOUBLE,PARAMETER_MODE_DEFAULT);
  add_parameter_to_list(parameter_list,"sigma_8",     SID_DOUBLE,PARAMETER_MODE_DEFAULT);
  add_parameter_to_list(parameter_list,"n_spectral",  SID_DOUBLE,PARAMETER_MODE_DEFAULT);
  add_parameter_to_list(parameter_list,"f_gas",       SID_DOUBLE,PARAMETER_MODE_OPTIONAL);

  // Read the cosmology from file
  char filename[MAX_FILENAME_LENGTH];
  if(!strcmp(filename_in,"default"))
     sprintf(filename,"%s/default_cosmology.txt",GBP_DATA_DIR);
  else
     sprintf(filename,"%s",filename_in);
  read_parameter_file(filename,parameter_list);
  fetch_parameter_data(parameter_list,"Omega_Lambda",&Omega_Lambda); 
  fetch_parameter_data(parameter_list,"Omega_M",     &Omega_M);
  fetch_parameter_data(parameter_list,"Omega_k",     &Omega_k);
  fetch_parameter_data(parameter_list,"Omega_b",     &Omega_b);
  fetch_parameter_data(parameter_list,"h_Hubble",    &h_Hubble);
  fetch_parameter_data(parameter_list,"sigma_8",     &sigma_8);
  fetch_parameter_data(parameter_list,"n_spectral",  &n_spectral);

  // Set defaults for optional paramaters
  char *Name_pass=NULL;
  if(fetch_parameter_data(parameter_list,"Name",Name))
     Name_pass=Name;
  if(!fetch_parameter_data(parameter_list,"f_gas",&f_gas))
     f_gas=Omega_b/Omega_M;

  // Perform initialization
  init_cosmo(cosmo,
             Name_pass,
             Omega_Lambda,
             Omega_M,
             Omega_k,
             Omega_b,
             f_gas,
             h_Hubble,
             sigma_8,
             n_spectral);

  //Clean-up
  free_parameter_list(&parameter_list);
}

