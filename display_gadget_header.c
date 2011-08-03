#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpSPH.h>

void display_gadget_header(plist_info  *plist){
  FILE   *fp;
  char  **pname;
  int     i,j,k; 
  int     counter;
  size_t  n_of_type[N_GADGET_TYPE];
  int     n_of_type_tmp[N_GADGET_TYPE];
  int     flag_used[N_GADGET_TYPE];
  size_t  n_particles;
  int     unused[256];
  size_t  n_all[N_GADGET_TYPE];
  int     n_all_tmp[N_GADGET_TYPE];
  int     n_files;
  int     junk;
  double  d_value;
  double *d_array;
  double *d1_array;
  double *d2_array;
  double *d3_array;
  int     min_i;
  int     max_i;
  double  mean;
  double  min;
  double  max;
  double  median;
  double  std_dev;
  double  mass_array;
  float   f_temp;
  float   f1_temp;
  float   f2_temp;
  float   f3_temp;
  int     i_value;
  int    *i_array;
  int     i_temp;
  int     i1_temp;
  int     i2_temp;
  int     n_type_used;
  long    record_length;
  int     n_return;
  int     s_load;
  int     flag_alloc_d1_array;
  size_t  n_all_species;
  char    var_name[256];
  char    var_name2[256];

  double  redshift;
  double  h_Hubble;
  double  Omega_M;
  double  Omega_Lambda;
  double  rho_crit;

  pname=plist->species;
  
  // Determine which species are present
  n_particles=0;
  for(i=0,n_type_used=0;i<plist->n_species;i++) {
    if(ADaPS_exist(plist->data,"n_%s",pname[i])){
      n_of_type[i]=((size_t *)ADaPS_fetch(plist->data,"n_%s",pname[i]))[0];
      if(n_of_type[i]>0){
        n_particles+=n_of_type[i];
        flag_used[i]=TRUE;
        n_type_used++;
      }
      else{
        n_of_type[i]=0;
        flag_used[i]=FALSE;
      }
    }
    else{
      n_of_type[i]=0;
      flag_used[i]=FALSE;
    }
  }

  h_Hubble=((double *)ADaPS_fetch(plist->data,"h_Hubble"))[0];
  fprintf(stderr,"\n");
  
  // Expansion factor (or time)
  if(ADaPS_exist(plist->data,"expansion_factor"))
    fprintf(stderr,"%20s = %le\n","Expansion factor",((double *)ADaPS_fetch(plist->data,"expansion_factor"))[0]);
  else if(ADaPS_exist(plist->data,"time"))
    fprintf(stderr,"%20s = %le Myrs\n","Time",((double *)ADaPS_fetch(plist->data,"time"))[0]/(S_PER_MYR));
  else
    fprintf(stderr,"time/expansion factor not set!\n");

  // Redshift
  if(ADaPS_exist(plist->data,"redshift")){
    redshift=((double *)ADaPS_fetch(plist->data,"redshift"))[0];
    fprintf(stderr,"%20s = %le\n","Redshift",redshift);
  }
  else{
    fprintf(stderr,"redshift not set!\n");
    redshift=0.;
  }

  // Some flags
  if(ADaPS_exist(plist->data,"flag_Sfr"))
    fprintf(stderr,"%20s = %d\n","flag_Sfr",((int *)ADaPS_fetch(plist->data,"flag_Sfr"))[0]);
  else
    fprintf(stderr,"flag_Sfr not set!\n");
  if(ADaPS_exist(plist->data,"flag_feedback"))
    fprintf(stderr,"%20s = %d\n","flag_feedback",((int *)ADaPS_fetch(plist->data,"flag_feedback"))[0]);
  else
    fprintf(stderr,"flag_feedback not set!\n");

  // Another flag
  if(ADaPS_exist(plist->data,"flag_cooling"))
    fprintf(stderr,"%20s = %d\n","flag_cooling",((int *)ADaPS_fetch(plist->data,"flag_cooling"))[0]);
  else
    fprintf(stderr,"flag_cooling not set!\n");

  // Number of files per snapshot
  if(ADaPS_exist(plist->data,"n_files"))
    fprintf(stderr,"%20s = %d\n","files per snap",((int *)ADaPS_fetch(plist->data,"n_files"))[0]);
  else
    fprintf(stderr,"files per snapshot not set!\n");

  // Box size
  if(ADaPS_exist(plist->data,"box_size"))
    fprintf(stderr,"%20s = %le [kpc/h]\n","box size",((double *)ADaPS_fetch(plist->data,"box_size"))[0]/(M_PER_KPC/h_Hubble));
  else
    fprintf(stderr,"box size not set!\n");

  // Cosmology

  if(ADaPS_exist(plist->data,"h_Hubble")){
    fprintf(stderr,"%20s = %le\n","h_Hubble",h_Hubble);
  }
  else{
    h_Hubble=1.;
    fprintf(stderr,"h_Hubble not set!\n");
  }

  if(ADaPS_exist(plist->data,"Omega_M")){
    Omega_M=((double *)ADaPS_fetch(plist->data,"Omega_M"))[0];
    fprintf(stderr,"%20s = %le\n","Omega_M",Omega_M);
  }
  else
    fprintf(stderr,"Omega_M not set!\n");
  if(Omega_M<=0.) Omega_M=0.3;

  if(ADaPS_exist(plist->data,"Omega_Lambda")){
    Omega_Lambda=((double *)ADaPS_fetch(plist->data,"Omega_Lambda"))[0];
    fprintf(stderr,"%20s = %le\n","Omega_Lambda",Omega_Lambda);
  }
  else
    fprintf(stderr,"Omega_Lambda not set!\n");

  /* Number of particles */
  for(i=0;i<plist->n_species;i++){
    sprintf(var_name,"n_%s",pname[i]);
    if(ADaPS_exist(plist->data,var_name)){
      sprintf(var_name2,"n_[%s]",pname[i]);
      fprintf(stderr,"%20s = %zd\n",var_name2,((size_t *)ADaPS_fetch(plist->data,var_name))[0]);
    }
  }
  for(i=0;i<plist->n_species;i++){
    sprintf(var_name,"n_all_%s",pname[i]);
    if(ADaPS_exist(plist->data,var_name)){
      sprintf(var_name2,"n_all_[%s]",pname[i]);
      fprintf(stderr,"%20s = %zd\n",var_name2,((size_t *)ADaPS_fetch(plist->data,var_name))[0]);
    }
  }
    
  /* Particle mass array */
  for(i=0;i<plist->n_species;i++){
    if(flag_used[i]){
      sprintf(var_name,"mass_array_%s",pname[i]);
      if(ADaPS_exist(plist->data,var_name)){
        sprintf(var_name2,"mass_array_[%s]",pname[i]);
        fprintf(stderr,"%20s = %le [M_sol/h]\n",var_name2,((double *)ADaPS_fetch(plist->data,var_name))[0]/(M_SOL/h_Hubble));
      }
    }
  }

};
