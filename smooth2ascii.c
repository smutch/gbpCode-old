#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>
int main(int argc, char *argv[]){
  char        filename_out[256];
  char        filename_smooth[256];
  char        filename_snapshot[256];
  char       *species_name;
  double      h_Hubble;
  plist_info  plist;
  size_t      i_particle;
  int         i_species;
  int         j_species;
  int         i_rank;
  size_t      n_particles;
  REAL       *x_array;
  REAL       *y_array;
  REAL       *z_array;
  REAL       *r_smooth_array;
  REAL       *rho_array;
  REAL       *sigma_v_array;
  FILE       *fp_out;

  SID_init(&argc,&argv,NULL);

  strcpy(filename_snapshot,argv[1]);
  strcpy(filename_smooth,  argv[2]);
  strcpy(filename_out,     argv[3]);

  SID_log("Creating ascii file {%s} from smmoth files {%s} and snapshot {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_out,filename_smooth,filename_snapshot);

  // Read snapshot files
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  read_gadget_binary(filename_snapshot,&plist,READ_GADGET_DEFAULT);
  read_smooth(&plist,filename_smooth,0,SMOOTH_DEFAULT);
  h_Hubble=((double *)ADaPS_fetch(plist.data,"h_Hubble",species_name))[0];

  // Loop over each species
  for(i_species=0,j_species=0;i_species<N_GADGET_TYPE;i_species++){
    species_name=plist.species[i_species];
    if(ADaPS_exist(plist.data,"n_all_%s",species_name))
      n_particles=((size_t *)ADaPS_fetch(plist.data,"n_all_%s",species_name))[0];
    else
     n_particles=0;
    // If at least one rank has particles for this species ...
    if(n_particles>0){
    SID_log("Writting %s particles...",SID_LOG_OPEN,species_name);
    // ... then fetch arrays ...
    n_particles=((size_t *)ADaPS_fetch(plist.data,"n_%s",species_name))[0];
    x_array    =(REAL *)ADaPS_fetch(plist.data,"x_%s",species_name);
    y_array    =(REAL *)ADaPS_fetch(plist.data,"y_%s",species_name);
    z_array    =(REAL *)ADaPS_fetch(plist.data,"z_%s",species_name);
    if(ADaPS_exist(plist.data,"r_smooth_%s",species_name))
      r_smooth_array=(REAL *)ADaPS_fetch(plist.data,"r_smooth_%s",species_name);
    else
      r_smooth_array=NULL;
    if(ADaPS_exist(plist.data,"rho_%s",species_name))
      rho_array=(REAL *)ADaPS_fetch(plist.data,"rho_%s",species_name);
    else
      rho_array=NULL;
    if(ADaPS_exist(plist.data,"sigma_v_%s",species_name))
      sigma_v_array=(REAL *)ADaPS_fetch(plist.data,"sigma_v_%s",species_name);
    else
      sigma_v_array=NULL;
    
    // ... and write this species' particles
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
      if(SID.My_rank==i_rank){
        if(j_species==0 && i_rank==0)
          fp_out=fopen(filename_out,"w");
        else
          fp_out=fopen(filename_out,"a");
        for(i_particle=0;i_particle<n_particles;i_particle++){
          fprintf(fp_out,"%2d %11.4le %11.4le %11.4le",
                  i_species,
                  (double)x_array[i_particle]*h_Hubble/M_PER_MPC,
                  (double)y_array[i_particle]*h_Hubble/M_PER_MPC,
                  (double)z_array[i_particle]*h_Hubble/M_PER_MPC);
          if(r_smooth_array!=NULL)
            fprintf(fp_out," %10.4le",(double)r_smooth_array[i_particle]*h_Hubble/M_PER_MPC);
          if(rho_array!=NULL)
            fprintf(fp_out," %10.4le",(double)rho_array[i_particle]/(M_SOL*pow(h_Hubble/M_PER_MPC,3.)));
          if(sigma_v_array!=NULL)
            fprintf(fp_out," %10.4le",(double)sigma_v_array[i_particle]*1e-3);
          fprintf(fp_out,"\n");
        }
        fclose(fp_out);
      }
      SID_Barrier(SID.COMM_WORLD);
    }
    j_species++;
    SID_log("Done.",SID_LOG_CLOSE);
    }  
  }    
  
  // Clean-up 
  free_plist(&plist);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
