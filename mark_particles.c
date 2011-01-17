#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>

size_t mark_particles(plist_info *plist,
                      int         run_mode,
                      double     *input_vals,
                      char       *mark_name){
  size_t  n_particles;
  size_t  n_particles_local;
  size_t  i_particle;
  int     i_species;
  size_t *id;
  REAL   *x;
  REAL   *y;
  REAL   *z;
  REAL    r;
  int    *mark_array;
  int     flag_volume;
  int     flag_volume_sphere;
  size_t  n_marked_local=0;
  size_t  n_marked=0;
  
  // Interpret run-mode
  flag_volume=
    check_mode_for_flag(run_mode,VOLUME_BOX)||
    check_mode_for_flag(run_mode,VOLUME_SPHERE);
  flag_volume_sphere=check_mode_for_flag(run_mode,VOLUME_SPHERE);

  // Loop over all species
  for(i_species=0;i_species<N_GADGET_TYPE;i_species++){
    if(ADaPS_exist(plist->data,"n_all_%s",plist->species[i_species])){
      n_particles      =((size_t *)ADaPS_fetch(plist->data,"n_all_%s",plist->species[i_species]))[0];
      n_particles_local=((size_t *)ADaPS_fetch(plist->data,"n_%s",plist->species[i_species]))[0];
      // If this species has local particles
      if(n_particles_local>0){
        mark_array=(int *)SID_malloc(sizeof(int)*n_particles_local);
        // Mark particles in a volume
        if(flag_volume){
          x=(REAL *)ADaPS_fetch(plist->data,"x_%s",plist->species[i_species]);
          y=(REAL *)ADaPS_fetch(plist->data,"y_%s",plist->species[i_species]);
          z=(REAL *)ADaPS_fetch(plist->data,"z_%s",plist->species[i_species]);
          // Loop over all particles
          for(i_particle=0;i_particle<n_particles_local;i_particle++){
            mark_array[i_particle]=FALSE;
            switch(flag_volume_sphere){
            case TRUE:
              if(add_quad(3,
                          (double)(x[i_particle])-input_vals[0],
                          (double)(y[i_particle])-input_vals[1],
                          (double)(z[i_particle])-input_vals[2])<=input_vals[3])
                mark_array[i_particle]=TRUE;
              break;
            case FALSE:
              if(x[i_particle]>=(REAL)input_vals[0] && x[i_particle]<=(REAL)input_vals[1]){
                if(y[i_particle]>=(REAL)input_vals[2] && y[i_particle]<=(REAL)input_vals[3]){
                  if(z[i_particle]>=(REAL)input_vals[4] && z[i_particle]<=(REAL)input_vals[5]){
                    mark_array[i_particle]=TRUE;
                  }
                }
              }
              break;
            }
          }
        }
        // Mark particles by property
        else{
        }
        for(i_particle=0;i_particle<n_particles_local;i_particle++)
          if(mark_array[i_particle]) n_marked_local++;
        ADaPS_store(&(plist->data),(void *)mark_array,"%s_%s",ADaPS_DEFAULT,plist->species[i_species]);
      }
    }
  }
#ifdef USE_MPI
  MPI_Allreduce(&n_marked_local,&n_marked,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
#else
  n_marked=n_marked_local;
#endif
  return(n_marked);
}
