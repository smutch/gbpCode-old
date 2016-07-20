#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpSPH.h>

void allocate_gadget_particles(plist_info *plist, 
                               size_t     *n_particles_type_local,
                               size_t     *n_particles_type,
                               int         flag_long_IDs){
   SID_log("Allocating for particles...",SID_LOG_OPEN);
   for(int i_type=0;i_type<N_GADGET_TYPE;i_type++){
      if(n_particles_type[i_type]>0){
         GBPREAL *x_array =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_particles_type_local[i_type]);
         GBPREAL *y_array =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_particles_type_local[i_type]);
         GBPREAL *z_array =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_particles_type_local[i_type]);
         GBPREAL *vx_array=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_particles_type_local[i_type]);
         GBPREAL *vy_array=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_particles_type_local[i_type]);
         GBPREAL *vz_array=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_particles_type_local[i_type]);
         size_t  *id_array=(size_t  *)SID_malloc(sizeof(size_t) *n_particles_type_local[i_type]);
         ADaPS_store(&(plist->data),(void *)x_array, "x_%s", ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)y_array, "y_%s", ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)z_array, "z_%s", ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)vx_array,"vx_%s",ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)vy_array,"vy_%s",ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)vz_array,"vz_%s",ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)id_array,"id_%s",ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)(&(n_particles_type_local[i_type])),"n_%s",    ADaPS_SCALAR_SIZE_T,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)(&(n_particles_type[i_type])),      "n_all_%s",ADaPS_SCALAR_SIZE_T,plist->species[i_type]);
      }
   }
   ADaPS_store(&(plist->data),(void *)(&flag_long_IDs),"flag_long_IDs",ADaPS_SCALAR_INT);
   SID_log("Done.",SID_LOG_CLOSE);
}

