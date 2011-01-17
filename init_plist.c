#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>

void init_plist(plist_info *plist,
                slab_info  *slab,
                double      length_unit,
                double      mass_unit,
                double      velocity_unit){

  // Data library
  ADaPS_init(&(plist->data));

  // Species types and names
  plist->n_species=N_GADGET_TYPE;
  prep_types(&(plist->species),
             plist->n_species,
             GADGET_TYPE_GAS,  "gas",
             GADGET_TYPE_DARK, "dark",
             GADGET_TYPE_STAR, "star",
             GADGET_TYPE_STAR2,"star2",
             GADGET_TYPE_STAR3,"star3",
             GADGET_TYPE_STAR4,"star4");

  // Units
  plist->length_unit=length_unit;
  plist->mass_unit  =mass_unit;
  if(velocity_unit>0.){
    plist->velocity_unit=velocity_unit;
    plist->time_unit    =plist->length_unit/plist->velocity_unit;
  }
  else{
    plist->time_unit    =sqrt(pow(plist->length_unit,3.0)/(G_NEWTON*plist->mass_unit));
    plist->velocity_unit=plist->length_unit/plist->time_unit;
  }

  // Offsets and rotations
  plist->d_x    =0.0;
  plist->d_y    =0.0;
  plist->d_z    =0.0;
  plist->d_vx   =0.0;
  plist->d_vy   =0.0;
  plist->d_vz   =0.0;
  plist->d_alpha=0.0;
  plist->d_beta =0.0;
  plist->d_gamma=0.0;

  // Zone decomposition (only supports x-slabs right now)
  if(slab!=NULL){
    ADaPS_store(&(plist->data),(void *)(&(slab->x_min_local)),"rank_x_min",ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&(plist->data),(void *)(&(slab->x_max_local)),"rank_x_max",ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&(plist->data),(void *)(&(slab->x_max)),      "x_period",  ADaPS_SCALAR_DOUBLE);
  }

}
