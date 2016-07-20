#include <gbpLib.h>
#include <gbpSPH.h>

void free_plist(plist_info *plist){
  int status;
  ADaPS_free(SID_FARG plist->data);
  free_types(&(plist->species),plist->n_species);
  plist->length_unit=0.0;
  plist->mass_unit  =0.0;
  plist->time_unit  =0.0;
  plist->d_x    =0.0;
  plist->d_y    =0.0;
  plist->d_z    =0.0;
  plist->d_vx   =0.0;
  plist->d_vy   =0.0;
  plist->d_vz   =0.0;
  plist->d_alpha=0.0;
  plist->d_beta =0.0;
  plist->d_gamma=0.0;
}
