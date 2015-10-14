#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHalos.h>

void write_group_analysis(FILE                 *fp_properties,
                          FILE                 *fp_profiles,
                          FILE                 *fp_indices,
                          halo_properties_info *properties,
                          halo_profile_info    *profile,
                          size_t               *R_index,
                          int                   n_particles){ // Needed because properties may not be defined
  int i_bin;

  // Write properties
  if(fp_properties!=NULL)
     fwrite(properties,sizeof(halo_properties_info),1,fp_properties);

  // Write profiles
  if(fp_profiles!=NULL){
     fwrite(&(profile->n_bins),sizeof(int),                  1,              fp_profiles);
     fwrite(profile->bins,     sizeof(halo_profile_bin_info),profile->n_bins,fp_profiles);
  }

  // Write the sort indices
  if(fp_indices!=NULL){
     fwrite(&n_particles,sizeof(int),   1,          fp_indices);
     fwrite(R_index,     sizeof(size_t),n_particles,fp_indices);
  }

}

