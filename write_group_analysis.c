#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHalos.h>

void write_group_analysis(FILE                 *fp_properties,
                          FILE                 *fp_profiles,
                          halo_properties_info *properties,
                          halo_profile_info    *profile){
  int i_bin;

  // Write properties
  fwrite(properties,sizeof(halo_properties_info),1,fp_properties);

  // Write profiles
  //fwrite(&(profile->n_bins),sizeof(int),                  1,              fp_profiles);
  //fwrite(profile->bins,     sizeof(halo_profile_bin_info),profile->n_bins,fp_profiles);

}

