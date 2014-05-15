#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpTrees_build.h>

void read_halo_sizes(tree_horizontal_info **halos,
                     int     n_halos_i,
                     int    *match_id,
                     float  *match_score,
                     size_t *match_index,
                     int    *n_particles,
                     int     i_file,
                     int     i_read,
                     int     i_read_step,
                     int     n_wrap,
                     int     n_halos_max,
                     char   *filename_root_matches,
                     int     flag_match_subgroups){

    SID_log("Reading halo sizes...",SID_LOG_OPEN|SID_LOG_TIMER);
    SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);

    // Read back-matching (we assume here that this code
    //    is not called for the highest-numbered snapshot)
    int n_halos_1_matches;
    int n_halos_2_matches;
    read_matches(filename_root_matches,
                 i_read,i_read+i_read_step,n_halos_max,
                 flag_match_subgroups,
                 &n_halos_1_matches,
                 &n_halos_2_matches,
                 NULL,
                 n_particles,
                 NULL,
                 NULL,
                 match_id,
                 match_score,
                 match_index);

    // Sanity check
    if(n_halos_i!=n_halos_1_matches)
       SID_trap_error("Halo count mismatch in read_halo_sizes().",ERROR_LOGIC);

    // Store halo sizes for the current snapshot's halos
    for(int i_halo=0;i_halo<n_halos_1_matches;i_halo++)
       halos[i_file%n_wrap][i_halo].n_particles=n_particles[i_halo];

    SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
    SID_log("Done.",SID_LOG_CLOSE);
}

