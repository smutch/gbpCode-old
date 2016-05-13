#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void propagate_progenitor_info(int         *n_groups,
                               int         *n_subgroups,
                               int        **n_subgroups_group,
                               int          i_file_start,
                               int          i_write_last,
                               int          j_write_last,
                               int          l_write_last,
                               int          i_read_stop,
                               int          i_read_step,
                               int          max_tree_id_subgroup,
                               int          max_tree_id_group,
                               int          n_subgroups_max,
                               int          n_groups_max,
                               int          n_search,
                               int          n_files,
                               int          n_wrap,
                               int          n_k_match,
                               double      *a_list,
                               cosmo_info **cosmo,
                               char        *filename_output_dir,
                               int          flag_compute_fragmented){

  SID_log("Propagating progenitor information...",SID_LOG_OPEN|SID_LOG_TIMER);
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);

  // Sanity check
  if(n_wrap<=(2*n_search+1))
     SID_trap_error("n_skip is not large enough in propagate_progenitor_info() (ie. %d<=%d)",ERROR_LOGIC,n_wrap,(2*n_search+1));

  // Allocate new data structures for propagating information
  tree_horizontal_extended_info **subgroups_read;
  tree_horizontal_extended_info **groups_read;
  subgroups_read=(tree_horizontal_extended_info **)SID_malloc(sizeof(tree_horizontal_extended_info *)*n_wrap);
  groups_read   =(tree_horizontal_extended_info **)SID_malloc(sizeof(tree_horizontal_extended_info *)*n_wrap);
  for(int i_search=0;i_search<n_wrap;i_search++){
     subgroups_read[i_search]=(tree_horizontal_extended_info *)SID_calloc(sizeof(tree_horizontal_extended_info)*n_subgroups_max);
     groups_read[i_search]   =(tree_horizontal_extended_info *)SID_calloc(sizeof(tree_horizontal_extended_info)*n_groups_max);
  }

  // Loop over all the used snapshots (first set write counters to last-used values)
  int flag_init_write=TRUE;
  int i_read   =i_write_last;
  int j_read   =j_write_last;
  int l_read   =l_write_last;
  int i_process=i_read;
  int j_process=j_read;
  int l_process=l_read;
  int i_write  =i_read;
  int j_write  =j_read;
  int l_write  =l_read;
  int j_file   =0;
  for(;
      j_write<=i_read_stop;
      i_read++,j_read+=i_read_step,l_read--,j_file++){

     // Read temporary tree files
     if(j_read<=i_read_stop){
        int n_groups_in;
        int n_subgroups_in;
        int n_trees_subgroup_in =0;
        int n_trees_group_in    =0;
        read_trees_horizontal((void **)groups_read,   &n_groups_in,
                              (void **)subgroups_read,&n_subgroups_in,
                              n_subgroups_group[i_read%n_wrap],
                              &n_trees_subgroup_in,
                              &n_trees_group_in,
                              i_read,
                              j_read,
                              l_read,
                              n_wrap,
                              filename_output_dir,
                              TREE_HORIZONTAL_READ_EXTENDED|TREE_HORIZONTAL_STORE_EXTENDED);

        // *** Progenitor propagation that only references progenitors and not descendants can be done here ***

        // Sanity checks
        if(n_groups_in!=n_groups[l_read])
           SID_trap_error("n_groups!=what it should be (ie.%d!=%d).  This might be a problem with the snapshot indexing.",ERROR_LOGIC,
                          n_groups_in,n_groups[l_read]);
        if(n_subgroups_in!=n_subgroups[l_read])
           SID_trap_error("n_subgroups!=what it should be (ie.%d!=%d).  This might be a problem with the snapshot indexing.",ERROR_LOGIC,
                          n_subgroups_in,n_subgroups[l_read]);
     }

     // Wait until we've loaded n_search+1 snapshots
     //    before starting to process or write results.
     if(j_file>n_search){
        if(j_process<=i_read_stop){
           // Propagate fragmented halo flags
           if(flag_compute_fragmented)
              propagate_fragmented_info(groups_read,   n_groups,
                                        subgroups_read,n_subgroups,
                                        n_subgroups_group,
                                        i_process, // tree index
                                        j_process, // actual snapshot number
                                        l_process,
                                        i_read_step,
                                        n_wrap);
           // Set the identy of the dominant substructures of each group.
           //    This needs to be complete for all descendants before propagate_n_particles_peak() 
           //    is called, so we have to do it here.
           propagate_dominant_substructures(groups_read,   n_groups,
                                            subgroups_read,n_subgroups,
                                            n_subgroups_group,
                                            i_process, // tree index
                                            j_process, // actual snapshot number
                                            l_process,
                                            i_read_step,
                                            n_wrap);
           // Propagate the peak halo size info. This 
           //    needs to be done for all halos before
           //    the merger information can be propagated.
           //    Fragmented flags for the descendant snapshot
           //    needs to have been set, but not for 
           //    subsequent ones.
           propagate_n_particles_peak(groups_read,   n_groups,
                                      subgroups_read,n_subgroups,
                                      n_subgroups_group,
                                      i_process, // tree index
                                      j_process, // actual snapshot number
                                      l_process,
                                      i_read_step,
                                      n_wrap);
           // Decide which halos are primaries and secondaries in merger events.
           propagate_merger_info(groups_read,   n_groups,
                                 subgroups_read,n_subgroups,
                                 n_subgroups_group,
                                 i_process, // tree index
                                 j_process, // actual snapshot number
                                 l_process,
                                 i_read_step,
                                 n_wrap);
           // Propagate tree IDs along the progenitor lines with the highest peak particle count.
           //    This needs to be done after mergers have been processed because we use
           //    the tree case merger flags to do this.
           propagate_tree_ids(groups_read,   n_groups,
                              subgroups_read,n_subgroups,
                              n_subgroups_group,
                              i_process, // tree index
                              j_process, // actual snapshot number
                              l_process,
                              i_read_step,
                              n_wrap);
           // Propagate bridged halo information...
           //    (n.b.: mergers need to be done before this because we will use those flags here)
           if(flag_compute_fragmented)
              // ... propagate bridged halo information 
              propagate_bridge_info(groups_read,   n_groups,
                                    subgroups_read,n_subgroups,
                                    n_subgroups_group,
                                    i_process, // tree index
                                    j_process, // actual snapshot number
                                    l_process,
                                    i_read_step,
                                    n_wrap);
        }
        i_process++;
        l_process--;
        j_process+=i_read_step;

        // Write the results.  Wait until another n_search snapshots
        //    have been read/processed before performing the write, to make sure
        //    that any changes that need to be made -- not just to descendants --
        //    but also to progenitors of descendants (such as what happens in
        //    propagate_merger_info(), where progenitors are changed but
        //    nothing can be done until peak halo sizes are set) are complete.
        if(j_file>(2*n_search)){
           write_trees_horizontal((void **)groups_read,
                                  (void **)subgroups_read,
                                  n_groups[l_write],   n_groups_max,   
                                  n_subgroups[l_write],n_subgroups_max,
                                  n_subgroups_group,
                                  max_tree_id_subgroup,
                                  max_tree_id_group,
                                  i_write,
                                  j_write,
                                  l_write,
                                  i_read_step,
                                  n_search,
                                  n_wrap,
                                  i_file_start,
                                  NULL,
                                  filename_output_dir,
                                  a_list,
                                  cosmo,
                                  n_k_match,
                                  flag_init_write,
                                  TREE_HORIZONTAL_WRITE_DEFAULT);
           i_write++;
           l_write--;
           j_write+=i_read_step;
           flag_init_write=FALSE;
        }
     }
  }

  // Clean-up
  for(int i_search=0;i_search<n_wrap;i_search++){
     SID_free(SID_FARG subgroups_read[i_search]);
     SID_free(SID_FARG groups_read[i_search]);
  }
  SID_free(SID_FARG subgroups_read);
  SID_free(SID_FARG groups_read);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  SID_log("Done.",SID_LOG_CLOSE);
}

