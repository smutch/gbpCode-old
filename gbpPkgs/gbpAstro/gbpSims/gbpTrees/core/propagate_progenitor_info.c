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
  int i_read;
  int j_read;
  int l_read;
  int j_file;
  int i_write;
  int j_write;
  int l_write;
  int flag_init_write;
  i_write=i_write_last;
  j_write=j_write_last;
  l_write=l_write_last;
  for(i_read=i_write,j_read=j_write,l_read=l_write,j_file=0,flag_init_write=TRUE;
      j_write<=i_read_stop;
      i_read++,j_read+=i_read_step,l_read--,j_file++){

     // Open temporary horizontal tree file and read its header
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

        // Set the identy of the dominant substructures of each group.
        //    This needs to be complete for all descendants before propagate_n_particles_peak() 
        //    is called, so we have to do it here.
        propagate_dominant_substructures(groups_read,   n_groups,
                                         subgroups_read,n_subgroups,
                                         n_subgroups_group,
                                         i_read, // tree index
                                         j_read, // actual snapshot number
                                         l_read,
                                         i_read_step,
                                         n_wrap);
        if(n_groups_in!=n_groups[l_read])
           SID_trap_error("n_groups!=what it should be (ie.%d!=%d).  This might be a problem with the snapshot indexing.",ERROR_LOGIC,
                          n_groups_in,n_groups[l_read]);
        if(n_subgroups_in!=n_subgroups[l_read])
           SID_trap_error("n_subgroups!=what it should be (ie.%d!=%d).  This might be a problem with the snapshot indexing.",ERROR_LOGIC,
                          n_subgroups_in,n_subgroups[l_read]);
     }

     // Write updated tree files
     if(j_file>n_search){
        // If bridge fixing is on...
        if(flag_compute_fragmented){
           // ... propagate bridged halo information 
           propagate_bridge_info(groups_read,   n_groups,
                                 subgroups_read,n_subgroups,
                                 n_subgroups_group,
                                 i_write, // tree index
                                 j_write, // actual snapshot number
                                 l_write,
                                 i_read_step,
                                 n_wrap);
           // ... propagate fragmented halo flags
           propagate_fragmented_info(groups_read,   n_groups,
                                     subgroups_read,n_subgroups,
                                     n_subgroups_group,
                                     i_write, // tree index
                                     j_write, // actual snapshot number
                                     l_write,
                                     i_read_step,
                                     n_wrap);
        }
        // Propigate the peak halo size info. This 
        //    needs to be done for all halos before
        //    the merger information can be propagated.
        propagate_n_particles_peak(groups_read,   n_groups,
                                   subgroups_read,n_subgroups,
                                   n_subgroups_group,
                                   i_write, // tree index
                                   j_write, // actual snapshot number
                                   l_write,
                                   i_read_step,
                                   n_wrap);
        // Decide which halos are primaries and secondaries in merger events.
        propagate_merger_info(groups_read,   n_groups,
                              subgroups_read,n_subgroups,
                              n_subgroups_group,
                              i_write, // tree index
                              j_write, // actual snapshot number
                              l_write,
                              i_read_step,
                              n_wrap);
        // Write the results
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

