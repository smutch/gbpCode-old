#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees.h>

void compute_trees_horizontal_ghosts(int         *n_groups,
                                     int         *n_subgroups,
                                     int        **n_subgroups_group,
                                     int          i_read_start,
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
                                     char        *filename_cat_root_in,
                                     char        *filename_output_dir){

  // Build ghost-halo version of horizontal tree files
  SID_log("Creating ghost-populated trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Allocate storage for the ghost halos
  SID_log("Allocating arrays...",SID_LOG_OPEN|SID_LOG_TIMER);
  int  i_search;
  int  i_group;
  int  i_subgroup;
  int  j_subgroup;
  int *n_groups_ghost_alloc;
  int *n_subgroups_ghost_alloc;
  tree_horizontal_ghost_group_info    **groups_ghost;
  tree_horizontal_ghost_subgroup_info **subgroups_ghost;
  groups_ghost   =(tree_horizontal_ghost_group_info    **)SID_malloc(sizeof(tree_horizontal_ghost_group_info    *)*n_wrap);
  subgroups_ghost=(tree_horizontal_ghost_subgroup_info **)SID_malloc(sizeof(tree_horizontal_ghost_subgroup_info *)*n_wrap);
  for(i_search=0;i_search<n_wrap;i_search++){
     groups_ghost[i_search]   =(tree_horizontal_ghost_group_info    *)SID_calloc(sizeof(tree_horizontal_ghost_group_info)   *n_groups_max);
     subgroups_ghost[i_search]=(tree_horizontal_ghost_subgroup_info *)SID_calloc(sizeof(tree_horizontal_ghost_subgroup_info)*n_subgroups_max);
  }
  n_groups_ghost_alloc   =(int *)SID_calloc(sizeof(int)*n_files);
  n_subgroups_ghost_alloc=(int *)SID_calloc(sizeof(int)*n_files);
  SID_log("Done.",SID_LOG_CLOSE);

  // Determine a conservative estimate for the array allocation sizes needed for the 
  //   ghost-populated trees.  This may not be the final number but must at least be
  //   larger than the final number.
  SID_log("Counting ghosts...",SID_LOG_OPEN|SID_LOG_TIMER);
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
  int  i_read;
  int  j_read;
  int  l_read;
  int  i_file;
  int  j_file;
  for(i_read   =i_read_stop,
        i_file =i_file_start,
        j_file =0;
      i_read>=i_read_start;
      i_read-=i_read_step,
        i_file--,
        j_file++){

     // Read horizontal tree files (make sure we stay ahead by n_search snapshots)
     int n_groups_in;
     int n_subgroups_in;
     int n_trees_subgroup_in =0;
     int n_trees_group_in    =0;
     SID_log("Reading tree file for snapshot #%03d...",SID_LOG_OPEN,i_read);
     read_trees_horizontal((void **)groups_ghost,   &n_groups_in,
                           (void **)subgroups_ghost,&n_subgroups_in,
                           n_subgroups_group[i_file%n_wrap],
                           &n_trees_subgroup_in,
                           &n_trees_group_in,
                           i_file,
                           i_read,
                           j_file,
                           n_wrap,
                           filename_output_dir,
                           TREE_HORIZONTAL_STORE_GHOSTS);
     if(n_groups_in!=n_groups[j_file])
        SID_trap_error("n_groups!=what it should be (ie.%d!=%d).  This might be a problem with the snapshot indexing.",ERROR_LOGIC,
                       n_groups_in,n_groups[j_file]);
     if(n_subgroups_in!=n_subgroups[j_file])
        SID_trap_error("n_subgroups!=what it should be (ie.%d!=%d).  This might be a problem with the snapshot indexing.",ERROR_LOGIC,
                       n_subgroups_in,n_subgroups[j_file]);
     SID_log("Done.",SID_LOG_CLOSE);

     // Create ghosts and write updated tree files and ghost catalogs
     SID_log("Creating ghost halos for snapshot #%03d...",SID_LOG_OPEN,i_read);
     process_ghosts(groups_ghost,
                    subgroups_ghost,
                    n_groups,
                    n_subgroups,
                    n_subgroups_group,
                    n_groups_ghost_alloc,
                    n_subgroups_ghost_alloc,
                    NULL,
                    NULL,
                    i_file,
                    i_read,
                    j_file,
                    i_file_start,
                    n_search,
                    n_wrap,
                    a_list,
                    cosmo);
     SID_log("Done.",SID_LOG_CLOSE);
  }

  // Create an array holding the total number of halos and
  //   find the largest element for memory allocation purposes.
  int n_groups_ghost_alloc_max;
  int n_subgroups_ghost_alloc_max;
  for(i_file=0;i_file<n_files;i_file++){
     n_groups_ghost_alloc[i_file]   +=n_groups[n_files-i_file-1];
     n_subgroups_ghost_alloc[i_file]+=n_subgroups[n_files-i_file-1];
  }
  calc_max(n_groups_ghost_alloc,   &n_groups_ghost_alloc_max,   n_files,SID_INT,CALC_MODE_DEFAULT);
  calc_max(n_subgroups_ghost_alloc,&n_subgroups_ghost_alloc_max,n_files,SID_INT,CALC_MODE_DEFAULT);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  SID_log("Done.",SID_LOG_CLOSE);

  // Allocate storage for the ghost halos
  SID_log("Reallocating/allocating arrays...",SID_LOG_OPEN|SID_LOG_TIMER);
  halo_properties_info ***group_properties;
  halo_properties_info ***subgroup_properties;
  int                    *n_groups_ghost_used;
  int                    *n_subgroups_ghost_used;
  for(i_search=0;i_search<n_wrap;i_search++){
     SID_free(SID_FARG groups_ghost[i_search]);
     SID_free(SID_FARG subgroups_ghost[i_search]);
  }
  SID_free(SID_FARG groups_ghost);
  SID_free(SID_FARG subgroups_ghost);
  SID_free(SID_FARG n_groups_ghost_alloc);
  SID_free(SID_FARG n_subgroups_ghost_alloc);
  groups_ghost       =(tree_horizontal_ghost_group_info    **)SID_malloc(sizeof(tree_horizontal_ghost_group_info     *)*n_wrap);
  subgroups_ghost    =(tree_horizontal_ghost_subgroup_info **)SID_malloc(sizeof(tree_horizontal_ghost_subgroup_info  *)*n_wrap);
  group_properties   =(halo_properties_info               ***)SID_malloc(sizeof(halo_properties_info                **)*n_wrap);
  subgroup_properties=(halo_properties_info               ***)SID_malloc(sizeof(halo_properties_info                **)*n_wrap);
  for(i_search=0;i_search<n_wrap;i_search++){
     groups_ghost[i_search]       =(tree_horizontal_ghost_group_info    *)SID_calloc(sizeof(tree_horizontal_ghost_group_info)   *n_groups_ghost_alloc_max);
     subgroups_ghost[i_search]    =(tree_horizontal_ghost_subgroup_info *)SID_calloc(sizeof(tree_horizontal_ghost_subgroup_info)*n_subgroups_ghost_alloc_max);
     group_properties[i_search]   =(halo_properties_info               **)SID_malloc(sizeof(halo_properties_info *)             *n_groups_max);
     subgroup_properties[i_search]=(halo_properties_info               **)SID_malloc(sizeof(halo_properties_info *)             *n_subgroups_max);
     for(i_group=0;   i_group<n_groups_max;      i_group++)    group_properties[i_search][i_group]      =NULL;
     for(i_subgroup=0;i_subgroup<n_subgroups_max;i_subgroup++) subgroup_properties[i_search][i_subgroup]=NULL;
  }
  n_groups_ghost_used   =(int *)SID_calloc(sizeof(int)*n_files);
  n_subgroups_ghost_used=(int *)SID_calloc(sizeof(int)*n_files);
  SID_log("Done.",SID_LOG_CLOSE);

  // Loop over each snapshot
  SID_log("Generating ghost-populated trees...",SID_LOG_OPEN|SID_LOG_TIMER);
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
  int i_write;
  int j_write;
  int l_write;
  for(i_read   =i_read_stop,
        i_file =i_file_start,
        j_file =0,
        i_write=i_file_start,
        j_write=i_read_stop,
        l_write=0;
      j_write>=i_read_start;
        i_write--,
        j_write-=i_read_step,
        l_write++){

     // Read horizontal tree files (make sure we stay ahead by n_search snapshots)
     int n_groups_in;
     int n_subgroups_in;
     int n_trees_subgroup_in =0;
     int n_trees_group_in    =0;
     while(i_read>=i_read_start && i_file>=(i_write-n_search)){
        SID_log("Reading tree file for snapshot #%03d...",SID_LOG_OPEN,i_read);
        read_trees_horizontal((void **)groups_ghost,   &n_groups_in,
                              (void **)subgroups_ghost,&n_subgroups_in,
                              n_subgroups_group[i_file%n_wrap],
                              &n_trees_subgroup_in,
                              &n_trees_group_in,
                              i_file,
                              i_read,
                              j_file,
                              n_wrap,
                              filename_output_dir,
                              TREE_HORIZONTAL_STORE_GHOSTS);
        if(n_groups_in!=n_groups[j_file])
           SID_trap_error("n_groups!=what it should be (ie.%d!=%d).  This might be a problem with the snapshot indexing.",ERROR_LOGIC,
                          n_groups_in,n_groups[j_file]);
        if(n_subgroups_in!=n_subgroups[j_file])
           SID_trap_error("n_subgroups!=what it should be (ie.%d!=%d).  This might be a problem with the snapshot indexing.",ERROR_LOGIC,
                          n_subgroups_in,n_subgroups[j_file]);
        SID_log("Done.",SID_LOG_CLOSE);

        SID_log("Reading catalog properties for snapshot #%03d...",SID_LOG_OPEN,i_read);
        read_catalog_ghost_interpolation(groups_ghost[i_file%n_wrap],
                                         group_properties[i_file%n_wrap],
                                         n_groups[j_file],
                                         subgroups_ghost[i_file%n_wrap],
                                         subgroup_properties[i_file%n_wrap],
                                         n_subgroups[j_file],
                                         filename_cat_root_in,
                                         i_file,
                                         i_read,
                                         n_wrap);
        SID_log("Done.",SID_LOG_CLOSE);

        // Create ghosts and write updated tree files and ghost catalogs
        SID_log("Creating ghost halos for snapshot #%03d...",SID_LOG_OPEN,i_read);
        process_ghosts(groups_ghost,
                       subgroups_ghost,
                       n_groups,
                       n_subgroups,
                       n_subgroups_group,
                       NULL,
                       NULL,
                       n_groups_ghost_used,
                       n_subgroups_ghost_used,
                       i_file,
                       i_read,
                       j_file,
                       i_file_start,
                       n_search,
                       n_wrap,
                       a_list,
                       cosmo);
        SID_log("Done.",SID_LOG_CLOSE);

        i_read-=i_read_step;
        i_file--;
        j_file++;
     }

     // Set and check validity of total halo counts
     int n_groups_ghost;
     int n_subgroups_ghost;
     n_groups_ghost   =n_groups[n_files-i_write-1]   +n_groups_ghost_used[i_write];
     n_subgroups_ghost=n_subgroups[n_files-i_write-1]+n_subgroups_ghost_used[i_write];
     if(n_groups_ghost>n_groups_ghost_alloc_max)
        SID_trap_error("The group count exceeds the group allocation size (ie. %d>%d).",ERROR_LOGIC,n_groups_ghost,n_groups_ghost_alloc_max);
     if(n_subgroups_ghost>n_subgroups_ghost_alloc_max)
        SID_trap_error("The subgroup count exceeds the subgroup allocation size (ie. %d>%d).",ERROR_LOGIC,n_subgroups_ghost,n_subgroups_ghost_alloc_max);

     // Write the ghost-populated trees
     write_trees_horizontal((void **)groups_ghost,
                            (void **)subgroups_ghost,
                            n_groups_ghost,   n_groups_ghost_alloc_max,   
                            n_subgroups_ghost,n_subgroups_ghost_alloc_max,
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
                            filename_cat_root_in,
                            filename_output_dir,
                            a_list,
                            cosmo,
                            n_k_match,
                            TREE_HORIZONTAL_WRITE_GHOSTS|TREE_HORIZONTAL_WRITE_NOCASES);

     // Write the ghost catalog files
     write_ghost_catalog(groups_ghost[i_write%n_wrap],
                         group_properties,
                         n_groups_ghost-n_groups[n_files-i_write-1],
                         n_groups[n_files-i_write-1],
                         subgroups_ghost[i_write%n_wrap],
                         subgroup_properties,
                         n_subgroups_ghost-n_subgroups[n_files-i_write-1],
                         n_subgroups[n_files-i_write-1],
                         filename_output_dir,
                         filename_cat_root_in,
                         i_write,
                         j_write,
                         l_write,
                         n_wrap,
                         a_list,
                         cosmo);
  }
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  SID_log("Freeing arrays...",SID_LOG_OPEN);
  for(i_search=0;i_search<n_wrap;i_search++){
     SID_free(SID_FARG groups_ghost[i_search]);
     SID_free(SID_FARG subgroups_ghost[i_search]);
     for(i_group=0;   i_group<n_groups_max;      i_group++)    SID_free(SID_FARG group_properties[i_search][i_group]);
     for(i_subgroup=0;i_subgroup<n_subgroups_max;i_subgroup++) SID_free(SID_FARG subgroup_properties[i_search][i_subgroup]);
     SID_free(SID_FARG group_properties[i_search]);
     SID_free(SID_FARG subgroup_properties[i_search]);
  }
  SID_free(SID_FARG groups_ghost);
  SID_free(SID_FARG subgroups_ghost);
  SID_free(SID_FARG n_groups_ghost_used);
  SID_free(SID_FARG n_subgroups_ghost_used);
  SID_free(SID_FARG group_properties);
  SID_free(SID_FARG subgroup_properties);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
}

