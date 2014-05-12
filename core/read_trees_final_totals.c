#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void read_trees_final_totals(char *filename_output_dir_horizontal_trees,
                             int   i_read_start,
                             int   i_read_stop, 
                             int   i_read_step,
                             int  *i_read_last,
                             int  *n_snap,
                             int  *n_groups_max_in,
                             int  *n_subgroups_max_in,
                             int  *n_progenitors_max,
                             int  *n_trees_subgroup,
                             int  *n_trees_group){

  SID_log("Getting final halo/tree/snapshot totals...",SID_LOG_OPEN);

  // Count the number of snaps used by these trees
  int i_read;
  for(i_read=i_read_stop,(*n_snap)=0;i_read>=i_read_start;i_read-=i_read_step) (*n_snap)++;

  // Figure-out what the last-used snapshot is so we can read its header
  (*i_read_last)=i_read_stop;
  for(i_read=i_read_stop;i_read>=i_read_start;i_read-=i_read_step) (*i_read_last)=i_read;

  // Determining iso-tree id mappings (onto tree structures, ranks and files).  Read header info
  //    from the i_read_start'th file since that is where the final tree tally gets written.
  char   filename_in[MAX_FILENAME_LENGTH];
  SID_fp fp_in;
  int    n_step_in;
  int    n_search_in;
  int    n_groups;
  int    n_subgroups;
  sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_output_dir_horizontal_trees,(*i_read_last));
  SID_fopen(filename_in,"r",&fp_in);
  SID_fread_all(&n_step_in,         sizeof(int),1,&fp_in);
  SID_fread_all(&n_search_in,       sizeof(int),1,&fp_in);
  SID_fread_all(&n_groups,          sizeof(int),1,&fp_in);
  SID_fread_all(&n_subgroups,       sizeof(int),1,&fp_in);
  SID_fread_all( n_groups_max_in,   sizeof(int),1,&fp_in);
  SID_fread_all( n_subgroups_max_in,sizeof(int),1,&fp_in);
  SID_fread_all( n_trees_subgroup,  sizeof(int),1,&fp_in);
  SID_fread_all( n_trees_group,     sizeof(int),1,&fp_in);
  SID_fclose(&fp_in);
  (*n_progenitors_max)=MAX((*n_groups_max_in),(*n_subgroups_max_in));
  if(n_step_in!=i_read_step) SID_trap_error("Snapshot step sizes don't match (ie. %d!=%d)",ERROR_LOGIC,n_step_in,i_read_step);

  SID_log("Done.",SID_LOG_CLOSE);
}
