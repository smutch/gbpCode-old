#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){
  char        filename_SSimPL_dir[256];
  char        filename_SSimPL_base[256];
  char        filename_trees_name[256];
  char        filename_halo_version_root[256];
  char        filename_halo_root_in[256];
  char        filename_snap_list_in[256];
  char        filename_cat_root_in[256];
  char        filename_root_out[256];
  char        filename_output_dir_horizontal_trees[256];
  char        filename_output_dir[256];
  char        filename_output_dir_horizontal[256];
  char        filename_root_matches[256];
  char        filename_snap_list_out[256];
  char        filename_output_file_root[256];
  int         i_read_start;
  int         i_read_stop;
  int         i_read_step;
  int         n_search;
  int         n_dim_files;
  int         n_k_match=2;
  int         flag_clean=FALSE;
  FILE       *fp_in;
  FILE       *fp_out;
  char       *line=NULL;
  size_t      line_length=0;
  int         n_snaps,i_read,i_next,i_write,n_keep;
  double     *a_list_in;
  double     *a_list_out;
  cosmo_info *cosmo;
  int         flag_fix_bridges;

  SID_init(&argc,&argv,NULL,NULL);

  // Initialize cosmology
  init_cosmo_std(&cosmo);

  // Fetch user inputs
  double box_size;
  if(argc!=11)
    SID_trap_error("Incorrect syntax",ERROR_SYNTAX);
  strcpy(filename_SSimPL_dir,       argv[1]);
  strcpy(filename_halo_version_root,argv[2]);
  strcpy(filename_trees_name,       argv[3]);
  i_read_start                =atoi(argv[4]);
  i_read_stop                 =atoi(argv[5]);
  i_read_step                 =atoi(argv[6]);
  n_search                    =atoi(argv[7]);
  flag_fix_bridges            =atoi(argv[8]);
  box_size            =(double)atof(argv[9]);
  n_dim_files                 =atoi(argv[10]);

  SID_log("Constructing merger trees for snapshots #%d->#%d (step=%d, n_search=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,
          i_read_start,i_read_stop,i_read_step,n_search);
  sprintf(filename_SSimPL_base,"%s",filename_SSimPL_dir);
  strip_path(filename_SSimPL_base);
  sprintf(filename_halo_root_in,"%s/halos/%s",     filename_SSimPL_dir,filename_halo_version_root);
  sprintf(filename_cat_root_in, "%s/catalogs/%s",  filename_SSimPL_dir,filename_halo_version_root);
  sprintf(filename_root_matches,"%s/trees/matches",filename_SSimPL_dir);
  sprintf(filename_snap_list_in,"%s/run/%s.a_list",filename_SSimPL_dir,filename_SSimPL_base);
  sprintf(filename_root_out,    "%s/trees/%s",     filename_SSimPL_dir,filename_trees_name);
  compute_trees_horizontal(filename_halo_root_in,
                           filename_cat_root_in,
                           filename_snap_list_in,
                           filename_root_matches,
                           filename_root_out,
                           &cosmo,
                           i_read_start,
                           i_read_stop,
                           i_read_step,
                           n_search,
                           flag_fix_bridges,
                           &flag_clean);
  compute_trees_vertical(filename_SSimPL_dir,
                         filename_halo_version_root,
                         filename_trees_name,
                         box_size,
                         n_dim_files);
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  ADaPS_free(SID_FARG cosmo);  

  SID_exit(ERROR_NONE);
}

