#define  _MAIN
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>

int main(int argc, char *argv[]){
  double  a_start;
  int     n_species;
  int     i,j,i_x,i_y,i_z;
  int     i_jack;
  int     n_jack_total;
  char    species_name[256];
  double  h_Hubble;
  double  Omega_Lambda;
  double  Omega_M;
  double  Omega_b;
  double  f_gas;
  double  Omega_k;
  double  sigma_8;
  double  n_spec;
  double  redshift;
  int     i_rank;
  int     i_group;
  int     n_groups_all;
  int     n_subgroups_all;
  int     n_groups;
  size_t  n_group;
  size_t  n_group_local;
  char    filename_groups_root[MAX_FILENAME_LENGTH];
  char    filename_groups_in[MAX_FILENAME_LENGTH];
  char    filename_subgroups_in[MAX_FILENAME_LENGTH];
  char    filename_cat_root[MAX_FILENAME_LENGTH];
  char    filename_out[MAX_FILENAME_LENGTH];
  char    filename_out_root[MAX_FILENAME_LENGTH];
  char    group_text_prefix[4];
  int     i_snap;
  int     n_groupings;
  int     n_subgroups;
  int     cfunc_mode;
  int     i_compute;
  char    group_name[6];
  char    filename_TF[256];
  char    n_string[64];
  int             n[3];
  double          box_size;
  double          L[3];
  size_t          n_all;
  cosmo_info     *cosmo;
  field_info      FFT;
  plist_info      plist;
  FILE           *fp;
  FILE           *fp_groups;
  FILE           *fp_subgroups;
  int     flag_write_header;
  int     i_temp;
  int     n_temp;
  double *r_temp;
  double *kmin_temp;
  double *kmax_temp;
  int     n_jack;
  int     i_subhalo,i_grouping;
  halo_properties_info  group_properties;
  halo_properties_info  subgroup_properties;
  int   n_halos_per_grouping,n_halos,n_halos_all,i_halo,j_halo;
  int   n_particles_min;
  int  *n_subgroups_group;
  REAL *x_halos;
  REAL *y_halos;
  REAL *z_halos;
  REAL *r_halos;
  REAL *vx_halos_FoF;
  REAL *vy_halos_FoF;
  REAL *vz_halos_FoF;
  REAL *vx_halos_sub;
  REAL *vy_halos_sub;
  REAL *vz_halos_sub;
  double *M_halos;
  size_t *V_halos_index;
  REAL *x_group;
  REAL *y_group;
  REAL *z_group;
  REAL *vx_group;
  REAL *vy_group;
  REAL *vz_group;
  REAL *V_halos;
  REAL  x_min,x_max;
  REAL  y_min,y_max;
  REAL  z_min,z_max;
  int     n_particles;
  double  M_min,M_max,M_med;
  double  V_min,V_max,V_med;
  char   *line=NULL;
  int     line_length=0;
  int     n_spline,n_model;
  double *r_spline;
  double *P_spline;
  double *r_model;
  double *P_model;
  int     i_bin,j_bin;
  double r_o,r_X;
  double *x_interp;
  double *y_interp;
  int    i_file,n_files;
  int    mass_assignment_scheme;
  double r_min_l1D;
  double r_max_1D;
  double dr_1D;
  double dlr_l1D;
  double r_min_2D;
  double r_max_2D;
  double dr_2D;
  double *CFUNC_l1D=NULL;
  double *dCFUNC_l1D=NULL;
  double *COVMTX_l1D=NULL;
  double *CFUNC_1D=NULL;
  double *dCFUNC_1D=NULL;
  double *COVMTX_1D=NULL;
  double *CFUNC_2D=NULL;
  double *dCFUNC_2D=NULL;
  int     n_1D,n_2D;
  int     F_random;
  size_t  n_random,n_random_local;
  int     i_random,j_random;
  RNG_info  RNG;
  int       seed=1327621;
  REAL     *x_random;
  REAL     *y_random;
  REAL     *z_random;
  interp_info *r_o_interp;
  int          flag_log_1D;
  int          flag_compute_RR;
  long long **DD_l1D;
  long long **DR_l1D;
  long long **RR_l1D;
  long long **DD_1D;
  long long **DR_1D;
  long long **RR_1D;
  long long **DD_2D;
  long long **DR_2D;
  long long **RR_2D;
  long long   n_random_ll;
  long long   n_data_ll;
  size_t      n_rank;
  REAL        V_grouping,V_max_lo,V_max_hi,delta_V_max;

  // Initialization -- MPI etc.
  SID_init(&argc,&argv,NULL);
  strcpy(filename_groups_root,argv[1]);
  strcpy(filename_cat_root,   argv[2]);
  strcpy(filename_out_root,   argv[3]);
  i_snap              =(int) atoi(argv[4]);
  n_halos_per_grouping=(int) atoi(argv[5]);
  n_particles_min     =(int) atoi(argv[6]);
  V_max_lo            =(REAL)atof(argv[7]);
  V_max_hi            =(REAL)atof(argv[8]);
  n_groupings         =(int) atoi(argv[9]);
  delta_V_max         =(V_max_hi-V_max_lo)/(double)(n_groupings-1);

  SID_log("Producing %d groupings of halos...",SID_LOG_OPEN,n_groupings);

  init_cosmo_std(&cosmo);
  init_RNG(&seed,&RNG,RNG_DEFAULT);

  // Read group catalogs
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  read_groups(filename_groups_root,i_snap,READ_GROUPS_ALL|READ_GROUPS_NOIDS,&plist,"groups");
  n_subgroups      =((int *)ADaPS_fetch(plist.data,"n_subgroups_all_groups"))[0];
  n_groups         =((int *)ADaPS_fetch(plist.data,"n_groups_all_groups"))[0];
  n_subgroups_group= (int *)ADaPS_fetch(plist.data,"n_subgroups_group_groups");
  SID_log("n_groups   =%d",SID_LOG_COMMENT,n_groups);
  SID_log("n_subgroups=%d",SID_LOG_COMMENT,n_subgroups);

  // Read catalog ...
  SID_log("Reading catalog...",SID_LOG_OPEN);

  //   ... open files and skip headers ...
  sprintf(filename_groups_in,   "%s_%03d.catalog_groups_properties",   filename_cat_root,i_snap);
  sprintf(filename_subgroups_in,"%s_%03d.catalog_subgroups_properties",filename_cat_root,i_snap);
  fp_groups=fopen(filename_groups_in,"r");
  fread(&i_file,      sizeof(int),1,fp_groups);
  fread(&n_files,     sizeof(int),1,fp_groups);
  fread(&n_groups,    sizeof(int),1,fp_groups);
  fread(&n_groups_all,sizeof(int),1,fp_groups);
  fp_subgroups=fopen(filename_subgroups_in,"r");
  fread(&i_file,         sizeof(int),1,fp_subgroups);
  fread(&n_files,        sizeof(int),1,fp_subgroups);
  fread(&n_subgroups,    sizeof(int),1,fp_subgroups);
  fread(&n_subgroups_all,sizeof(int),1,fp_subgroups);
  SID_log("%d groups and %d subgroups...",SID_LOG_CONTINUE,n_groups,n_subgroups);

  //   ... allocate arrays ...
  n_group     =(size_t)n_halos;
  x_halos     =(REAL   *)SID_malloc(sizeof(REAL)*n_subgroups_all);
  y_halos     =(REAL   *)SID_malloc(sizeof(REAL)*n_subgroups_all);
  z_halos     =(REAL   *)SID_malloc(sizeof(REAL)*n_subgroups_all);
  r_halos     =(REAL   *)SID_malloc(sizeof(REAL)*n_subgroups_all);
  vx_halos_FoF=(REAL   *)SID_malloc(sizeof(REAL)*n_subgroups_all);
  vy_halos_FoF=(REAL   *)SID_malloc(sizeof(REAL)*n_subgroups_all);
  vz_halos_FoF=(REAL   *)SID_malloc(sizeof(REAL)*n_subgroups_all);
  vx_halos_sub=(REAL   *)SID_malloc(sizeof(REAL)*n_subgroups_all);
  vy_halos_sub=(REAL   *)SID_malloc(sizeof(REAL)*n_subgroups_all);
  vz_halos_sub=(REAL   *)SID_malloc(sizeof(REAL)*n_subgroups_all);
  M_halos     =(double *)SID_malloc(sizeof(double)*n_subgroups_all);
  V_halos     =(REAL   *)SID_malloc(sizeof(REAL)*n_subgroups_all);
  x_min       =1e10;
  x_max       =0.;
  y_min       =1e10;
  y_max       =0.;
  z_min       =1e10;
  z_max       =0.;

  //   ... read global properties ...
  for(i_halo=0,n_halos=0;i_halo<n_groups_all;i_halo++){
    fread(&group_properties,sizeof(halo_properties_info),1,fp_groups);
    for(i_subhalo=0;i_subhalo<n_subgroups_group[i_halo];i_subhalo++){
      fread(&subgroup_properties,sizeof(halo_properties_info),1,fp_subgroups);
      if(subgroup_properties.n_particles>=n_particles_min && subgroup_properties.M_vir>0.){
        x_halos[n_halos]   =(REAL)subgroup_properties.position_MBP[0];
        y_halos[n_halos]   =(REAL)subgroup_properties.position_MBP[1];
        z_halos[n_halos]   =(REAL)subgroup_properties.position_MBP[2];
        r_halos[n_halos]   =(REAL)sqrt(pow(group_properties.position_MBP[0]-subgroup_properties.position_MBP[0],2.0)+
                                       pow(group_properties.position_MBP[1]-subgroup_properties.position_MBP[1],2.0)+
                                       pow(group_properties.position_MBP[2]-subgroup_properties.position_MBP[2],2.0));
        vx_halos_FoF[n_halos]=(REAL)group_properties.velocity_COM[0];
        vy_halos_FoF[n_halos]=(REAL)group_properties.velocity_COM[1];
        vz_halos_FoF[n_halos]=(REAL)group_properties.velocity_COM[2];
        vx_halos_sub[n_halos]=(REAL)subgroup_properties.velocity_COM[0];
        vy_halos_sub[n_halos]=(REAL)subgroup_properties.velocity_COM[1];
        vz_halos_sub[n_halos]=(REAL)subgroup_properties.velocity_COM[2];
        M_halos[n_halos]     =(double)subgroup_properties.M_vir;
        V_halos[n_halos]     =(REAL)subgroup_properties.V_max;
        x_min                =MIN(x_min,x_halos[n_halos]);
        x_max                =MAX(x_max,x_halos[n_halos]);
        y_min                =MIN(y_min,y_halos[n_halos]);
        y_max                =MAX(y_max,y_halos[n_halos]);
        z_min                =MIN(z_min,z_halos[n_halos]);
        z_max                =MAX(z_max,z_halos[n_halos]);
        n_halos++;
      }
    }
  }
  SID_log("x_range=%le->%le",SID_LOG_COMMENT,x_min,x_max);
  SID_log("y_range=%le->%le",SID_LOG_COMMENT,y_min,y_max);
  SID_log("z_range=%le->%le",SID_LOG_COMMENT,z_min,z_max);

  //   ... sort halos by V_max and close files.
  merge_sort(V_halos,(size_t)n_halos,&V_halos_index,SID_REAL,SORT_COMPUTE_INDEX,FALSE);
  fclose(fp_groups);
  fclose(fp_subgroups);
  SID_log("Done.",SID_LOG_CLOSE);

  // Write each halo grouping in turn ...
  SID_log("Writing %d groupings of halos...",SID_LOG_OPEN|SID_LOG_TIMER,n_groups);
  for(i_grouping=0,V_grouping=V_max_lo;i_grouping<n_groupings;i_grouping++,V_grouping+=delta_V_max){
    SID_log("Generating grouping %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_grouping+1,n_groupings);
    sprintf(filename_out,"%s_grouping_%03d.dat",filename_out_root,i_grouping);
    fp=fopen(filename_out,"w");

    //   ... set grouping ...
    i_halo=0;
    while(V_halos[V_halos_index[i_halo]]<V_grouping && i_halo<(n_halos-1)) i_halo++;
    i_halo=MAX(0,i_halo-n_halos_per_grouping/2);

    //   ... write grouping to file.
    SID_log("Halo index range=%d -> %d (out of %d)",SID_LOG_COMMENT,i_halo,i_halo+n_halos_per_grouping-1,n_halos);
    V_min=V_halos[V_halos_index[i_halo]];
    V_med=V_halos[V_halos_index[i_halo+n_halos_per_grouping/2-1]];
    V_max=V_halos[V_halos_index[i_halo+n_halos_per_grouping-1]];
    SID_log("n_halos=%d", SID_LOG_COMMENT,n_halos_per_grouping);
    SID_log("offset =%d", SID_LOG_COMMENT,i_halo);
    SID_log("V_min  =%le",SID_LOG_COMMENT,V_min);
    SID_log("V_med  =%le",SID_LOG_COMMENT,V_med);
    SID_log("V_max  =%le",SID_LOG_COMMENT,V_max);
    fprintf(fp,"# Grouping #%03d for %s\n",i_grouping,filename_out_root);
    fprintf(fp,"# V_min  =%le km/s\n",V_min);
    fprintf(fp,"# V_med  =%le km/s\n",V_med);
    fprintf(fp,"# V_max  =%le km/s\n",V_max);
    fprintf(fp,"# Columns:  1) Mass    [M_sol]\n");
    fprintf(fp,"#           2) V_max   [km/s]\n");
    fprintf(fp,"#           3) x       [Mpc/h]\n");
    fprintf(fp,"#           4) y       [Mpc/h]\n");
    fprintf(fp,"#           5) z       [Mpc/h]\n");
    fprintf(fp,"#           6) r_sys   [Mpc/h]\n");
    fprintf(fp,"#           7) v_x_sub [km/s]\n");
    fprintf(fp,"#           8) v_y_sub [km/s]\n");
    fprintf(fp,"#           9) v_z_sub [km/s]\n");
    fprintf(fp,"#          10) v_x_FoF [km/s]\n");
    fprintf(fp,"#          11) v_y_FoF [km/s]\n");
    fprintf(fp,"#          12) v_z_FoF [km/s]\n");
    fprintf(fp,"#          13) v_x_sys [km/s]\n");
    fprintf(fp,"#          14) v_y_sys [km/s]\n");
    fprintf(fp,"#          15) v_z_sys [km/s]\n");
    for(j_halo=0;j_halo<n_halos_per_grouping && i_halo<n_halos;i_halo++,j_halo++)
      fprintf(fp,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
                 M_halos[V_halos_index[i_halo]],
                 V_halos[V_halos_index[i_halo]],
                 x_halos[V_halos_index[i_halo]],
                 y_halos[V_halos_index[i_halo]],
                 z_halos[V_halos_index[i_halo]],
                 r_halos[V_halos_index[i_halo]],
                 vx_halos_sub[V_halos_index[i_halo]],
                 vy_halos_sub[V_halos_index[i_halo]],
                 vz_halos_sub[V_halos_index[i_halo]],
                 vx_halos_FoF[V_halos_index[i_halo]],
                 vy_halos_FoF[V_halos_index[i_halo]],
                 vz_halos_FoF[V_halos_index[i_halo]],
                 vx_halos_sub[V_halos_index[i_halo]]-vx_halos_FoF[V_halos_index[i_halo]],
                 vy_halos_sub[V_halos_index[i_halo]]-vy_halos_FoF[V_halos_index[i_halo]],
                 vz_halos_sub[V_halos_index[i_halo]]-vz_halos_FoF[V_halos_index[i_halo]]);

    fclose(fp);
    SID_log("Done.",SID_LOG_CLOSE);
  }
  SID_log("Done.",SID_LOG_CLOSE);
  
  // Clean-up
  SID_log("Cleaning-up ...",SID_LOG_OPEN);
  SID_free(SID_FARG x_halos);
  SID_free(SID_FARG y_halos);
  SID_free(SID_FARG z_halos);
  SID_free(SID_FARG r_halos);
  SID_free(SID_FARG vx_halos_FoF);
  SID_free(SID_FARG vy_halos_FoF);
  SID_free(SID_FARG vz_halos_FoF);
  SID_free(SID_FARG vx_halos_sub);
  SID_free(SID_FARG vy_halos_sub);
  SID_free(SID_FARG vz_halos_sub);
  SID_free(SID_FARG M_halos);
  SID_free(SID_FARG V_halos);
  SID_free(SID_FARG V_halos_index);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
