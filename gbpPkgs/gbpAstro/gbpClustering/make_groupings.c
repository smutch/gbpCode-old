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

#define MAKE_GROUPINGS_CENTRALS_ONLY  1
#define MAKE_GROUPINGS_USE_BIAS_MODEL 2

int main(int argc, char *argv[]){
  double  a_start;
  int     n_species;
  int     i,j,i_x,i_y,i_z;
  int     i_jack;
  int     n_jack_total;
  char    species_name[256];
  double  Omega_Lambda;
  double  Omega_M;
  double  Omega_b;
  double  f_gas;
  double  Omega_k;
  double  sigma_8;
  double  n_spec;
  int     i_rank;
  int     i_group;
  int     n_groups_all;
  int     n_subgroups_all;
  int     n_groups;
  size_t  n_group_local;
  char    filename_groups_in[MAX_FILENAME_LENGTH];
  char    filename_subgroups_in[MAX_FILENAME_LENGTH];
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
  field_info      FFT;
  plist_info      plist;
  FILE           *fp;
  FILE           *fp_groups;
  FILE           *fp_subgroups;
  FILE           *fp_stats;
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
  GBPREAL *x_halos;
  GBPREAL *y_halos;
  GBPREAL *z_halos;
  GBPREAL *r_halos;
  GBPREAL *vx_halos_FoF;
  GBPREAL *vy_halos_FoF;
  GBPREAL *vz_halos_FoF;
  GBPREAL *vx_halos_sub;
  GBPREAL *vy_halos_sub;
  GBPREAL *vz_halos_sub;
  double *M_halos;
  double *M_FoF;
  size_t *V_halos_index;
  GBPREAL *x_group;
  GBPREAL *y_group;
  GBPREAL *z_group;
  GBPREAL *vx_group;
  GBPREAL *vy_group;
  GBPREAL *vz_group;
  GBPREAL *V_halos;
  GBPREAL *V_FoF;
  int  *group_index;
  int  *subgroup_index;
  int  *i_FoF;
  int  *n_FoF;
  GBPREAL  x_min,x_max;
  GBPREAL  y_min,y_max;
  GBPREAL  z_min,z_max;
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
  size_t  n_random,n_random_local;
  int     i_random,j_random;
  RNG_info  RNG;
  int       seed=1327621;
  GBPREAL     *x_random;
  GBPREAL     *y_random;
  GBPREAL     *z_random;
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
  GBPREAL     V_grouping,V_max_lo,V_max_hi,delta_V_max;
  int         mode;

  // Initialization -- MPI etc.
  char filename_SSimPL_root[MAX_FILENAME_LENGTH];
  char filename_halos_version[MAX_FILENAME_LENGTH];
  SID_init(&argc,&argv,NULL,NULL);
  strcpy(filename_SSimPL_root,       argv[1]);
  strcpy(filename_halos_version,     argv[2]);
  strcpy(filename_out_root,          argv[3]);
  i_snap              =(int)    atoi(argv[4]);
  n_halos_per_grouping=(int)    atoi(argv[5]);
  n_particles_min     =(int)    atoi(argv[6]);
  V_max_lo            =(GBPREAL)atof(argv[7]);
  V_max_hi            =(GBPREAL)atof(argv[8]);
  n_groupings         =(int)    atoi(argv[9]);
  mode                =(int)    atoi(argv[10]);
  delta_V_max         =(V_max_hi-V_max_lo)/(double)(n_groupings-1);

  // Set some filename roots
  char filename_groups_root[MAX_FILENAME_LENGTH];
  char filename_cat_root[MAX_FILENAME_LENGTH];
  sprintf(filename_groups_root,"%s/groups/%s",  filename_SSimPL_root,filename_halos_version);
  sprintf(filename_cat_root,   "%s/catalogs/%s",filename_SSimPL_root,filename_halos_version);

  SID_log("Producing (up to) %d groupings of halos...",SID_LOG_OPEN,n_groupings);

  // Parse mode flag and report results
  double redshift;
  double redshift_norm;
  int    flag_centrals_only;
  int    flag_use_bias_model;
  flag_centrals_only =check_mode_for_flag(mode,MAKE_GROUPINGS_CENTRALS_ONLY);
  flag_use_bias_model=check_mode_for_flag(mode,MAKE_GROUPINGS_USE_BIAS_MODEL);
  if(flag_centrals_only)
     SID_log("USING CENTRALS ONLY.",SID_LOG_COMMENT);
  else
     SID_log("USING SUBSTRUCTURE.",SID_LOG_COMMENT);
  if(flag_use_bias_model){
     redshift     =(double)atof(argv[11]);
     redshift_norm=(double)atof(argv[12]);
     SID_log("USING THK BIAS MODEL TO SET NUMBER DENSITIES.",SID_LOG_COMMENT);
  }

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

  //   ... open files ...
  fp_catalog_info fp_properties_groups;
  fp_catalog_info fp_properties_subgroups;
  fopen_catalog(filename_cat_root,
                i_snap,
                READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES,
                &fp_properties_groups);
  fopen_catalog(filename_cat_root,
                i_snap,
                READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES,
                &fp_properties_subgroups);
  
  n_groups_all   =fp_properties_groups.n_halos_total;
  n_subgroups_all=fp_properties_subgroups.n_halos_total;

  if(n_groups_all   !=n_groups)    SID_trap_error("There's a mismatch in the number of groups (ie. %d!=%d)",   ERROR_LOGIC,n_groups_all,n_groups);
  if(n_subgroups_all!=n_subgroups) SID_trap_error("There's a mismatch in the number of subgroups (ie. %d!=%d)",ERROR_LOGIC,n_subgroups_all,n_subgroups);

  SID_log("%d groups and %d subgroups...",SID_LOG_CONTINUE,n_groups_all,n_subgroups_all);

  //   ... allocate arrays ...
  int n_halos_max;
  n_halos_max   =MAX(n_groups_all,n_subgroups_all);
  x_halos       =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_max);
  y_halos       =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_max);
  z_halos       =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_max);
  r_halos       =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_max);
  vx_halos_FoF  =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_max);
  vy_halos_FoF  =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_max);
  vz_halos_FoF  =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_max);
  vx_halos_sub  =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_max);
  vy_halos_sub  =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_max);
  vz_halos_sub  =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_max);
  M_halos       =(double  *)SID_malloc(sizeof(double) *n_halos_max);
  M_FoF         =(double  *)SID_malloc(sizeof(double) *n_halos_max);
  V_halos       =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_max);
  V_FoF         =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_max);
  i_FoF         =(int     *)SID_malloc(sizeof(int)    *n_halos_max);
  group_index   =(int     *)SID_malloc(sizeof(int)    *n_halos_max);
  subgroup_index=(int     *)SID_malloc(sizeof(int)    *n_halos_max);
  n_FoF         =(int     *)SID_malloc(sizeof(int)    *n_halos_max);
  x_min         =1e10;
  x_max         =0.;
  y_min         =1e10;
  y_max         =0.;
  z_min         =1e10;
  z_max         =0.;

  //   ... read global properties ...
  int k_subgroup=0;
  for(i_halo=0,n_halos=0;i_halo<n_groups_all;i_halo++){
    fread_catalog_raw(&fp_properties_groups,&group_properties,NULL,i_halo);
    if(flag_centrals_only){
      if(group_properties.n_particles>=n_particles_min && group_properties.M_vir>0.){
        if(!flag_centrals_only || (flag_centrals_only && i_subhalo==0)){
           group_index[n_halos]   =i_halo;
           group_index[n_halos]   =i_halo;
           i_FoF[n_halos]         =i_subhalo;
           n_FoF[n_halos]         =n_subgroups_group[i_halo];
           x_halos[n_halos]       =(GBPREAL)group_properties.position_MBP[0];
           y_halos[n_halos]       =(GBPREAL)group_properties.position_MBP[1];
           z_halos[n_halos]       =(GBPREAL)group_properties.position_MBP[2];
           r_halos[n_halos]       =(GBPREAL)sqrt(pow(group_properties.position_MBP[0]-group_properties.position_MBP[0],2.0)+
                                                 pow(group_properties.position_MBP[1]-group_properties.position_MBP[1],2.0)+
                                                 pow(group_properties.position_MBP[2]-group_properties.position_MBP[2],2.0));
           vx_halos_FoF[n_halos]=(GBPREAL)group_properties.velocity_COM[0];
           vy_halos_FoF[n_halos]=(GBPREAL)group_properties.velocity_COM[1];
           vz_halos_FoF[n_halos]=(GBPREAL)group_properties.velocity_COM[2];
           vx_halos_sub[n_halos]=(GBPREAL)group_properties.velocity_COM[0];
           vy_halos_sub[n_halos]=(GBPREAL)group_properties.velocity_COM[1];
           vz_halos_sub[n_halos]=(GBPREAL)group_properties.velocity_COM[2];
           M_halos[n_halos]     =(double)group_properties.M_vir;
           M_FoF[n_halos]       =(double)group_properties.M_vir;
           V_halos[n_halos]     =(GBPREAL)group_properties.V_max;
           V_FoF[n_halos]       =(GBPREAL)group_properties.V_max;
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
    else{
      for(i_subhalo=0;i_subhalo<n_subgroups_group[i_halo];i_subhalo++,k_subgroup++){
        fread_catalog_raw(&fp_properties_subgroups,&subgroup_properties,NULL,k_subgroup);
        if(subgroup_properties.n_particles>=n_particles_min && subgroup_properties.M_vir>0.){
          if(!flag_centrals_only || (flag_centrals_only && i_subhalo==0)){
             group_index[n_halos]   =i_halo;
             subgroup_index[n_halos]=k_subgroup;
             i_FoF[n_halos]         =i_subhalo;
             n_FoF[n_halos]         =n_subgroups_group[i_halo];
             x_halos[n_halos]       =(GBPREAL)subgroup_properties.position_MBP[0];
             y_halos[n_halos]       =(GBPREAL)subgroup_properties.position_MBP[1];
             z_halos[n_halos]       =(GBPREAL)subgroup_properties.position_MBP[2];
             r_halos[n_halos]       =(GBPREAL)sqrt(pow(group_properties.position_MBP[0]-subgroup_properties.position_MBP[0],2.0)+
                                                   pow(group_properties.position_MBP[1]-subgroup_properties.position_MBP[1],2.0)+
                                                   pow(group_properties.position_MBP[2]-subgroup_properties.position_MBP[2],2.0));
             vx_halos_FoF[n_halos]=(GBPREAL)group_properties.velocity_COM[0];
             vy_halos_FoF[n_halos]=(GBPREAL)group_properties.velocity_COM[1];
             vz_halos_FoF[n_halos]=(GBPREAL)group_properties.velocity_COM[2];
             vx_halos_sub[n_halos]=(GBPREAL)subgroup_properties.velocity_COM[0];
             vy_halos_sub[n_halos]=(GBPREAL)subgroup_properties.velocity_COM[1];
             vz_halos_sub[n_halos]=(GBPREAL)subgroup_properties.velocity_COM[2];
             M_halos[n_halos]     =(double)subgroup_properties.M_vir;
             M_FoF[n_halos]       =(double)group_properties.M_vir;
             V_halos[n_halos]     =(GBPREAL)subgroup_properties.V_max;
             V_FoF[n_halos]       =(GBPREAL)group_properties.V_max;
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
    }
  }
  SID_log("x_range=%le->%le",SID_LOG_COMMENT,x_min,x_max);
  SID_log("y_range=%le->%le",SID_LOG_COMMENT,y_min,y_max);
  SID_log("z_range=%le->%le",SID_LOG_COMMENT,z_min,z_max);

  //   ... sort halos by V_max and close files.
  merge_sort(V_halos,(size_t)n_halos,&V_halos_index,SID_REAL,SORT_COMPUTE_INDEX,FALSE);
  fclose_catalog(&fp_properties_groups);
  fclose_catalog(&fp_properties_subgroups);
  SID_log("Done.",SID_LOG_CLOSE);

  // Initialize cosmology if we're using the bias model to set number densities
  cosmo_info *cosmo=NULL;
  if(flag_use_bias_model){
     char filename_cosmology[MAX_FILENAME_LENGTH];
     sprintf(filename_cosmology,"%s/run",filename_SSimPL_root);
     read_gbpCosmo_file(&cosmo,filename_cosmology);
     init_sigma_M(&cosmo,PSPEC_LINEAR_TF,PSPEC_ALL_MATTER);
  }

  // Write each halo grouping in turn ...
  int i_column;
  SID_log("Writing %d groupings of halos...",SID_LOG_OPEN|SID_LOG_TIMER,n_groupings);
  sprintf(filename_out,"%s_stats.dat",filename_out_root);
  fp_stats=fopen(filename_out,"w");
  for(i_grouping=0,V_grouping=V_max_lo;i_grouping<=n_groupings;i_grouping++,V_grouping+=delta_V_max){

    // Set the number density for this grouping
    int    n_halos_this_grouping;
    double bias;
    if(flag_use_bias_model){
       double delta_c=1.686;
       double b_z_norm;
       b_z_norm=linear_growth_factor(redshift_norm,cosmo)/linear_growth_factor(redshift,cosmo);
       bias    =bias_model(V_grouping*1e3,delta_c,redshift,&cosmo,BIAS_MODEL_TRK|BIAS_MODEL_VMAX_ORDINATE);
       n_halos_this_grouping=(int)((b_z_norm*(double)n_halos_per_grouping)/(double)bias);
    }
    else
       n_halos_this_grouping=n_halos_per_grouping;

    //   ... set grouping.  If this is the last iteration, do all halos ...
    fp=NULL;
    if(i_grouping<n_groupings){
       i_halo=0;
       while(V_halos[V_halos_index[i_halo]]<V_grouping && i_halo<(n_halos-1)) i_halo++;
       if(i_halo>0){
          if((V_halos[V_halos_index[i_halo-1]]-V_grouping)<(V_halos[V_halos_index[i_halo]]-V_grouping)) i_halo--;
       }
       i_halo=MAX(0,i_halo-n_halos_this_grouping/2);

       // If this grouping pushes past the end of the list, then 
       //    skip this grouping and go ahead with the final 'all' list
       if((i_halo+n_halos_this_grouping)>=n_halos){
          SID_log("Out of halos.  Grouping %d not generated.",SID_LOG_COMMENT,i_grouping+1);
          i_grouping=n_groupings;
          i_halo=0;
          n_halos_this_grouping=n_halos;
       }
    }
    else{
       i_halo=0;
       n_halos_this_grouping=n_halos;
    }

    if(i_grouping==n_groupings){
       SID_log("Generating catalog of all halos...",SID_LOG_OPEN|SID_LOG_TIMER,i_grouping+1,n_groupings);
       sprintf(filename_out,"%s_all.dat",filename_out_root);
       fp=fopen(filename_out,"w");
    }
    else{
       SID_log("Generating grouping %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_grouping+1,n_groupings);
       if(flag_use_bias_model)
          SID_log("Using bias=%.3lf for this grouping (V_max=%.3lf km/s).",SID_LOG_COMMENT,bias,V_grouping);
       sprintf(filename_out,"%s_grouping_%03d.dat",filename_out_root,i_grouping);
       fp=fopen(filename_out,"w");
    }

    //   ... write halos to file.
    SID_log("Halo index range=%d -> %d (out of %d)",SID_LOG_COMMENT,i_halo,i_halo+n_halos_this_grouping-1,n_halos);
    V_min=V_halos[V_halos_index[i_halo]];
    V_med=V_halos[V_halos_index[i_halo+n_halos_this_grouping/2]];
    V_max=V_halos[V_halos_index[i_halo+n_halos_this_grouping-1]];
    M_min=M_halos[V_halos_index[i_halo]];
    M_med=M_halos[V_halos_index[i_halo+n_halos_this_grouping/2]];
    M_max=M_halos[V_halos_index[i_halo+n_halos_this_grouping-1]];

    SID_log("n_halos=%d", SID_LOG_COMMENT,n_halos_this_grouping);
    SID_log("offset =%d", SID_LOG_COMMENT,i_halo);
    SID_log("V_min  =%le",SID_LOG_COMMENT,V_min);
    SID_log("V_med  =%le",SID_LOG_COMMENT,V_med);
    SID_log("V_max  =%le",SID_LOG_COMMENT,V_max);
    fprintf(fp,"# Grouping #%03d for %s\n",i_grouping,filename_out_root);
    fprintf(fp,"# V_min  =%le km/s\n",V_min);
    fprintf(fp,"# V_med  =%le km/s\n",V_med);
    fprintf(fp,"# V_max  =%le km/s\n",V_max);
    i_column=1;
    fprintf(fp,"# Columns: (%02d) Snapshot Number\n",             i_column++);
    fprintf(fp,"#          (%02d) Group Index\n",                 i_column++);
    fprintf(fp,"#          (%02d) Subgroup Index\n",              i_column++);
    fprintf(fp,"#          (%02d) Rank order in FoF group\n",     i_column++);
    fprintf(fp,"#          (%02d) Number of halos in FoF group\n",i_column++);
    fprintf(fp,"#          (%02d) Mass sub  [M_sol/h]\n",         i_column++);
    fprintf(fp,"#          (%02d) V_max sub [km/s]\n",            i_column++);
    fprintf(fp,"#          (%02d) Mass FoF  [M_sol/h]\n",         i_column++);
    fprintf(fp,"#          (%02d) V_max FoF [km/s]\n",            i_column++);
    fprintf(fp,"#          (%02d) x         [Mpc/h]\n",i_column++);
    fprintf(fp,"#          (%02d) y         [Mpc/h]\n",i_column++);
    fprintf(fp,"#          (%02d) z         [Mpc/h]\n",i_column++);
    fprintf(fp,"#          (%02d) r_sub     [Mpc/h]\n",i_column++);
    fprintf(fp,"#          (%02d) v_x_sub   [km/s]\n", i_column++);
    fprintf(fp,"#          (%02d) v_y_sub   [km/s]\n", i_column++);
    fprintf(fp,"#          (%02d) v_z_sub   [km/s]\n", i_column++);
    fprintf(fp,"#          (%02d) v_x_FoF   [km/s]\n", i_column++);
    fprintf(fp,"#          (%02d) v_y_FoF   [km/s]\n", i_column++);
    fprintf(fp,"#          (%02d) v_z_FoF   [km/s]\n", i_column++);
    fprintf(fp,"#          (%02d) v_x_vir   [km/s]\n", i_column++);
    fprintf(fp,"#          (%02d) v_y_vir   [km/s]\n", i_column++);
    fprintf(fp,"#          (%02d) v_z_vir   [km/s]\n", i_column++);
    for(j_halo=0;j_halo<n_halos_this_grouping && i_halo<n_halos;i_halo++,j_halo++)
      fprintf(fp,"%d %d %d %d %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
                 i_snap,
                 group_index[V_halos_index[i_halo]],
                 subgroup_index[V_halos_index[i_halo]],
                 i_FoF[V_halos_index[i_halo]],
                 n_FoF[V_halos_index[i_halo]],
                 M_halos[V_halos_index[i_halo]],
                 V_halos[V_halos_index[i_halo]],
                 M_FoF[V_halos_index[i_halo]],
                 V_FoF[V_halos_index[i_halo]],
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

    // ... write statistics to file ...
    if(i_grouping==0){
       i_column=0;
       fprintf(fp_stats,"# Stats for groupings {%s*}\n",filename_out_root);
       fprintf(fp_stats,"# Column: (%d) grouping ID\n",   i_column++);
       fprintf(fp_stats,"#         (%d) n_halos\n",       i_column++);
       fprintf(fp_stats,"#         (%d) M_sub h^-1 [M_sol] (min)\n",   i_column++);
       fprintf(fp_stats,"#         (%d) M_sub h^-1 [M_sol] (median)\n",i_column++);
       fprintf(fp_stats,"#         (%d) M_sub h^-1 [M_sol] (max)\n",   i_column++);
       fprintf(fp_stats,"#         (%d) V_sub [km/s]       (min)\n",   i_column++);
       fprintf(fp_stats,"#         (%d) V_sub [km/s]       (median)\n",i_column++);
       fprintf(fp_stats,"#         (%d) V_sub [km/s]       (max)\n",   i_column++);
    }
    if(i_grouping<n_groupings)
       fprintf(fp_stats,"%03d %8d %le %le %le %le %le %le\n",i_grouping,n_halos_this_grouping,M_min,M_med,M_max,V_min,V_med,V_max);
    else
       fprintf(fp_stats,"all %8d %le %le %le %le %le %le\n",n_halos_this_grouping,M_min,M_med,M_max,V_min,V_med,V_max);
    SID_log("Done.",SID_LOG_CLOSE);
  }
  fclose(fp_stats);
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
  SID_free(SID_FARG M_FoF);
  SID_free(SID_FARG V_halos);
  SID_free(SID_FARG V_FoF);
  SID_free(SID_FARG i_FoF);
  SID_free(SID_FARG group_index);
  SID_free(SID_FARG subgroup_index);
  SID_free(SID_FARG n_FoF);
  SID_free(SID_FARG V_halos_index);
  if(cosmo!=NULL)
     free_cosmo(&cosmo);

  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
