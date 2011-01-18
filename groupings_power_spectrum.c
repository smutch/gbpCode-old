#define  _MAIN
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <common.h>
#include <groups.h>
#include <sph.h>
#include <clustering.h>

int main(int argc, char *argv[]){
  double  a_start;
  int     n_species;
  int     i,j,i_x,i_y,i_z;
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
  int     i_group;
  int     n_groups;
  size_t  n_group;
  char    filename_in[MAX_FILE_LENGTH];
  char    filename_in_model[MAX_FILE_LENGTH];
  char    filename_out_1D[MAX_FILE_LENGTH];
  char    filename_out_2D[MAX_FILE_LENGTH];
  char    filename_out_stats[MAX_FILE_LENGTH];
  char    filename_out_root[MAX_FILE_LENGTH];
  char    species_name_modifier[16];
  int     pspec_mode;
  int     i_compute;
  char    group_name[6];
  char    filename_TF[256];
  char    n_string[64];
  int             n[3];
  double          box_size;
  double          L[3];
  size_t          n_all;
  FILE           *fp_in;
  cosmo_info     *cosmo;
  FFT_info        FFT;
  plist_info      plist_header;
  plist_info      plist;
  FILE           *fp;
  FILE           *fp_1D;
  FILE           *fp_2D;
  int     flag_write_header;
  int     i_temp;
  int     n_temp;
  double *k_temp;
  double *kmin_temp;
  double *kmax_temp;
  int      grid_size;
  halo_properties_info  properties;
  int   n_halos_per_group,n_halos,n_halos_all,i_halo,j_halo;
  int   n_particles_min;
  REAL *x_halos;
  REAL *y_halos;
  REAL *z_halos;
  REAL *vx_halos;
  REAL *vy_halos;
  REAL *vz_halos;
  double *M_halos;
  size_t *M_halos_index;
  REAL *x_group;
  REAL *y_group;
  REAL *z_group;
  REAL *vx_group;
  REAL *vy_group;
  REAL *vz_group;
  REAL  x_min,x_max;
  REAL  y_min,y_max;
  REAL  z_min,z_max;
  int     n_particles;
  double  M_min,M_max,M_med;
  char   *line=NULL;
  int     line_length=0;
  int     n_spline,n_model;
  double *k_spline;
  double *P_spline;
  double *k_model;
  double *P_model;
  int     i_k;
  interp_info *interp_model;
  interp_info *interp_spline;
  double bias,bias_norm;
  int    i_file,n_files;
  int    mass_assignment_scheme;
  double k_min_1D;
  double k_max_1D;
  double dk_1D;
  double k_min_2D;
  double k_max_2D;
  double dk_2D;
  double *P_k_1D=NULL;
  double *dP_k_1D=NULL;
  int    *n_modes_1D=NULL;
  double *P_k_2D=NULL;
  double *dP_k_2D=NULL;
  int    *n_modes_2D=NULL;
  double *k_1D=NULL;
  int     n_k_1D,n_k_2D;

  // Initialization -- MPI etc.
  SID_init(&argc,&argv,NULL);
  strcpy(filename_in,      argv[1]);
  strcpy(filename_in_model,argv[2]);
  strcpy(filename_out_root,argv[3]);
  box_size         =(double)atof(argv[4]);
  grid_size        =(int)   atoi(argv[5]);
  n_halos_per_group=(int)   atoi(argv[6]);
  n_particles_min  =(int)   atoi(argv[7]);
  n_groups         =(int)   atoi(argv[8]);
  sprintf(filename_out_stats,"%s.stats",      filename_out_root);
  sprintf(filename_out_1D,   "%s.1D_pow_spec",filename_out_root);
  sprintf(filename_out_2D,   "%s.2D_pow_spec",filename_out_root);

  init_cosmo_std(&cosmo);

  mass_assignment_scheme=PSPEC_DIST_DWT20;
  //mass_assignment_scheme=PSPEC_DIST_NGP;

  // Read model
  SID_log("Reading model power spectrum from {%s}...",SID_LOG_OPEN,filename_in_model);
  fp     =fopen(filename_in_model,"r");
  n_model=count_lines_data(fp);
  k_model=(double  *)SID_malloc(sizeof(double)*n_model);
  P_model=(double  *)SID_malloc(sizeof(double)*n_model);
  for(i_k=0;i_k<n_model;i_k++){
    grab_next_line_data(fp,&line,&line_length);
    grab_double(line,1,&(k_model[i_k]));
    grab_double(line,2,&(P_model[i_k]));
  }
  fclose(fp);
  init_interpolate(k_model,P_model,(size_t)n_model,gsl_interp_cspline,&interp_model);
  SID_free(SID_FARG k_model);
  SID_free(SID_FARG P_model);
  SID_log("Done.",SID_LOG_CLOSE);

  // Set the k ranges
  k_min_1D=0.0;
  k_max_1D=0.3;
  dk_1D   =0.01;
  k_min_2D=0.0;
  k_max_2D=0.3;
  dk_2D   =0.02;

/*
  k_min_1D=0.2;
  k_max_1D=0.5;
  dk_1D   =0.01;
  k_min_2D=0.2;
  k_max_2D=1.0;
  dk_2D   =0.1;
*/

  // Create spline-fit of smooth P(k)
  SID_log("Generating smooth spline-fit to P(k)...",SID_LOG_OPEN);
  n_spline=12;
  k_spline=(double *)SID_malloc(n_spline*sizeof(double));  
  P_spline=(double *)SID_malloc(n_spline*sizeof(double));  
  k_spline[0] =0.01e-1;
  k_spline[1] =0.05e-1;
  k_spline[2] =0.10e-1;
  k_spline[3] =0.20e-1;
  k_spline[4] =0.25e-1;
  k_spline[5] =0.75e-1;
  k_spline[6] =1.25e-1;
  k_spline[7] =1.75e-1;
  k_spline[8] =2.25e-1;
  k_spline[9] =2.75e-1;
  k_spline[10]=3.25e-1;
  k_spline[11]=3.75e-1;
  for(i_k=0;i_k<n_spline;i_k++)
    P_spline[i_k]=interpolate(interp_model,k_spline[i_k]);
  init_interpolate(k_spline,P_spline,(size_t)n_spline,gsl_interp_cspline,&interp_spline);
  SID_log("Done.",SID_LOG_CLOSE);

  n[0]=grid_size;
  n[1]=n[0];
  n[2]=n[0];
  L[0]=box_size;
  L[1]=L[0];
  L[2]=L[0];

  // Initialization
  init_FFT(3,n,L,&FFT);
  init_plist(&plist,&(FFT.slab),GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);

  SID_log("k_min=%le k_max=%le dk=%le dk_in=%le",SID_LOG_COMMENT,FFT.k_field[0],FFT.k_field[FFT.n[0]-1],FFT.dk[0],dk_2D);

  // Read catalog
  SID_log("Reading catalog...",SID_LOG_OPEN);
  fp=fopen(filename_in,"r");
  fread(&i_file,     sizeof(int),1,fp);
  fread(&n_files,    sizeof(int),1,fp);
  fread(&n_halos,    sizeof(int),1,fp);
  fread(&n_halos_all,sizeof(int),1,fp);
  SID_log("%d halos...",SID_LOG_CONTINUE,n_halos);
  n_group =(size_t)n_halos;
  x_halos =(REAL  *)SID_malloc(sizeof(REAL)*n_halos_all);
  y_halos =(REAL  *)SID_malloc(sizeof(REAL)*n_halos_all);
  z_halos =(REAL  *)SID_malloc(sizeof(REAL)*n_halos_all);
  vx_halos=(REAL  *)SID_malloc(sizeof(REAL)*n_halos_all);
  vy_halos=(REAL  *)SID_malloc(sizeof(REAL)*n_halos_all);
  vz_halos=(REAL  *)SID_malloc(sizeof(REAL)*n_halos_all);
  x_group =(REAL  *)SID_malloc(sizeof(REAL)*n_halos_all);
  y_group =(REAL  *)SID_malloc(sizeof(REAL)*n_halos_all);
  z_group =(REAL  *)SID_malloc(sizeof(REAL)*n_halos_all);
  vx_group=(REAL  *)SID_malloc(sizeof(REAL)*n_halos_all);
  vy_group=(REAL  *)SID_malloc(sizeof(REAL)*n_halos_all);
  vz_group=(REAL  *)SID_malloc(sizeof(REAL)*n_halos_all);
  M_halos=(double *)SID_malloc(sizeof(double)*n_halos_all);
  x_min=1e10;
  x_max=0.;
  y_min=1e10;
  y_max=0.;
  z_min=1e10;
  z_max=0.;
  for(i_halo=0,n_halos=0;i_halo<n_halos_all;i_halo++){
    fread(&properties,sizeof(halo_properties_info),1,fp);
    if(properties.n_particles>=n_particles_min && properties.M_vir>0.){
      x_halos[n_halos] =(REAL)properties.position_MBP[0];
      y_halos[n_halos] =(REAL)properties.position_MBP[1];
      z_halos[n_halos] =(REAL)properties.position_MBP[2];
      vx_halos[n_halos]=(REAL)properties.velocity_COM[0];
      vy_halos[n_halos]=(REAL)properties.velocity_COM[1];
      vz_halos[n_halos]=(REAL)properties.velocity_COM[2];
      M_halos[n_halos] =(double)properties.M_vir;
      x_min=MIN(x_min,x_halos[n_halos]);
      x_max=MAX(x_max,x_halos[n_halos]);
      y_min=MIN(y_min,y_halos[n_halos]);
      y_max=MAX(y_max,y_halos[n_halos]);
      z_min=MIN(z_min,z_halos[n_halos]);
      z_max=MAX(z_max,z_halos[n_halos]);
      n_halos++;
    }
  }
  SID_log("x_range=%le->%le",SID_LOG_COMMENT,x_min,x_max);
  SID_log("y_range=%le->%le",SID_LOG_COMMENT,y_min,y_max);
  SID_log("z_range=%le->%le",SID_LOG_COMMENT,z_min,z_max);
  ADaPS_store(&(plist.data),(void *)x_group,  "x_halos", ADaPS_DEFAULT);
  ADaPS_store(&(plist.data),(void *)y_group,  "y_halos", ADaPS_DEFAULT);
  ADaPS_store(&(plist.data),(void *)z_group,  "z_halos", ADaPS_DEFAULT);
  ADaPS_store(&(plist.data),(void *)vx_group, "vx_halos",ADaPS_DEFAULT);
  ADaPS_store(&(plist.data),(void *)vy_group, "vy_halos",ADaPS_DEFAULT);
  ADaPS_store(&(plist.data),(void *)vz_group, "vz_halos",ADaPS_DEFAULT);
  fclose(fp);
  SID_log("Done.",SID_LOG_CLOSE);

  // Sort halos by mass
  SID_log("Sorting %d halos...",SID_LOG_OPEN,n_halos);
  merge_sort(M_halos,(size_t)n_halos,&M_halos_index,ADaM_DOUBLE,SORT_COMPUTE_INDEX,FALSE);
  SID_log("Done.",SID_LOG_CLOSE);

  // Initialize arrays used for the 1D power spectrum
  SID_log("Initializing arrays...",SID_LOG_OPEN);
  (n_k_1D)    =(int)(0.5+(k_max_1D-k_min_1D)/dk_1D);
  (k_1D)      =(double *)SID_malloc(sizeof(double)*(n_k_1D)); 
  (P_k_1D)    =(double *)SID_malloc(sizeof(double)*(n_k_1D)); 
  (dP_k_1D)   =(double *)SID_malloc(sizeof(double)*(n_k_1D)); 
  (n_modes_1D)=(int *)SID_malloc(sizeof(int)*(n_k_1D)); 

  // Initialize arrays used for the 2D power spectrum
  (n_k_2D)    =(int)(0.5+(k_max_2D-k_min_2D)/dk_1D);
  (P_k_2D)    =(double *)SID_malloc(sizeof(double)*(n_k_2D)*(n_k_2D));
  (dP_k_2D)   =(double *)SID_malloc(sizeof(double)*(n_k_2D)*(n_k_2D));
  (n_modes_2D)=(int *)SID_malloc(sizeof(int)*(n_k_2D)*(n_k_2D));
  SID_log("Done.",SID_LOG_CLOSE);

  // Process each grouping in turn
  SID_log("Computing power spectra of %d groupings of halos...",SID_LOG_OPEN|SID_LOG_TIMER,n_groups);
  fp   =fopen(filename_out_stats,"w");
  fp_1D=fopen(filename_out_1D,"w");
  fp_2D=fopen(filename_out_2D,"w");
  fprintf(fp,"# Power spectrum stats for %s\n",filename_out_1D);
  fprintf(fp,"# n_groupings= %d\n",n_groups);
  fprintf(fp,"# Columns:\n");
  fprintf(fp,"#   (1) # of halos\n");
  fprintf(fp,"#   (2) Minimum mass [M_sol/h]\n");
  fprintf(fp,"#   (3) Median  mass [M_sol/h]\n");
  fprintf(fp,"#   (4) Maximum mass [M_sol/h]\n");
  fprintf(fp,"#   (5) Bias (real-space)\n");
  fprintf(fp,"#   (6) Bias (v_x redshift-space)\n");
  fprintf(fp,"#   (7) Bias (v_y redshift-space)\n");
  fprintf(fp,"#   (8) Bias (v_z redshift-space)\n");
  for(i_group=0,flag_write_header=TRUE;i_group<n_groups;i_group++){
    SID_log("Processing grouping %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_group+1,n_groups);

    // Set group
    SID_log("Assembling grouping...",SID_LOG_OPEN);
    if(i_group==n_groups-1)
      i_halo=n_halos-1-n_halos_per_group;
    else
      i_halo =(int)((float)i_group*(float)(n_halos-n_halos_per_group)/(float)(n_groups));
    n_group=n_halos_per_group;
    SID_log("n_halos=%d", SID_LOG_COMMENT,n_group);
    M_max  =M_halos[M_halos_index[n_halos-i_halo-1]];
    M_med  =M_halos[M_halos_index[n_halos-i_halo-n_halos_per_group/2-1]];
    SID_log("offset =%d", SID_LOG_COMMENT,i_halo);
    for(j_halo=0;j_halo<n_halos_per_group;j_halo++,i_halo++){
      x_group[j_halo] =x_halos[M_halos_index[n_halos-i_halo-1]];
      y_group[j_halo] =y_halos[M_halos_index[n_halos-i_halo-1]];
      z_group[j_halo] =z_halos[M_halos_index[n_halos-i_halo-1]];
      vx_group[j_halo]=vx_halos[M_halos_index[n_halos-i_halo-1]];
      vy_group[j_halo]=vy_halos[M_halos_index[n_halos-i_halo-1]];
      vz_group[j_halo]=vz_halos[M_halos_index[n_halos-i_halo-1]];
    }
    M_min=M_halos[M_halos_index[n_halos-i_halo-1]];
    SID_log("M_min  =%le",SID_LOG_COMMENT,M_min);
    SID_log("M_med  =%le",SID_LOG_COMMENT,M_med);
    SID_log("M_max  =%le",SID_LOG_COMMENT,M_max);
    SID_log("Done.",SID_LOG_CLOSE);

    // Store group size
    ADaPS_store(&(plist.data),(void *)&n_group,"n_all_halos",ADaPS_SCALAR_SIZE_T);
    ADaPS_store(&(plist.data),(void *)&n_group,"n_halos",    ADaPS_SCALAR_SIZE_T);
    
    fprintf(fp,"%d %le %le %le",n_group,M_min,M_med,M_max);
    for(i_compute=0;i_compute<4;i_compute++){
      switch(i_compute){
      case 0:
        SID_log("Processing real-space ...",SID_LOG_OPEN|SID_LOG_TIMER);
        pspec_mode=PSPEC_DEFAULT;
        sprintf(species_name_modifier,"no_pec_%d", i_group);
        break;
      case 1:
        SID_log("Processing v_x redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
        pspec_mode=PSPEC_ADD_VX;
        sprintf(species_name_modifier,"with_vx_%d",i_group);
        break;
      case 2:
        SID_log("Processing v_y redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
        pspec_mode=PSPEC_ADD_VY;
        sprintf(species_name_modifier,"with_vy_%d",i_group);
        break;
      case 3:
        SID_log("Processing v_z redsift space...",SID_LOG_OPEN|SID_LOG_TIMER);
        pspec_mode=PSPEC_ADD_VZ;
        sprintf(species_name_modifier,"with_vz_%d",i_group);
        break;
      }
      
      // Compute power spectrum
      compute_power_spectrum(&plist,
                             &FFT,
                             mass_assignment_scheme,
                             pspec_mode,
                             cosmo,
                             "halos",
                             redshift,
                             n_k_1D,
                             k_min_1D,
                             k_max_1D,
                             n_k_2D,
                             k_min_2D,
                             k_max_2D,
                             k_1D,
                             P_k_1D,
                             dP_k_1D,
                             n_modes_1D,
                             P_k_2D,
                             dP_k_2D,
                             n_modes_2D);

      // Write power spectra
      if(SID.I_am_Master){
        // Write 1D power spectra
        if(flag_write_header){
          fwrite(&n_groups,sizeof(int),   1,     fp_1D);
          fwrite(&n_k_1D,  sizeof(int),   1,     fp_1D);
          fwrite(&k_min_1D,sizeof(double),1,     fp_1D);
          fwrite(&k_max_1D,sizeof(double),1,     fp_1D);
          fwrite(&dk_1D,   sizeof(double),1,     fp_1D);
          fwrite(k_1D,     sizeof(double),n_k_1D,fp_1D);
        }
        fwrite(P_k_1D,    sizeof(double),n_k_1D,fp_1D);
        fwrite(dP_k_1D,   sizeof(double),n_k_1D,fp_1D);
        fwrite(n_modes_1D,sizeof(int),   n_k_1D,fp_1D);

        // Write 2D power spectra
        if(flag_write_header){
          fwrite(&n_groups,sizeof(int),   1,fp_2D);
          fwrite(&n_k_2D,  sizeof(int),   1,fp_2D);          
          fwrite(&k_min_2D,sizeof(double),1,fp_2D);
          fwrite(&k_max_2D,sizeof(double),1,fp_2D);
          fwrite(&dk_2D,   sizeof(double),1,fp_2D);
        }
        fwrite(P_k_2D,    sizeof(double),n_k_2D*n_k_2D,fp_2D);
        fwrite(dP_k_2D,   sizeof(double),n_k_2D*n_k_2D,fp_2D);
        fwrite(n_modes_2D,sizeof(int),   n_k_2D*n_k_2D,fp_2D);

        flag_write_header=FALSE;
    
        // Compute (and print) bias
        SID_log("Computing bias (%d bins)...",SID_LOG_OPEN,n_k_1D);
        for(i_k=0,bias=0.,bias_norm=0.;i_k<n_k_1D;i_k++){
          bias     +=P_k_1D[i_k]/(dP_k_1D[i_k]*interpolate(interp_spline,k_1D[i_k]));
          bias_norm+=1./dP_k_1D[i_k];
        }
        bias=sqrt(bias/bias_norm);
        SID_log("bias   =%lf",SID_LOG_COMMENT,bias);
        fprintf(fp," %le",bias);
        SID_log("Done.",SID_LOG_CLOSE);
      }
      SID_log("Done.",SID_LOG_CLOSE);
    }
    fprintf(fp,"\n");

    SID_log("Done.",SID_LOG_CLOSE);
  }
  fclose(fp);
  fclose(fp_1D);
  fclose(fp_2D);
  SID_log("Done.",SID_LOG_CLOSE);
  
  // Clean-up
  SID_free(SID_FARG k_1D);
  SID_free(SID_FARG P_k_1D);
  SID_free(SID_FARG dP_k_1D);
  SID_free(SID_FARG n_modes_1D);
  SID_free(SID_FARG P_k_2D);
  SID_free(SID_FARG dP_k_2D);
  SID_free(SID_FARG n_modes_2D);
  free_plist(&plist);
  free_interpolate(&interp_model);
  free_interpolate(&interp_spline);

  SID_exit(ERROR_NONE);
}
