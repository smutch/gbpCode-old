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
  int     i,j,i_x,i_y,i_z;
  int     n_2D_total;
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
  int     i_grouping;
  int     i_grouping_start;
  int     i_grouping_stop;
  int     n_groupings;
  size_t  n_halos_local;
  char    filename_in[MAX_FILENAME_LENGTH];
  char    filename_in_root[MAX_FILENAME_LENGTH];
  char    filename_in_model[MAX_FILENAME_LENGTH];
  char    filename_out_1D[MAX_FILENAME_LENGTH];
  char    filename_out_2D[MAX_FILENAME_LENGTH];
  char    filename_out_root[MAX_FILENAME_LENGTH];
  int     i_compute;
  char    grouping_name[6];
  char    filename_TF[256];
  char    n_string[64];
  int             n[3];
  double          x_in,y_in,z_in,vx_in,vy_in,vz_in;
  double          box_size;
  double          L[3];
  size_t          n_all;
  FILE           *fp_in;
  cosmo_info     *cosmo;
  plist_info      plist_header;
  plist_info      plist;
  int     flag_write_header=TRUE;
  int   n_halos,i_halo,j_halo;
  int   n_particles_min;
  GBPREAL *x_halos;
  GBPREAL *y_halos;
  GBPREAL *z_halos;
  GBPREAL *vx_halos_sub;
  GBPREAL *vy_halos_sub;
  GBPREAL *vz_halos_sub;
  GBPREAL *vx_halos_FoF;
  GBPREAL *vy_halos_FoF;
  GBPREAL *vz_halos_FoF;
  GBPREAL *vx_halos_sys;
  GBPREAL *vy_halos_sys;
  GBPREAL *vz_halos_sys;
  GBPREAL *V_halos;
  double *M_halos;
  size_t *M_halos_index;
  GBPREAL  x_min,x_max;
  GBPREAL  y_min,y_max;
  GBPREAL  z_min,z_max;
  int     n_particles;
  double  M_min,M_max,M_med;
  double  V_min,V_max,V_med;
  char   *line=NULL;
  size_t  line_length=0;
  int     i_bin,j_bin;
  int    i_file,n_files;
  int    mass_assignment_scheme;
  int     n_1D,n_2D;
  int       seed=1327621;
  size_t      n_rank;
  int         flag_init_store;
  
  int        grid_size; 
  double     k_min_1D;
  double     k_max_1D; 
  double     dk_1D; 
  double     k_min_2D; 
  double     k_max_2D; 
  double     dk_2D; 
  field_info   FFT; 
  int        n_k_1D; 
  int        n_k_2D; 
  double    *k_1D;
  int       *n_modes_1D;
  int       *n_modes_2D;
  double   **P_k_1D_array;
  double   **dP_k_1D_array;
  double   **P_k_2D_array;
  double   **dP_k_2D_array;
  double    *P_k_1D; 
  double    *dP_k_1D; 
  double    *P_k_2D; 
  double    *dP_k_2D; 

  // Initialization -- MPI etc.
  SID_init(&argc,&argv,NULL);
  if(argc!=8)
    SID_trap_error("Incorrect syntax.",ERROR_SYNTAX);
  strcpy(filename_in_root, argv[1]);
  strcpy(filename_out_root,argv[2]);
  redshift         =(double)atof(argv[3]);
  grid_size        =(int)   atoi(argv[4]);
  box_size         =(double)atof(argv[5]);
  i_grouping_start =(int)   atoi(argv[6]);
  i_grouping_stop  =(int)   atoi(argv[7]);
  n_groupings      =i_grouping_stop-i_grouping_start+1;
  if(n_groupings<1)
      SID_trap_error("No groupings have been selected (you chose start=%d, stop=%d).",ERROR_LOGIC,i_grouping_start,i_grouping_stop);

  SID_log("Producing P(k)'s...",SID_LOG_OPEN);

  // Set the k ranges
  double k_Nyq;
  k_Nyq   =(TWO_PI*(double)grid_size/box_size)/2.;
  k_min_1D=0.02;
  dk_1D   =0.02;
  k_max_1D=dk_1D*(float)((int)(k_Nyq/dk_1D));
  k_min_2D=0.0;
  dk_2D   =0.02;
  k_max_2D=dk_2D*(float)((int)(k_Nyq/dk_2D));

  // Set mass assignment scheme
  mass_assignment_scheme=MAP2GRID_DIST_DWT20;

  // Initialization
  n[0]=grid_size;
  n[1]=n[0];
  n[2]=n[0];
  L[0]=box_size;
  L[1]=L[0];
  L[2]=L[0];
  init_cosmo_std(&cosmo);
  init_field(3,n,L,&FFT);

  // Initialize arrays
  SID_log("Initializing arrays...",SID_LOG_OPEN);
  n_k_1D       =(int)(0.5+(k_max_1D-k_min_1D)/dk_1D);
  n_k_2D       =(int)(0.5+(k_max_2D-k_min_2D)/dk_1D);
  k_1D         =(double  *)SID_malloc(sizeof(double)*(n_k_1D)); 
  n_modes_1D   =(int     *)SID_malloc(sizeof(int   )*(n_k_1D)); 
  n_modes_2D   =(int     *)SID_malloc(sizeof(int   )*(n_k_2D)*(n_k_2D));
  P_k_1D_array =(double **)SID_malloc(4*sizeof(double *)); 
  dP_k_1D_array=(double **)SID_malloc(4*sizeof(double *)); 
  P_k_2D_array =(double **)SID_malloc(4*sizeof(double *));
  dP_k_2D_array=(double **)SID_malloc(4*sizeof(double *));
  int i_array;
  for(i_array=0;i_array<4;i_array++){
      P_k_1D_array[i_array]    =(double *)SID_malloc(sizeof(double)*(n_k_1D)); 
      dP_k_1D_array[i_array]   =(double *)SID_malloc(sizeof(double)*(n_k_1D)); 
      P_k_2D_array[i_array]    =(double *)SID_malloc(sizeof(double)*(n_k_2D)*(n_k_2D));
      dP_k_2D_array[i_array]   =(double *)SID_malloc(sizeof(double)*(n_k_2D)*(n_k_2D));
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Process each grouping in turn
  for(i_grouping=i_grouping_start,flag_init_store=TRUE;i_grouping<=i_grouping_stop;i_grouping++){
    SID_log("Processing grouping #%03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_grouping);

    // Loop over ithe real-space and 3 redshift-space frames
    for(i_compute=0;i_compute<4;i_compute++){

      // Read catalog
      int pspec_mode;
      init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
      switch(i_compute){
      case 0:
        SID_log("Processing real-space ...",SID_LOG_OPEN|SID_LOG_TIMER);
        read_groupings(filename_in_root,i_grouping,&plist,READ_GROUPING_DEFAULT,&(FFT.slab));
        pspec_mode=READ_GROUPING_DEFAULT;
        P_k_1D    =P_k_1D_array[0];
        dP_k_1D   =dP_k_1D_array[0];
        P_k_2D    =P_k_2D_array[0];
        dP_k_2D   =dP_k_2D_array[0];
        break;
      case 1:
        SID_log("Processing v_x redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
        read_groupings(filename_in_root,i_grouping,&plist,READ_GROUPING_DEFAULT|READ_GROUPING_ADD_VX,&(FFT.slab),box_size,redshift,cosmo);
        pspec_mode=READ_GROUPING_DEFAULT|READ_GROUPING_ADD_VX;
        P_k_1D    =P_k_1D_array[1];
        dP_k_1D   =dP_k_1D_array[1];
        P_k_2D    =P_k_2D_array[1];
        dP_k_2D   =dP_k_2D_array[1];
        break;
      case 2:
        SID_log("Processing v_y redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
        read_groupings(filename_in_root,i_grouping,&plist,READ_GROUPING_DEFAULT|READ_GROUPING_ADD_VY,&(FFT.slab),box_size,redshift,cosmo);
        pspec_mode=READ_GROUPING_DEFAULT|READ_GROUPING_ADD_VY;
        P_k_1D    =P_k_1D_array[2];
        dP_k_1D   =dP_k_1D_array[2];
        P_k_2D    =P_k_2D_array[2];
        dP_k_2D   =dP_k_2D_array[2];
        break;
      case 3:
        SID_log("Processing v_z redsift space...",SID_LOG_OPEN|SID_LOG_TIMER);
        read_groupings(filename_in_root,i_grouping,&plist,READ_GROUPING_DEFAULT|READ_GROUPING_ADD_VZ,&(FFT.slab),box_size,redshift,cosmo);
        pspec_mode=READ_GROUPING_DEFAULT|READ_GROUPING_ADD_VZ;
        P_k_1D    =P_k_1D_array[3];
        dP_k_1D   =dP_k_1D_array[3];
        P_k_2D    =P_k_2D_array[3];
        dP_k_2D   =dP_k_2D_array[3];
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

      // Clean-up
      free_plist(&plist);

      SID_log("Done.",SID_LOG_CLOSE);
    } // Loop over 4 P(k)'s

    // Now that all 4 runs are done, let's write the results
    SID_log("Writing results...",SID_LOG_OPEN|SID_LOG_TIMER);
    if(SID.I_am_Master){
       // Set output filenames
       sprintf(filename_out_1D,"%s_grouping_%03d_1D_pspec.dat",filename_out_root,i_grouping);
       sprintf(filename_out_2D,"%s_grouping_%03d_2D_pspec.dat",filename_out_root,i_grouping);

       // Write 1-D Results
       FILE *fp_out;
       int   i_k;
       char  grouping_name[MAX_FILENAME_LENGTH];
       strcpy(grouping_name,filename_in_root);
       strip_path(grouping_name);
       fp_out=fopen(filename_out_1D,"w");
       fprintf(fp_out,"# 1D Power spectra for grouping #%03d of {%s}\n",i_grouping,grouping_name);
       fprintf(fp_out,"#\n");
       switch(mass_assignment_scheme){
          case MAP2GRID_DIST_DWT20:
             fprintf(fp_out,"# Mass assignment scheme: DWT20\n",i_grouping,grouping_name);
             break;
          case MAP2GRID_DIST_DWT12:
             fprintf(fp_out,"# Mass assignment scheme: DWT12\n",i_grouping,grouping_name);
             break;
          case MAP2GRID_DIST_NGP:
             fprintf(fp_out,"# Mass assignment scheme: NGP\n",i_grouping,grouping_name);
             break;
          case MAP2GRID_DIST_TSC:
             fprintf(fp_out,"# Mass assignment scheme: TSC\n",i_grouping,grouping_name);
             break;
       }
       fprintf(fp_out,"# Redshift:               %5.3lf\n",        redshift);
       fprintf(fp_out,"# Box size:               %9.3le [Mpc/h]\n",box_size);
       fprintf(fp_out,"# Grid size:              %d^3\n",          grid_size);
       int i_column;
       int i_run;
       fprintf(fp_out,"#\n");
       fprintf(fp_out,"# Column: (%02d) k [h/Mpc]\n",i_column++);
       fprintf(fp_out,"#         (%02d) # of Fourier modes in this bin\n",i_column++);
       for(i_run=0;i_run<4;i_run++){
          char run_name[128];
          switch(i_run){
             case 0:
                sprintf(run_name,"real");
                break;
             case 1:
                sprintf(run_name,"x-projected z");
                break;
             case 2:
                sprintf(run_name,"y-projected z");
                break;
             case 3:
                sprintf(run_name,"z-projected z");
                break;
          }
          fprintf(fp_out,"#         (%02d) P(k)  [h^2/Mpc^2];  %s-space\n",i_column++,run_name);
          fprintf(fp_out,"#         (%02d) dP(k) [h^2/Mpc^2];  %s-space\n",i_column++,run_name);
       }
       fprintf(fp_out,"#\n");
       for(i_k=0;i_k<n_k_1D;i_k++){
          fprintf(fp_out,"%9.3le %8d ",k_1D[i_k],n_modes_1D[i_k]);
          for(i_run=0;i_run<4;i_run++)
             fprintf(fp_out,"  %le %le",P_k_1D_array[i_run][i_k],dP_k_1D_array[i_run][i_k]);
          fprintf(fp_out,"\n");
       }
       fclose(fp_out);

       // Write 2D power spectra
       fp_out=fopen(filename_out_2D,"w");
       fwrite(&n_k_2D,   sizeof(int),   1,fp_out);
       fwrite(&k_min_2D, sizeof(double),1,fp_out);
       fwrite(&k_max_2D, sizeof(double),1,fp_out);
       fwrite(&dk_2D,    sizeof(double),1,fp_out);
       fwrite(n_modes_2D,sizeof(int),   n_k_2D*n_k_2D,fp_out);
       for(i_run=0;i_run<4;i_run++){
          fwrite(P_k_2D_array[i_run],    sizeof(double),n_k_2D*n_k_2D,fp_out);
          fwrite(dP_k_2D_array[i_run],   sizeof(double),n_k_2D*n_k_2D,fp_out);
       }
       fclose(fp_out);
    }
    SID_log("Done.",SID_LOG_CLOSE);

    SID_log("Done.",SID_LOG_CLOSE);
  } // Loop over groupings

  // Clean-up
  SID_log("Cleaning-up...",SID_LOG_OPEN);
  free_cosmo(&cosmo);
  free_field(&FFT);
  SID_free(SID_FARG k_1D);
  SID_free(SID_FARG n_modes_1D);
  SID_free(SID_FARG n_modes_2D);
  for(i_array=0;i_array<4;i_array++){
      SID_free(SID_FARG P_k_1D_array[i_array]);
      SID_free(SID_FARG dP_k_1D_array[i_array]);
      SID_free(SID_FARG P_k_2D_array[i_array]);
      SID_free(SID_FARG dP_k_2D_array[i_array]);
  }
  SID_free(SID_FARG P_k_1D_array);
  SID_free(SID_FARG dP_k_1D_array);
  SID_free(SID_FARG P_k_2D_array);
  SID_free(SID_FARG dP_k_2D_array);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
