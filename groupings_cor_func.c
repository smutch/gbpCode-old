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
  int     i_group;
  int     n_groups=1;
  size_t  n_group_local;
  char    filename_in[MAX_FILENAME_LENGTH];
  char    filename_in_model[MAX_FILENAME_LENGTH];
  char    filename_out_1D[MAX_FILENAME_LENGTH];
  char    filename_out_2D[MAX_FILENAME_LENGTH];
  char    filename_out_stats[MAX_FILENAME_LENGTH];
  char    filename_out_root[MAX_FILENAME_LENGTH];
  int     cfunc_mode;
  int     i_compute;
  char    group_name[6];
  char    filename_TF[256];
  char    n_string[64];
  int             n[3];
  int             x_column;
  int             y_column;
  int             z_column;
  int             vx_column;
  int             vy_column;
  int             vz_column;
  double          x_in,y_in,z_in,vx_in,vy_in,vz_in;
  double          box_size;
  double          L[3];
  size_t          n_all;
  FILE           *fp_in;
  cosmo_info     *cosmo;
  plist_info      plist_header;
  plist_info      plist;
  FILE           *fp;
  FILE           *fp_1D;
  FILE           *fp_2D;
  int     flag_write_header=TRUE;
  int     i_temp;
  size_t  n_group;
  int     n_temp;
  double *r_temp;
  double *kmin_temp;
  double *kmax_temp;
  int      n_jack;
  int   n_halos,i_halo,j_halo;
  int   n_particles_min;
  REAL *x_halos;
  REAL *y_halos;
  REAL *z_halos;
  REAL *vx_halos_sub;
  REAL *vy_halos_sub;
  REAL *vz_halos_sub;
  REAL *V_halos;
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
  double  V_min,V_max,V_med;
  char   *line=NULL;
  size_t  line_length=0;
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
  double *COVMTX_2D=NULL;
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
  int        *rank_own;
  int        *rank_own_index;
  size_t      n_rank;

  // Initialization -- MPI etc.
  SID_init(&argc,&argv,NULL);
  if(argc!=12)
    SID_trap_error("Incorrect syntax.",ERROR_SYNTAX);
  strcpy(filename_in,      argv[1]);
  strcpy(filename_out_root,argv[2]);
  redshift         =(double)atof(argv[3]);
  box_size         =(double)atof(argv[4]);
  n_jack           =(int)   atoi(argv[5]);
  i_grouping_start =(int)   atoi(argv[6]);
  i_grouping_stop  =(int)   atoi(argv[7]);
  sprintf(filename_out_stats,"%s.stats_cor_func",filename_out_root);
  sprintf(filename_out_1D,   "%s.1D_cor_func",   filename_out_root);
  sprintf(filename_out_2D,   "%s.2D_cor_func",   filename_out_root);

  SID_log("Producing correlation function...",SID_LOG_OPEN);

  // Set the k ranges
  flag_log_1D=FALSE;
  n_1D       =   20;
  r_min_l1D  =  0.1;
  r_max_1D   =200.0;
  dr_1D      =r_max_1D/(double)(n_1D-1);
  dlr_l1D     =(take_log10(r_max_1D)-take_log10(r_min_l1D))/(double)(n_1D-1);
  n_2D       =   25;
  r_min_2D   =  0.0;
  r_max_2D   = 50.0;
  dr_2D      =(r_max_2D-r_min_2D)/(double)(n_2D-1);
  F_random   =5;

  // Initialization
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  init_cosmo_std(&cosmo);
  init_RNG(&seed,&RNG,RNG_DEFAULT);

  // Initialize arrays used for the 1D correlation function
  SID_log("Initializing arrays...",SID_LOG_OPEN);
  CFUNC_l1D =(double *)SID_malloc(sizeof(double)*(n_1D)); 
  dCFUNC_l1D=(double *)SID_malloc(sizeof(double)*(n_1D)); 
  COVMTX_l1D=(double *)SID_malloc(sizeof(double)*(n_1D*n_1D)); 
  CFUNC_1D  =(double *)SID_malloc(sizeof(double)*(n_1D)); 
  dCFUNC_1D =(double *)SID_malloc(sizeof(double)*(n_1D)); 
  COVMTX_1D =(double *)SID_malloc(sizeof(double)*(n_1D*n_1D)); 

  // Initialize arrays used for the 2D correlation function
  n_2D_total=n_2D*n_2D;
  CFUNC_2D  =(double *)SID_malloc(sizeof(double)*(n_2D)*(n_2D));
  dCFUNC_2D =(double *)SID_malloc(sizeof(double)*(n_2D)*(n_2D));
  COVMTX_2D =(double *)SID_malloc(sizeof(double)*(n_2D_total*n_2D_total)); 

  flag_compute_RR=TRUE;
  n_jack_total   =n_jack*n_jack*n_jack;
  DD_l1D=(long long **)SID_malloc(sizeof(long long *)*(n_jack_total+1));
  DR_l1D=(long long **)SID_malloc(sizeof(long long *)*(n_jack_total+1));
  RR_l1D=(long long **)SID_malloc(sizeof(long long *)*(n_jack_total+1));
  DD_1D =(long long **)SID_malloc(sizeof(long long *)*(n_jack_total+1));
  DR_1D =(long long **)SID_malloc(sizeof(long long *)*(n_jack_total+1));
  RR_1D =(long long **)SID_malloc(sizeof(long long *)*(n_jack_total+1));
  DD_2D =(long long **)SID_malloc(sizeof(long long *)*(n_jack_total+1));
  DR_2D =(long long **)SID_malloc(sizeof(long long *)*(n_jack_total+1));
  RR_2D =(long long **)SID_malloc(sizeof(long long *)*(n_jack_total+1));
  for(i_jack=0;i_jack<=n_jack_total;i_jack++){
    DD_l1D[i_jack]=(long long *)SID_malloc(sizeof(long long)*n_1D);
    DR_l1D[i_jack]=(long long *)SID_malloc(sizeof(long long)*n_1D);
    RR_l1D[i_jack]=(long long *)SID_malloc(sizeof(long long)*n_1D);
    DD_1D[i_jack] =(long long *)SID_malloc(sizeof(long long)*n_1D);
    DR_1D[i_jack] =(long long *)SID_malloc(sizeof(long long)*n_1D);
    RR_1D[i_jack] =(long long *)SID_malloc(sizeof(long long)*n_1D);
    DD_2D[i_jack] =(long long *)SID_malloc(sizeof(long long)*n_2D*n_2D);
    DR_2D[i_jack] =(long long *)SID_malloc(sizeof(long long)*n_2D*n_2D);
    RR_2D[i_jack] =(long long *)SID_malloc(sizeof(long long)*n_2D*n_2D);
    for(i_bin=0;i_bin<n_1D;i_bin++){
      DD_l1D[i_jack][i_bin]=0;
      DR_l1D[i_jack][i_bin]=0;
      RR_l1D[i_jack][i_bin]=0;
      DD_1D[i_jack][i_bin] =0;
      DR_1D[i_jack][i_bin] =0;
      RR_1D[i_jack][i_bin] =0;
    }
    for(i_bin=0;i_bin<n_2D*n_2D;i_bin++){
      DD_2D[i_jack][i_bin]=0;
      DR_2D[i_jack][i_bin]=0;
      RR_2D[i_jack][i_bin]=0;
    }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Initialize files
  if(SID.I_am_Master){
    fp   =fopen(filename_out_stats,"w");
    fp_1D=fopen(filename_out_1D,"w");
    fp_2D=fopen(filename_out_2D,"w");
    fprintf(fp,"# Correlation function stats for %s\n",filename_out_1D);
    fprintf(fp,"# Columns:\n");
    fprintf(fp,"#   (1) # of objects\n");
    fprintf(fp,"#   (2) r_o [Mpc/h] (real-space)\n");
    fprintf(fp,"#   (3) r_X [Mpc/h] (real-space)\n");
    fprintf(fp,"#   (4) r_o [Mpc/h] (v_x redshift-space)\n");
    fprintf(fp,"#   (5) r_X [Mpc/h] (v_x redshift-space)\n");
    fprintf(fp,"#   (6) r_o [Mpc/h] (v_y redshift-space)\n");
    fprintf(fp,"#   (7) r_X [Mpc/h] (v_y redshift-space)\n");
    fprintf(fp,"#   (8) r_o [Mpc/h] (v_z redshift-space)\n");
    fprintf(fp,"#   (9) r_X [Mpc/h] (v_z redshift-space)\n");
  }

  // Process each grouping in turn
  for(i_grouping=i_grouping_start;i_grouping<=i_grouping_stop;i_grouping++){
    SID_log("Processing grouping #%03d...",SID_LOG_OPEN,i_grouping);

    // Read catalog
    SID_log("Reading catalog...",SID_LOG_OPEN);
    if(SID.I_am_Master){
      fp=fopen(filename_in,"r");
      n_group=(size_t)count_lines_data(fp);
      SID_log("%lld halos...",SID_LOG_CONTINUE,n_group);
      x_halos     =(REAL  *)SID_malloc(sizeof(REAL)*n_group);
      y_halos     =(REAL  *)SID_malloc(sizeof(REAL)*n_group);
      z_halos     =(REAL  *)SID_malloc(sizeof(REAL)*n_group);
      vx_halos_sub=(REAL  *)SID_malloc(sizeof(REAL)*n_group);
      vy_halos_sub=(REAL  *)SID_malloc(sizeof(REAL)*n_group);
      vz_halos_sub=(REAL  *)SID_malloc(sizeof(REAL)*n_group);
      vx_halos_FoF=(REAL  *)SID_malloc(sizeof(REAL)*n_group);
      vy_halos_FoF=(REAL  *)SID_malloc(sizeof(REAL)*n_group);
      vz_halos_FoF=(REAL  *)SID_malloc(sizeof(REAL)*n_group);
      vx_halos_sys=(REAL  *)SID_malloc(sizeof(REAL)*n_group);
      vy_halos_sys=(REAL  *)SID_malloc(sizeof(REAL)*n_group);
      vz_halos_sys=(REAL  *)SID_malloc(sizeof(REAL)*n_group);
      rank_own    =(int   *)SID_malloc(sizeof(int) *n_group);
      x_min=1e10;
      x_max=0.;
      y_min=1e10;
      y_max=0.;
      z_min=1e10;
      z_max=0.;
      for(i_halo=0;i_halo<n_group;i_halo++){
        grab_next_line_data(fp,&line,&line_length);
        grab_double(line,3, &x_in);
        grab_double(line,4, &y_in);
        grab_double(line,5, &z_in);
        x_halos[i_halo] =(REAL)x_in;
        y_halos[i_halo] =(REAL)y_in;
        z_halos[i_halo] =(REAL)z_in;
        x_min=MIN(x_min,x_halos[i_halo]);
        x_max=MAX(x_max,x_halos[i_halo]);
        y_min=MIN(y_min,y_halos[i_halo]);
        y_max=MAX(y_max,y_halos[i_halo]);
        z_min=MIN(z_min,z_halos[i_halo]);
        z_max=MAX(z_max,z_halos[i_halo]);
        grab_double(line,7,&vx_in);
        grab_double(line,8,&vy_in);
        grab_double(line,9,&vz_in);
        vx_halos_sub[i_halo]=(REAL)vx_in;
        vy_halos_sub[i_halo]=(REAL)vy_in;
        vz_halos_sub[i_halo]=(REAL)vz_in;
        grab_double(line,10,&vx_in);
        grab_double(line,11,&vy_in);
        grab_double(line,12,&vz_in);
        vx_halos_FoF[i_halo]=(REAL)vx_in;
        vy_halos_FoF[i_halo]=(REAL)vy_in;
        vz_halos_FoF[i_halo]=(REAL)vz_in;
        grab_double(line,13,&vx_in);
        grab_double(line,14,&vy_in);
        grab_double(line,15,&vz_in);
        vx_halos_sys[i_halo]=(REAL)vx_in;
        vy_halos_sys[i_halo]=(REAL)vy_in;
        vz_halos_sys[i_halo]=(REAL)vz_in;
      }
      SID_log("x_range=%le->%le",SID_LOG_COMMENT,x_min,x_max);
      SID_log("y_range=%le->%le",SID_LOG_COMMENT,y_min,y_max);
      SID_log("z_range=%le->%le",SID_LOG_COMMENT,z_min,z_max);
      fclose(fp);
    }
    SID_Bcast(&n_group,sizeof(size_t),MASTER_RANK);
    ADaPS_store(&(plist.data),(void *)x_group, "x_halos", ADaPS_DEFAULT);
    ADaPS_store(&(plist.data),(void *)y_group, "y_halos", ADaPS_DEFAULT);
    ADaPS_store(&(plist.data),(void *)z_group, "z_halos", ADaPS_DEFAULT);
    ADaPS_store(&(plist.data),(void *)vx_group,"vx_halos_sub",ADaPS_DEFAULT);
    ADaPS_store(&(plist.data),(void *)vy_group,"vy_halos_sub",ADaPS_DEFAULT);
    ADaPS_store(&(plist.data),(void *)vz_group,"vz_halos_sub",ADaPS_DEFAULT);
    SID_log("Done.",SID_LOG_CLOSE);

    // GENERATE RANDOMS
    n_random=(size_t)F_random*n_group;
    SID_log("Generating %d random halos...",SID_LOG_OPEN,n_random);
    for(i_rank=0,i_random=0,j_random=0;i_rank<SID.n_proc;i_rank++){
      if(i_rank==SID.n_proc-1)
        j_random=n_random-i_random;
      else
        j_random=(n_random-i_random)/(SID.n_proc-i_rank);
      if(SID.My_rank==i_rank)
        n_random_local=j_random;
      i_random+=j_random;
    }
    x_random=(REAL *)SID_malloc(sizeof(REAL)*n_random_local);
    y_random=(REAL *)SID_malloc(sizeof(REAL)*n_random_local);
    z_random=(REAL *)SID_malloc(sizeof(REAL)*n_random_local);
    for(i_random=0;i_random<n_random_local;i_random++){
      x_random[i_random]=random_number(&RNG)*box_size;
      if(x_random[i_random]>=box_size)
        x_random[i_random]-=box_size;
      y_random[i_random]=random_number(&RNG)*box_size;
      if(y_random[i_random]>=box_size)
        y_random[i_random]-=box_size;
      z_random[i_random]=random_number(&RNG)*box_size;
      if(z_random[i_random]>=box_size)
        z_random[i_random]-=box_size;
    }
    ADaPS_store(&(plist.data),(void *)&n_random,      "n_all_random",ADaPS_SCALAR_SIZE_T);
    ADaPS_store(&(plist.data),(void *)&n_random_local,"n_random",    ADaPS_SCALAR_SIZE_T);
    ADaPS_store(&(plist.data),(void *)x_random,       "x_random",    ADaPS_DEFAULT);
    ADaPS_store(&(plist.data),(void *)y_random,       "y_random",    ADaPS_DEFAULT);
    ADaPS_store(&(plist.data),(void *)z_random,       "z_random",    ADaPS_DEFAULT);
    SID_log("Done.",SID_LOG_CLOSE);

    SID_log("Assign rank ownership...",SID_LOG_OPEN|SID_LOG_TIMER);
    if(SID.I_am_Master){
      for(j_halo=0;j_halo<n_group;j_halo++,i_halo++)
        rank_own[j_halo]=(int)(random_number(&RNG)*SID.n_proc);
    }

    // Communicate with other ranks
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
      if(i_rank!=MASTER_RANK){
        if(SID.I_am_Master){
          for(j_halo=0,n_rank=0;j_halo<n_group;j_halo++){
            if(rank_own[j_halo]==i_rank){
              x_group[n_rank] =x_halos[j_halo];
              y_group[n_rank] =y_halos[j_halo];
              z_group[n_rank] =z_halos[j_halo];
              vx_group[n_rank]=vx_halos_sub[j_halo];
              vy_group[n_rank]=vy_halos_sub[j_halo];
              vz_group[n_rank]=vz_halos_sub[j_halo];
              n_rank++;
            }
          }
          SID_Send(&n_rank, 1,          SID_SIZE_T,i_rank,136789,SID.COMM_WORLD);
          SID_Send(x_group, (int)n_rank,SID_REAL,  i_rank,136790,SID.COMM_WORLD);
          SID_Send(y_group, (int)n_rank,SID_REAL,  i_rank,136791,SID.COMM_WORLD);
          SID_Send(z_group, (int)n_rank,SID_REAL,  i_rank,136792,SID.COMM_WORLD);
          SID_Send(vx_group,(int)n_rank,SID_REAL,  i_rank,136793,SID.COMM_WORLD);
          SID_Send(vy_group,(int)n_rank,SID_REAL,  i_rank,136794,SID.COMM_WORLD);
          SID_Send(vz_group,(int)n_rank,SID_REAL,  i_rank,136795,SID.COMM_WORLD);
        }
        if(SID.My_rank==i_rank){
          SID_Recv(&n_group_local,1,                 SID_SIZE_T,MASTER_RANK,136789,SID.COMM_WORLD);
          SID_Recv(x_group,       (int)n_group_local,SID_REAL,  MASTER_RANK,136790,SID.COMM_WORLD);
          SID_Recv(y_group,       (int)n_group_local,SID_REAL,  MASTER_RANK,136791,SID.COMM_WORLD);
          SID_Recv(z_group,       (int)n_group_local,SID_REAL,  MASTER_RANK,136792,SID.COMM_WORLD);
          SID_Recv(vx_group,      (int)n_group_local,SID_REAL,  MASTER_RANK,136793,SID.COMM_WORLD);
          SID_Recv(vy_group,      (int)n_group_local,SID_REAL,  MASTER_RANK,136794,SID.COMM_WORLD);
          SID_Recv(vz_group,      (int)n_group_local,SID_REAL,  MASTER_RANK,136795,SID.COMM_WORLD);
        }
      }
    }

    // Take care of the master rank
    if(SID.I_am_Master){
      for(j_halo=0,n_group_local=0;j_halo<n_group;j_halo++){
        if(rank_own[j_halo]==MASTER_RANK){
          x_group[n_group_local] =x_halos[j_halo];
          y_group[n_group_local] =y_halos[j_halo];
          z_group[n_group_local] =z_halos[j_halo];
          vx_group[n_group_local]=vx_halos_sub[j_halo];
          vy_group[n_group_local]=vy_halos_sub[j_halo];
          vz_group[n_group_local]=vz_halos_sub[j_halo];
          n_group_local++;
        }
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);

    // Store group size
    ADaPS_store(&(plist.data),(void *)&n_group,      "n_all_halos",ADaPS_SCALAR_SIZE_T);
    ADaPS_store(&(plist.data),(void *)&n_group_local,"n_halos",    ADaPS_SCALAR_SIZE_T);

    // Needed to make writing more straight-forward
    n_data_ll  =(long long)n_group;
    n_random_ll=(long long)n_random;
    
    if(SID.I_am_Master)
      fprintf(fp,"%lld",n_group);
    for(i_compute=0;i_compute<4;i_compute++){
      switch(i_compute){
      case 0:
        SID_log("Processing real-space ...",SID_LOG_OPEN|SID_LOG_TIMER);
        cfunc_mode=CFUNC_DEFAULT;
        break;
      case 1:
        SID_log("Processing v_x redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
        cfunc_mode=CFUNC_ADD_VX;
        break;
      case 2:
        SID_log("Processing v_y redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
        cfunc_mode=CFUNC_ADD_VY;
        break;
      case 3:
        SID_log("Processing v_z redsift space...",SID_LOG_OPEN|SID_LOG_TIMER);
        cfunc_mode=CFUNC_ADD_VZ;
        break;
      }
      
      // Compute correlation function
      compute_cor_func(&plist,
                       cfunc_mode,
                       cosmo,
                       "halos",
                       "random",
                       redshift,
                       box_size,
                       r_min_l1D,
                       r_max_1D,
                       r_min_2D,
                       r_max_2D,
                       n_1D,
                       n_2D,
                       n_jack,
                       CFUNC_l1D,
                       dCFUNC_l1D,
                       COVMTX_l1D,
                       CFUNC_1D,
                       dCFUNC_1D,
                       COVMTX_1D,
                       CFUNC_2D,
                       dCFUNC_2D,
                       COVMTX_2D,
                       &flag_compute_RR,
                       DD_l1D,
                       DR_l1D,
                       RR_l1D,
                       DD_1D,
                       DR_1D,
                       RR_1D,
                       DD_2D,
                       DR_2D,
                       RR_2D);
      
      // Write correlation function
      if(SID.I_am_Master){
        // Write 1D correlation function
        if(flag_write_header){
          fwrite(&n_groups, sizeof(int),   1,fp_1D);
          fwrite(&n_1D,     sizeof(int),   1,fp_1D);
          fwrite(&n_jack,   sizeof(int),   1,fp_1D);
          fwrite(&r_min_l1D,sizeof(double),1,fp_1D);
          fwrite(&r_max_1D, sizeof(double),1,fp_1D);
        }
        fwrite(&n_data_ll,  sizeof(long long),   1,fp_1D);
        fwrite(&n_random_ll,sizeof(long long),   1,fp_1D);
        fwrite(CFUNC_l1D,   sizeof(double),   n_1D,fp_1D);
        fwrite(dCFUNC_l1D,  sizeof(double),   n_1D,fp_1D);
        fwrite(COVMTX_l1D,  sizeof(double),   n_1D*n_1D,fp_1D);
        for(i_jack=0;i_jack<=n_jack_total;i_jack++){
          fwrite(DD_l1D[i_jack],sizeof(long long),n_1D,fp_1D);
          fwrite(DR_l1D[i_jack],sizeof(long long),n_1D,fp_1D);
          fwrite(RR_l1D[i_jack],sizeof(long long),n_1D,fp_1D);
        }
        fwrite(CFUNC_1D,    sizeof(double),   n_1D,fp_1D);
        fwrite(dCFUNC_1D,   sizeof(double),   n_1D,fp_1D);
        fwrite(COVMTX_1D,   sizeof(double),   n_1D*n_1D,fp_1D);
        for(i_jack=0;i_jack<=n_jack_total;i_jack++){
          fwrite(DD_1D[i_jack],sizeof(long long),n_1D,fp_1D);
          fwrite(DR_1D[i_jack],sizeof(long long),n_1D,fp_1D);
          fwrite(RR_1D[i_jack],sizeof(long long),n_1D,fp_1D);
        }

        // Write 2D correlation function
        if(flag_write_header){
          fwrite(&n_groups,sizeof(int),   1,fp_2D);
          fwrite(&n_2D,    sizeof(int),   1,fp_2D);
          fwrite(&n_jack,   sizeof(int),  1,fp_2D);
          fwrite(&r_min_2D,sizeof(double),1,fp_2D);
          fwrite(&r_max_2D,sizeof(double),1,fp_2D);
        }
        fwrite(&n_data_ll,  sizeof(long long),        1,fp_2D);
        fwrite(&n_random_ll,sizeof(long long),        1,fp_2D);
        fwrite(CFUNC_2D,    sizeof(double),   n_2D*n_2D,fp_2D);
        fwrite(dCFUNC_2D,   sizeof(double),   n_2D*n_2D,fp_2D);
        fwrite(COVMTX_2D,   sizeof(double),   n_2D_total*n_2D_total,fp_2D);
        for(i_jack=0;i_jack<=n_jack_total;i_jack++){
          fwrite(DD_2D[i_jack],sizeof(long long),n_2D*n_2D,fp_2D);
          fwrite(DR_2D[i_jack],sizeof(long long),n_2D*n_2D,fp_2D);
          fwrite(RR_2D[i_jack],sizeof(long long),n_2D*n_2D,fp_2D);
        }

        flag_write_header=FALSE;
    
        // Compute (and print) r_o and r_X (the cross-over radius)
        r_o=0.;
        r_X=0.;
        x_interp=(double *)SID_malloc(sizeof(double)*n_1D);
        y_interp=(double *)SID_malloc(sizeof(double)*n_1D);
        i_bin=0;
        while(CFUNC_l1D[i_bin+1]>=CFUNC_l1D[i_bin] || DD_l1D[0][i_bin]<10){
          i_bin++;
          if(i_bin>=(n_1D-1)) break;
        }
        x_interp[0]=CFUNC_l1D[i_bin];
        y_interp[0]=take_log10(r_min_l1D)+((double)i_bin+0.5)*dlr_l1D;
        for(j_bin=1;i_bin<n_1D;i_bin++){
          if(CFUNC_l1D[i_bin]<x_interp[j_bin-1]){
            x_interp[j_bin]=CFUNC_l1D[i_bin];
            y_interp[j_bin]=take_log10(r_min_l1D)+((double)i_bin+0.5)*dlr_l1D;
            j_bin++;
          }          
        }

        if(x_interp[0]*x_interp[j_bin-1]<0.){
          init_interpolate(x_interp,y_interp,(size_t)j_bin,gsl_interp_linear,&r_o_interp);
          r_X=take_alog10(interpolate(r_o_interp,0.));
          free_interpolate(&r_o_interp);
        }
        else if(x_interp[0]>0.)
          r_X=999.999;
        else
          r_X=0.;
        if(x_interp[0]>1. && x_interp[j_bin-1]<1.){
          init_interpolate(x_interp,y_interp,(size_t)j_bin,gsl_interp_linear,&r_o_interp);
          r_o=take_alog10(interpolate(r_o_interp,1.));
          free_interpolate(&r_o_interp);
        }
        else if(x_interp[0]<=1.)
          r_o=0.;
        else
          r_o=999.999;
        SID_free(SID_FARG x_interp);
        SID_free(SID_FARG y_interp);          
        SID_log("r_o=%7.3lf [Mpc/h]",SID_LOG_COMMENT,r_o);
        SID_log("r_X=%7.3lf [Mpc/h]",SID_LOG_COMMENT,r_X);
        fprintf(fp," %le %le",r_o,r_X);
      }
      SID_log("Done.",SID_LOG_CLOSE);
    }
    if(SID.I_am_Master)
      fprintf(fp,"\n");

    SID_log("Done.",SID_LOG_CLOSE);
  }
  if(SID.I_am_Master){
    fclose(fp);
    fclose(fp_1D);
    fclose(fp_2D);
  }
  
  // Clean-up
  SID_log("Cleaning-up...",SID_LOG_OPEN);
  for(i_jack=0;i_jack<=n_jack_total;i_jack++){
    SID_free(SID_FARG DD_l1D[i_jack]);
    SID_free(SID_FARG DR_l1D[i_jack]);
    SID_free(SID_FARG RR_l1D[i_jack]);
    SID_free(SID_FARG DD_1D[i_jack]);
    SID_free(SID_FARG DR_1D[i_jack]);
    SID_free(SID_FARG RR_1D[i_jack]);
    SID_free(SID_FARG DD_2D[i_jack]);
    SID_free(SID_FARG DR_2D[i_jack]);
    SID_free(SID_FARG RR_2D[i_jack]);
  }
  SID_free(SID_FARG DD_1D);
  SID_free(SID_FARG DR_1D);
  SID_free(SID_FARG RR_1D);
  SID_free(SID_FARG DD_2D);
  SID_free(SID_FARG DR_2D);
  SID_free(SID_FARG RR_2D);
  SID_free(SID_FARG CFUNC_l1D);
  SID_free(SID_FARG dCFUNC_l1D);
  SID_free(SID_FARG COVMTX_l1D);
  SID_free(SID_FARG CFUNC_1D);
  SID_free(SID_FARG dCFUNC_1D);
  SID_free(SID_FARG COVMTX_1D);
  SID_free(SID_FARG CFUNC_2D);
  SID_free(SID_FARG dCFUNC_2D);
  SID_free(SID_FARG COVMTX_2D);
  free_plist(&plist);
  free_RNG(&RNG);
  if(SID.I_am_Master)
    SID_free(SID_FARG rank_own);
  SID_free(SID_FARG x_halos);
  SID_free(SID_FARG y_halos);
  SID_free(SID_FARG z_halos);
  SID_free(SID_FARG vx_halos_sub);
  SID_free(SID_FARG vy_halos_sub);
  SID_free(SID_FARG vz_halos_sub);
  SID_free(SID_FARG vx_halos_FoF);
  SID_free(SID_FARG vy_halos_FoF);
  SID_free(SID_FARG vz_halos_FoF);
  SID_free(SID_FARG vx_halos_sys);
  SID_free(SID_FARG vy_halos_sys);
  SID_free(SID_FARG vz_halos_sys);
  
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
