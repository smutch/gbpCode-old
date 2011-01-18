#define  _MAIN
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpCosmo.h>
#include <gbpClustering.h>

#define ALLOC_FACTOR_LOCAL 1.5

void read_gadget_binary_local(char       *filename_root_in,
                              int         snapshot_number,
                              plist_info *plist);
void read_gadget_binary_local(char       *filename_root_in,
                              int         snapshot_number,
                              plist_info *plist){
  char    **pname;
  char     *name_initpositions;
  char      filename[MAX_FILENAME_LENGTH];
  char     *read_catalog;
  size_t    i,j,k,l,jj;
  size_t    k_offset[N_GADGET_TYPE];
  size_t    n_particles_file;
  size_t    n_particles_kept;
  size_t    n_particles_kept_all;
  size_t    n_of_type_rank[N_GADGET_TYPE];
  size_t    n_allocate[N_GADGET_TYPE];
  size_t    n_of_type[N_GADGET_TYPE];
  size_t    i_particle;
  size_t    n_unallocated_all;
  size_t    n_unallocated[N_GADGET_TYPE];
  int       n_type_used;
  int       n_particles_all_in_groups;
  int       n_particles_in_groups;
  size_t    n_particles_rank;
  size_t    n_particles_all;
  size_t    i_p,i_p_temp;
  size_t    id_min,id_max;
  double    mass_array[N_GADGET_TYPE];
  int       unused[GADGET_HEADER_SIZE];
  size_t    n_all[N_GADGET_TYPE];
  int       n_all_tmp[N_GADGET_TYPE];
  int       flag_metals;
  int       flag_ages;
  int       flag_entropyICs;
  double    expansion_factor;
  double    h_Hubble;
  double    box_size;
  double    d_value;
  double    d1_value;
  double    d2_value;
  double    d3_value;
  float     f_value;
  float     f1_value;
  float     f2_value;
  float     f3_value;
  size_t    id_test;
  int       i_value;
  int       i1_value;
  int       i2_value;
  int       i3_value;
  int       n_warning;
  REAL     *x_array[N_GADGET_TYPE];
  REAL     *y_array[N_GADGET_TYPE];
  REAL     *z_array[N_GADGET_TYPE];
  REAL     *vx_array[N_GADGET_TYPE];
  REAL     *vy_array[N_GADGET_TYPE];
  REAL     *vz_array[N_GADGET_TYPE];
  size_t   *id_array[N_GADGET_TYPE];
  size_t   *id_list;
  size_t   *id_list_offset;
  size_t   *id_list_in;
  size_t   *id_list_index;
  size_t    n_id_list;
  int       record_length_open;
  int       record_length_close;
  int       n_return;
  size_t    s_load;
  char     *keep;
  int       n_keep[N_GADGET_TYPE];
  int       flag_keep_IDs      =TRUE;
  int       flag_multifile     =FALSE;
  int       flag_multimass     =FALSE;
  int       flag_file_type     =0;
  int       flag_filefound     =FALSE;
  int       flag_gas           =FALSE;
  int       flag_initpositions =FALSE;
  int       flag_no_velocities =FALSE;
  int       flag_LONGIDs       =FALSE;
  int       flag_read_marked   =FALSE;
  int       flag_read_catalog  =FALSE;
  int       flag_all_read_all  =FALSE;
  int       i_file;
  int       i_rank;
  int       n_files;
  int       n_files_in;
  double    x_scale;
  double    y_scale;
  double    z_scale;
  double    x_offset;
  double    y_offset;
  double    z_offset;
  double    x_period;
  double    y_period;
  double    z_period;
  int       flag_xscale;
  int       flag_yscale;
  int       flag_zscale;
  int       flag_xoffset;
  int       flag_yoffset;
  int       flag_zoffset;
  int       flag_xperiod;
  int       flag_yperiod;
  int       flag_zperiod;
  FILE     *fp;
  void     *buffer;
  size_t   *buffer_index;
  void     *buffer_i;
  size_t    id_search;
  size_t    id_value;
  size_t    n_buffer;
  int       n_non_zero;
  int       seed =1073743;
  int       seed2=1073743;
  int       scatter_rank;
  int       read_mode;
  int       mark_mode;
  double    x_min_read;
  double    x_max_read;
  double    y_min_read;
  double    y_max_read;
  double    z_min_read;
  double    z_max_read;
  double    x_min_bcast;
  double    x_max_bcast;
  double    y_min_bcast;
  double    y_max_bcast;
  double    z_min_bcast;
  double    z_max_bcast;
  gadget_header_info header;
  int                read_rank;
 
  // Determine file format
  n_files=1;
  for(i_file=0;i_file<3 && !flag_filefound;i_file++){  
    if(i_file==0)
      sprintf(filename,"%s/snapshot_%03d/snapshot_%03d",filename_root_in,snapshot_number,snapshot_number);
    else if(i_file==1)
      sprintf(filename,"%s/snapshot_%03d",filename_root_in,snapshot_number);
    else if(i_file==2)
      sprintf(filename,"%s",filename_root_in,snapshot_number);
    else if(i_file==3)
      sprintf(filename,"%s_%03d",filename_root_in,snapshot_number);
    fp=fopen(filename,"r");
    if(fp!=NULL){
      flag_filefound=TRUE;
      flag_multifile=FALSE;
      flag_file_type=i_file;
    }
    // ... if that doesn't work, check for multi-file
    else{
      strcat(filename,".0");
      fp=fopen(filename,"r");
      if(fp!=NULL){
        flag_filefound=TRUE;
        flag_multifile=TRUE;
        flag_file_type=i_file;
      }
    }
  }

  // A file was found ... 
  if(flag_filefound){
    SID_log("Reading GADGET binary file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_root_in);

    pname=plist->species;
    x_min_read=((double *)ADaPS_fetch(plist->data,"rank_x_min"))[0];
    x_max_read=((double *)ADaPS_fetch(plist->data,"rank_x_max"))[0];

    // Header record length
    SID_log("Reading header...",SID_LOG_OPEN);
    fread(&record_length_open,4,1,fp);
    if(record_length_open!=GADGET_HEADER_SIZE)
      SID_log_warning("Problem with GADGET record size (opening size of header is wrong)",ERROR_LOGIC);

    // Number of particles for each species in this file
    fread(&header,sizeof(gadget_header_info),1,fp); 
    for(i=0;i<N_GADGET_TYPE;i++)
      n_of_type[i]=(size_t)header.n_file[i];

    // Expansion factor (or time) 
    SID_log("a=%le...",SID_LOG_CONTINUE,expansion_factor);
    expansion_factor=header.time;
    ADaPS_store(&(plist->data),(void *)(&(header.time)),"expansion_factor",ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&(plist->data),(void *)(&(header.time)),"time",            ADaPS_SCALAR_DOUBLE);

    // Redshift
    d_value=(double)header.redshift;
    ADaPS_store(&(plist->data),(void *)(&d_value),"redshift",ADaPS_SCALAR_DOUBLE);

    // Number of particles for each species in all files
    for(i=0;i<N_GADGET_TYPE;i++)
      n_all[i]=(size_t)header.n_all[i];

    // Number of files in this snapshot 
    ADaPS_store(&(plist->data),(void *)(&(header.n_files)),"n_files",ADaPS_SCALAR_INT);
    n_files=header.n_files;

    // Cosmology 

    // Omega_o
    d_value=(double)header.Omega_M; 
    ADaPS_store(&(plist->data),(void *)(&d_value),"Omega_M",ADaPS_SCALAR_DOUBLE);

    // Omega_Lambda 
    d_value=(double)header.Omega_Lambda; 
    ADaPS_store(&(plist->data),(void *)(&d_value),"Omega_Lambda",ADaPS_SCALAR_DOUBLE);

    // Hubble parameter 
    box_size=header.box_size;
    ADaPS_store(&(plist->data),(void *)(&box_size),"box_size",ADaPS_SCALAR_DOUBLE);

    for(i=0;i<N_GADGET_TYPE;i++){
      if(header.n_all_high_word[i]>0)
        n_all[i]+=(((size_t)header.n_all_high_word[i]) << 32);
    }

    // Count the total number of particles
    n_particles_all=0;
    for(i=0,n_non_zero=0;i<N_GADGET_TYPE;i++){
      if(n_all[i]>0){
        n_particles_all+=n_all[i];
        n_non_zero++;
      }
    }

    // List numbers of particles in the log output
    SID_log("%lld",SID_LOG_CONTINUE,n_particles_all);
    if(n_non_zero>0)
      SID_log(" (",SID_LOG_CONTINUE,n_particles_all);
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_all[i]>0){
        if(i==n_non_zero-1){
          if(n_non_zero>1)
            SID_log("and %lld %s",SID_LOG_CONTINUE,n_all[i],pname[i]);
          else
            SID_log("%lld %s",SID_LOG_CONTINUE,n_all[i],pname[i]);
        }
        else{
          if(n_non_zero>1)
            SID_log("%lld %s, ",SID_LOG_CONTINUE,n_all[i],pname[i]);        
          else
            SID_log("%lld %s",SID_LOG_CONTINUE,n_all[i],pname[i]);        
        }
      }
    }
    if(n_non_zero>0)
      SID_log(") particles...",SID_LOG_CONTINUE);
    else
      SID_log(" particles...",SID_LOG_CONTINUE);

    // Check closing record length
    fread(&record_length_close,4,1,fp);
    if(record_length_open!=record_length_close)
      SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);

    // Close file
    fclose(fp);
    SID_log("Done.",SID_LOG_CLOSE);

    // Store mass array
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_all[i]>0){
        if(mass_array[i]!=0.){
          d_value=mass_array[i];
          ADaPS_store(&(plist->data),(void *)(&d_value),"mass_array_%s",ADaPS_SCALAR_DOUBLE,plist->species[i]);
        }
      }
      else
        mass_array[i]=0.;
    }

    // Is this a multimass snapshot?
    flag_multimass=FALSE;
    for(i=0;i<N_GADGET_TYPE;i++)
      if(n_all[i]>0 && mass_array[i]==0.)
        flag_multimass=TRUE;

    // Is sph being used?
    if(n_all[GADGET_TYPE_GAS]>0)
      flag_gas=TRUE;
    else
      flag_gas=FALSE;

    if(SID.n_proc>1){
      for(i=0;i<N_GADGET_TYPE;i++){
        n_of_type_rank[i]=0;
        n_allocate[i]     =ALLOC_FACTOR_LOCAL*(size_t)((float)n_all[i]/(float)SID.n_proc);
        n_unallocated[i]  =0;
      }
    }
    else{
      for(i=0;i<N_GADGET_TYPE;i++){
        n_of_type_rank[i]=0;
        n_allocate[i]     =(size_t)n_all[i];
        n_unallocated[i]  =0;
      }
    }

    // Allocate data arrays
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_allocate[i]>0){
        x_array[i] =(REAL   *)SID_malloc(sizeof(REAL)  *(size_t)n_allocate[i]);
        y_array[i] =(REAL   *)SID_malloc(sizeof(REAL)  *(size_t)n_allocate[i]);
        z_array[i] =(REAL   *)SID_malloc(sizeof(REAL)  *(size_t)n_allocate[i]);
        vx_array[i]=(REAL   *)SID_malloc(sizeof(REAL)  *(size_t)n_allocate[i]);
        vy_array[i]=(REAL   *)SID_malloc(sizeof(REAL)  *(size_t)n_allocate[i]);
        vz_array[i]=(REAL   *)SID_malloc(sizeof(REAL)  *(size_t)n_allocate[i]);
      }
    }
    
    // Set offsets
    for(i=0;i<N_GADGET_TYPE;i++)
      k_offset[i]=0;

    // Read data
    for(i_file=0,n_particles_kept=0;i_file<n_files;i_file++){
      read_rank=i_file%SID.n_proc;
      read_rank=0;
      if(n_files>1)
        SID_log("Reading file %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file+1,n_files);

      if(flag_file_type==0)
        sprintf(filename,"%s/snapshot_%03d/snapshot_%03d",filename_root_in,snapshot_number,snapshot_number);
      else if(flag_file_type==1)
        sprintf(filename,"%s/snapshot_%03d",filename_root_in,snapshot_number);
      else if(flag_file_type==2)
        sprintf(filename,"%s",filename_root_in,snapshot_number);
      else if(flag_file_type==3)
        sprintf(filename,"%s_%03d",filename_root_in,snapshot_number);
      if(flag_multifile)
        sprintf(filename,"%s.%d",filename,i_file);
      // Open file and read header
      if(SID.My_rank==read_rank){
        fp=fopen(filename,"r");
        fread(&record_length_open,4,1,fp);
        fread(&header,sizeof(gadget_header_info),1,fp);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
      }
      SID_Bcast(&header,(int)sizeof(gadget_header_info),read_rank);
      for(i=0,n_particles_file=0;i<N_GADGET_TYPE;i++)
        n_particles_file+=(size_t)header.n_file[i];

      // Initialize buffer
      if(SID.My_rank==read_rank)
        fread(&record_length_open,4,1,fp);
      SID_Bcast(&record_length_open,(int)sizeof(int),read_rank);
      buffer=SID_malloc((size_t)record_length_open);
      keep  =(char *)SID_malloc(sizeof(char)*n_particles_file);

      // Read positions
      SID_log("Reading positions...",SID_LOG_OPEN|SID_LOG_TIMER);
      if(SID.My_rank==read_rank){
        fread(buffer,(size_t)record_length_open,1,fp);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of positions)",ERROR_LOGIC);
      }
      SID_Barrier(SID.COMM_WORLD);
      SID_Bcast(buffer,(int)record_length_open,read_rank);
      for(i=0,jj=0;i<N_GADGET_TYPE;i++){
        n_keep[i]=0;
        for(j=0,k=0;j<header.n_file[i];j++,jj++){
          f1_value=((float *)buffer)[3*jj+0];
          f2_value=((float *)buffer)[3*jj+1];
          f3_value=((float *)buffer)[3*jj+2];
          d1_value=((double)(f1_value));
          if(d1_value>=x_min_read && d1_value<x_max_read)
            keep[jj]=TRUE;
          if(keep[jj]){
            if(k+k_offset[i]>=n_allocate[i])
              n_unallocated[i]++;
            else{
              x_array[i][k+k_offset[i]]=((REAL)(f1_value))*(REAL)(plist->length_unit/M_PER_MPC);
              y_array[i][k+k_offset[i]]=((REAL)(f2_value))*(REAL)(plist->length_unit/M_PER_MPC);
              z_array[i][k+k_offset[i]]=((REAL)(f3_value))*(REAL)(plist->length_unit/M_PER_MPC);
              k++;
              n_keep[i]++;
              n_particles_kept++;
            }
          }
        }
        if(k!=n_keep[i])
          SID_trap_error("Particle count mismatch (ie. %d!=%d) during positions read",ERROR_LOGIC,k,n_keep[i]);
      }
      for(i=0;i<N_GADGET_TYPE;i++){
        if(n_unallocated[i]>0)
          SID_log_warning("Rank %d has exceeded it's memory allocation (%lld) for particle type %d",ERROR_LOGIC,SID.My_rank,n_allocate[i],i);
      }
      n_unallocated_all=calc_sum_global(n_unallocated,N_GADGET_TYPE,SID_SIZE_T);
      if(n_unallocated_all>0)
        SID_trap_error("You need to increase ALLOC_FACTOR",ERROR_LOGIC);
      SID_log("Done.",SID_LOG_CLOSE);

      // Initialize buffer
      if(SID.My_rank==read_rank)
        fread(&record_length_open,4,1,fp);
      SID_Bcast(&record_length_open,(int)sizeof(int),read_rank);

      // Read velocities
      SID_log("Reading velocities...",SID_LOG_OPEN|SID_LOG_TIMER);
      if(SID.My_rank==read_rank){
        fread(buffer,(size_t)record_length_open,1,fp);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of positions)",ERROR_LOGIC);
      }
      SID_Barrier(SID.COMM_WORLD);
      SID_Bcast(buffer,(int)record_length_open,read_rank);
      for(i=0,jj=0;i<N_GADGET_TYPE;i++){
        for(j=0,k=0;j<header.n_file[i];j++,jj++){
          f1_value=((float *)buffer)[3*jj+0];
          f2_value=((float *)buffer)[3*jj+1];
          f3_value=((float *)buffer)[3*jj+2];
          if(keep[jj]){
            if(k+k_offset[i]>=n_allocate[i])
              n_unallocated[i]++;
            else{
              vx_array[i][k+k_offset[i]]=((REAL)(f1_value))*(REAL)(plist->velocity_unit*sqrt(expansion_factor)*1e-3);
              vy_array[i][k+k_offset[i]]=((REAL)(f2_value))*(REAL)(plist->velocity_unit*sqrt(expansion_factor)*1e-3);
              vz_array[i][k+k_offset[i]]=((REAL)(f3_value))*(REAL)(plist->velocity_unit*sqrt(expansion_factor)*1e-3);
              k++;
            }
          }
        }
        if(k!=n_keep[i])
          SID_trap_error("Particle count mismatch (ie. %d!=%d) during velocities read",ERROR_LOGIC,k,n_keep[i]);
      }
      for(i=0;i<N_GADGET_TYPE;i++){
        if(n_unallocated[i]>0)
          SID_log_warning("Rank %d has exceeded it's memory allocation (%lld) for particle type %d",ERROR_LOGIC,SID.My_rank,n_allocate[i],i);
      }
      n_unallocated_all=calc_sum_global(n_unallocated,N_GADGET_TYPE,SID_SIZE_T);
      if(n_unallocated_all>0)
        SID_trap_error("You need to increase ALLOC_FACTOR",ERROR_LOGIC);
      SID_log("Done.",SID_LOG_CLOSE);
        
      // Update offsets
      for(i=0;i<N_GADGET_TYPE;i++)
        k_offset[i]+=n_keep[i];
      
      // Clean-up
      if(SID.My_rank==read_rank)
        fclose(fp);
      SID_free(SID_FARG buffer);
      SID_free(SID_FARG keep);
      if(n_files>1)
        SID_log("Done.",SID_LOG_CLOSE);
    }
    
    for(i=0,n_particles_rank=0;i<N_GADGET_TYPE;i++){
      n_of_type_rank[i]=k_offset[i];
      n_particles_rank+=k_offset[i];
    }
    
    // Store everything in the data structure...

    //   ... particle counts ...
    for(i=0,n_particles_kept_all=0;i<N_GADGET_TYPE;i++){
#ifdef USE_MPI
      MPI_Allreduce(&(n_of_type_rank[i]),&(n_of_type[i]),1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
#else
      n_of_type[i]=n_of_type_rank[i];
#endif
      if(n_of_type[i]>0){
        n_particles_kept_all+=n_of_type[i];
        ADaPS_store(&(plist->data),(void *)(&(n_of_type_rank[i])),"n_%s",    ADaPS_SCALAR_SIZE_T,pname[i]);
        ADaPS_store(&(plist->data),(void *)(&(n_of_type[i])),     "n_all_%s",ADaPS_SCALAR_SIZE_T,pname[i]);
      }
    }
    ADaPS_store(&(plist->data),(void *)(&n_particles_all),"n_particles_all",ADaPS_SCALAR_SIZE_T);

    // Check that the right number of particles have been read
    if(n_particles_kept!=n_particles_rank)
      SID_log_warning("Rank %d did not receive the right number of particles (ie. %lld!=%lld)",ERROR_LOGIC,SID.My_rank,n_particles_kept,n_particles_rank);
    if(n_particles_kept_all!=n_particles_all)
      SID_log_warning("The right number of particles were not read (ie. %lld!=%lld)",ERROR_LOGIC,n_particles_kept_all,n_particles_all);

    //   ... positions ...
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type_rank[i]>0){
        ADaPS_store(&(plist->data),(void *)x_array[i],"x_%s",ADaPS_DEFAULT,pname[i]);
        ADaPS_store(&(plist->data),(void *)y_array[i],"y_%s",ADaPS_DEFAULT,pname[i]);
        ADaPS_store(&(plist->data),(void *)z_array[i],"z_%s",ADaPS_DEFAULT,pname[i]);
      }
    }
    //  ... velocities ...
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type_rank[i]>0){
        ADaPS_store(&(plist->data),(void *)vx_array[i],"vx_%s",ADaPS_DEFAULT,pname[i]);
        ADaPS_store(&(plist->data),(void *)vy_array[i],"vy_%s",ADaPS_DEFAULT,pname[i]);
        ADaPS_store(&(plist->data),(void *)vz_array[i],"vz_%s",ADaPS_DEFAULT,pname[i]);
      }
    }

    SID_log("Done.",SID_LOG_CLOSE);
  }
  else
    SID_trap_error("Could not find snapshot %d for file with root {%s}",ERROR_IO_OPEN,snapshot_number,filename_root_in);
}


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
  int     i_species;
  char    filename_root[MAX_FILENAME_LENGTH];
  char    filename_out[MAX_FILENAME_LENGTH];
  char    filename_out_1D[MAX_FILENAME_LENGTH];
  char    filename_out_2D[MAX_FILENAME_LENGTH];
  char    filename_out_root[MAX_FILENAME_LENGTH];
  char    filename_TF[256];
  char    n_string[64];
  int             n[3];
  double          L[3];
  size_t          n_all;
  FILE           *fp_1D;
  FILE           *fp_2D;
  cosmo_info     *cosmo;
  field_info      FFT;
  plist_info      plist_header;
  plist_info      plist;
  FILE           *fp;
  int     i_temp;
  int     n_temp;
  double *k_temp;
  double *kmin_temp;
  double *kmax_temp;
  double *P_temp;
  size_t *n_mode_temp;
  double *sigma_P_temp;
  double *shot_noise_temp;
  double *dP_temp;
  int     pspec_mode;
  int     snapshot_number;
  int     i_compute;
  int     mass_assignment_scheme;
  double  k_min_1D;
  double  k_max_1D;
  double  k_min_2D;
  double  k_max_2D;
  int     n_k_1D;
  int     n_k_2D;
  double *k_1D;
  double *P_k_1D;
  double *dP_k_1D;
  int    *n_modes_1D;
  double *P_k_2D;
  double *dP_k_2D;
  int    *n_modes_2D;
  int     n_groups=1;
  double  dk_1D;
  double  dk_2D;

  // Initialization -- MPI etc.
  SID_init(&argc,&argv,NULL);

  // Initialization -- user input
  strcpy(filename_root,argv[1]);
  snapshot_number=atoi(argv[2]);
  L[0]           =(double)atof(argv[3]);
  n[0]           =atoi(argv[4]);
  k_min_1D       =(double)atof(argv[5]);
  k_max_1D       =(double)atof(argv[6]);
  n_k_1D         =atoi(argv[7]);
  k_min_2D       =0.;
  k_max_2D       =(double)atof(argv[8]);
  n_k_2D         =atoi(argv[9]);
  strcpy(filename_out_root,argv[10]);

  dk_1D=(k_max_1D-k_min_1D)/(double)n_k_1D;
  dk_2D=(k_max_2D-k_min_2D)/(double)n_k_2D;

  SID_log("Processing the power spectra of %s, snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,filename_root,snapshot_number);

  n[1]=n[0];
  n[2]=n[0];
  L[1]=L[0];
  L[2]=L[0];

  // Initialization -- default cosmology
  init_cosmo_std(&cosmo);

  // Initialization -- FFT
  init_field(3,n,L,&FFT);

  // Initialization -- data structure which holds all   
  //                   the (local) particle information  
  init_plist(&plist,&(FFT.slab),GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);

  // Initialization -- read particle file and fetch header info
  read_gadget_binary_local(filename_root,snapshot_number,&plist);

  // Initialize arrays
  k_1D      =(double *)SID_malloc(sizeof(double)*n_k_1D);
  P_k_1D    =(double *)SID_malloc(sizeof(double)*n_k_1D);
  dP_k_1D   =(double *)SID_malloc(sizeof(double)*n_k_1D);
  n_modes_1D=(int    *)SID_malloc(sizeof(int)   *n_k_1D);
  P_k_2D    =(double *)SID_malloc(sizeof(double)*n_k_2D*n_k_2D);
  dP_k_2D   =(double *)SID_malloc(sizeof(double)*n_k_2D*n_k_2D);
  n_modes_2D=(int    *)SID_malloc(sizeof(int)   *n_k_2D*n_k_2D);    

  // Process each species in turn
  mass_assignment_scheme=MAP2GRID_DIST_DWT20;
  for(i_species=0;i_species<plist.n_species;i_species++){
    sprintf(n_string,"n_all_%s",plist.species[i_species]);
    if(ADaPS_exist(plist.data,n_string)){
      n_all=((size_t *)ADaPS_fetch(plist.data,n_string))[0];
      if(n_all>0){
        SID_log("Processing the power spectra of %lld %s particles...",SID_LOG_OPEN|SID_LOG_TIMER,n_all,plist.species[i_species]);
        if(SID.I_am_Master){
          sprintf(filename_out_1D,   "%s_%s.1D_pow_spec",filename_out_root,plist.species[i_species]);
          sprintf(filename_out_2D,   "%s_%s.2D_pow_spec",filename_out_root,plist.species[i_species]);
          SID_log("One-D filename is {%s}",SID_LOG_COMMENT,filename_out_1D);
          SID_log("Two-D filename is {%s}",SID_LOG_COMMENT,filename_out_2D);
          fp_1D=fopen(filename_out_1D,"wb");
          fp_2D=fopen(filename_out_2D,"wb");
        }
        for(i_compute=0;i_compute<4;i_compute++){
          // Compute power spectrum
          switch(i_compute){
          case 0:
            SID_log("Processing real-space ...",SID_LOG_OPEN|SID_LOG_TIMER);
            pspec_mode=PSPEC_DEFAULT;
            break;
          case 1:
            SID_log("Processing v_x redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
            pspec_mode=PSPEC_ADD_VX;
            break;
          case 2:
            SID_log("Processing v_y redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
            pspec_mode=PSPEC_ADD_VY;
            break;
          case 3:
            SID_log("Processing v_z redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
            pspec_mode=PSPEC_ADD_VZ;
            break;
          }
          compute_power_spectrum(&plist,
                                 &FFT,
                                 mass_assignment_scheme,
                                 pspec_mode,
                                 cosmo,
                                 plist.species[i_species],
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
          // Write results
          if(SID.I_am_Master){
            // Write 1D power spectra
            if(i_compute==0){
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
            if(i_compute==0){
              fwrite(&n_groups,sizeof(int),   1,fp_2D);
              fwrite(&n_k_2D,  sizeof(int),   1,fp_2D);          
              fwrite(&k_min_2D,sizeof(double),1,fp_2D);
              fwrite(&k_max_2D,sizeof(double),1,fp_2D);
              fwrite(&dk_2D,   sizeof(double),1,fp_2D);
            }
            fwrite(P_k_2D,    sizeof(double),n_k_2D*n_k_2D,fp_2D);
            fwrite(dP_k_2D,   sizeof(double),n_k_2D*n_k_2D,fp_2D);
            fwrite(n_modes_2D,sizeof(int),   n_k_2D*n_k_2D,fp_2D);
          }
          SID_log("Done.",SID_LOG_CLOSE);
        }
        // Close files
        if(SID.I_am_Master){
          fclose(fp_1D);
          fclose(fp_2D);
        }  
        SID_log("Done.",SID_LOG_CLOSE);
      }
    }
  }

    
  // Clean-up
  SID_log("Clean-up...",SID_LOG_OPEN);
  free_cosmo(&cosmo);
  free_plist(&plist);
  SID_free(SID_FARG k_1D);
  SID_free(SID_FARG P_k_1D);
  SID_free(SID_FARG dP_k_1D);
  SID_free(SID_FARG n_modes_1D);
  SID_free(SID_FARG P_k_2D);
  SID_free(SID_FARG dP_k_2D);
  SID_free(SID_FARG n_modes_2D);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}
