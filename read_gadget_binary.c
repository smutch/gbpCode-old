#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>

#define READ_GADGET_NONE         2
#define READ_GADGET_DEFAULT      4
#define READ_GADGET_RANDOM       8
#define READ_GADGET_MARKED      16
#define READ_GADGET_X_MIN       32
#define READ_GADGET_X_MAX       64
#define READ_GADGET_Y_MIN      128
#define READ_GADGET_Y_MAX      256
#define READ_GADGET_Z_MIN      512
#define READ_GADGET_Z_MAX     1024

int assign_particle_to_rank(plist_info *plist,
                            double      x_p,
                            double      y_p,
                            double      z_p,
                            size_t     *id_list,
                            size_t     *id_list_index,
                            size_t      n_id_list,
                            size_t      id_p,
                            size_t      i_p,
                            RNG_info   *RNG,
                            int         mode,
                            int         mode_mark,
                            double      x_min,
                            double      x_max,
                            double      y_min,
                            double      y_max,
                            double      z_min,
                            double      z_max,
                            size_t      id_min,
                            size_t      id_max);
int assign_particle_to_rank(plist_info *plist,
                            double      x_p,
                            double      y_p,
                            double      z_p,
                            size_t     *id_list,
                            size_t     *id_list_index,
                            size_t      n_id_list,
                            size_t      id_p,
                            size_t      i_p,
                            RNG_info   *RNG,
                            int         mode,
                            int         mode_mark,
                            double      x_min,
                            double      x_max,
                            double      y_min,
                            double      y_max,
                            double      z_min,
                            double      z_max,
                            size_t      id_min,
                            size_t      id_max){
  int      keep=TRUE;
  int      scatter_rank;

  if(mode_mark!=READ_GADGET_NONE){
    if(check_mode_for_flag(mode_mark,READ_GADGET_MARKED)){
      if(id_list[id_list_index[find_index(id_list,id_p,n_id_list,id_list_index)]]!=id_p)
        keep*=FALSE;    
    }
#ifdef USE_MPI
    else if(check_mode_for_flag(mode,READ_GADGET_RANDOM)){
      scatter_rank=MIN(SID.n_proc-1,(int)(random_number(RNG)*(REAL)SID.n_proc));
      if(scatter_rank!=SID.My_rank)
        keep*=FALSE;
    }
#endif
    else if(check_mode_for_flag(mode,READ_GADGET_DEFAULT)){
      if(i_p<id_min || i_p>id_max)
        keep*=FALSE;
    }
    if(check_mode_for_flag(mode,READ_GADGET_X_MIN)){
      if(x_p<x_min)
        keep*=FALSE;
    }
    if(check_mode_for_flag(mode,READ_GADGET_X_MAX)){
      if(x_p>=x_max)
        keep*=FALSE;
    }
    if(check_mode_for_flag(mode,READ_GADGET_Y_MIN)){
      if(y_p<y_min)
        keep*=FALSE;
    }
    if(check_mode_for_flag(mode,READ_GADGET_Y_MAX)){
      if(y_p>=y_max)
        keep*=FALSE;
    }
    if(check_mode_for_flag(mode,READ_GADGET_Z_MIN)){
      if(z_p<z_min)
        keep*=FALSE;
    }
    if(check_mode_for_flag(mode,READ_GADGET_Z_MAX)){
      if(z_p>=z_max)
        keep*=FALSE;
    }
  }
  else
    keep*=FALSE;
  return(keep);  
}

void read_gadget_binary(char       *filename_in,
                        plist_info *plist,
                        int         mode){
  char    **pname;
  char     *name_initpositions;
  char      filename[MAX_FILENAME_LENGTH];
  char     *read_catalog;
  size_t    i,j,k,l,jj;
  size_t    k_offset[N_GADGET_TYPE];
  size_t    k_offset_mass[N_GADGET_TYPE];
  size_t    n_of_type[N_GADGET_TYPE];
  size_t    n_of_type_file[N_GADGET_TYPE];
  size_t    n_particles_file;
  unsigned int n_of_type_tmp[N_GADGET_TYPE];
  size_t    n_of_type_rank[N_GADGET_TYPE];
  size_t   *n_of_type_file_rank[N_GADGET_TYPE];
  int       n_type_used;
  size_t   *n_rank;
  size_t  **list_file_rank[N_GADGET_TYPE];
  size_t    n_particles_rank;
  size_t    n_particles_bcast;
  size_t    n_particles_mass_rank;
  size_t    n_particles_mass;
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
  size_t    id_value_LL;
  int       id_value_I;
  int       i_value;
  int       i1_value;
  int       i2_value;
  int       i3_value;
  REAL     *x_array;
  REAL     *y_array;
  REAL     *z_array;
  REAL     *vx_array;
  REAL     *vy_array;
  REAL     *vz_array;
  size_t   *id_array;
  double   *M_array;
  REAL     *u_array;
  REAL     *rho_array;
  REAL     *T_array;
  REAL     *smooth_array;
  size_t   *id_list;
  size_t   *id_list_index;
  size_t    n_id_list;
  int       record_length_open;
  int       record_length_close;
  int       n_return;
  size_t    s_load;
  int       flag_keep_IDs      =TRUE;
  int       flag_multifile     =FALSE;
  int       flag_multimass     =FALSE;
  int       flag_filefound     =FALSE;
  int       flag_gas           =FALSE;
  int       flag_initpositions =FALSE;
  int       flag_no_velocities =FALSE;
  int       flag_LONGIDS       =FALSE;
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
  SID_fp    fp;
  SID_fp    fp_ids;
  void     *buffer;
  void     *buffer_ids;
  size_t    id_search;
  size_t    id_value;
  size_t    n_buffer;
  int       n_non_zero;
  int       seed=1073743;
  RNG_info *RNG;
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
 
  // First assume that this is not a multifile snapshot ...
  i_file =0;
  n_files=1;
  strcpy(filename,filename_in);
  if(SID_fopen(filename,"r",&fp)){
    flag_filefound=TRUE;
    flag_multifile=FALSE;
  }
  // ... if that doesn't work, check for multi-file
  else{
    sprintf(filename,"%s.%d",filename_in,i_file);
    if(SID_fopen(filename,"r",&fp)){
      flag_filefound=TRUE;
      flag_multifile=TRUE;
    }
    else
      fprintf(stderr,"ERROR: Can not open {%s}!\n",filename);
  }

  // A file was found ... 
  if(flag_filefound){
    SID_log("Reading GADGET binary file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in);
    // This flag is set when we are only reading in
    //   initial particle positions (eg. to set a glass)
    // If TRUE then only positions will be stored and they
    //   will be stored with '_init' at the end of the
    //   variable names
    if(ADaPS_exist(plist->data,"flag_initpositions")){
      flag_initpositions=((int  *)ADaPS_fetch(plist->data,"flag_initpositions"))[0];
      name_initpositions= (char *)ADaPS_fetch(plist->data,"name_initpositions");
    }
    else
      flag_initpositions=FALSE;

    // This flag is set when we are only reading in
    //   initial particle positions (eg. to set a glass)
    if(ADaPS_exist(plist->data,"flag_no_velocities")){
      flag_no_velocities=TRUE;
    }
    else
      flag_no_velocities=FALSE;

    // This flag is set if we want all cores to read all particles
    if(ADaPS_exist(plist->data,"flag_all_read_all"))
      flag_all_read_all=((int  *)ADaPS_fetch(plist->data,"flag_all_read_all"))[0];

    // Set this flag if you want to save RAM and don't need the particle IDs
    if(ADaPS_exist(plist->data,"flag_keep_IDs"))
      flag_keep_IDs=((int  *)ADaPS_fetch(plist->data,"flag_keep_IDs"))[0];
    else
      flag_keep_IDs=TRUE;

    // These values are set when you want to scale the positions by some factor
    if(ADaPS_exist(plist->data,"x_scale")){
      x_scale    =((double  *)ADaPS_fetch(plist->data,"x_scale"))[0];
      flag_xscale=TRUE;
    }
    else
      flag_xscale=FALSE;
    if(ADaPS_exist(plist->data,"y_scale")){
      y_scale    =((double  *)ADaPS_fetch(plist->data,"y_scale"))[0];
      flag_yscale=TRUE;
    }
    else
      flag_yscale=FALSE;
    if(ADaPS_exist(plist->data,"z_scale")){
      z_scale    =((double  *)ADaPS_fetch(plist->data,"z_scale"))[0];
      flag_zscale=TRUE;
    }
    else
      flag_zscale=FALSE;

    // These values are set when you want to apply offsets to the positions
    if(ADaPS_exist(plist->data,"x_offset")){
      x_offset    =((double  *)ADaPS_fetch(plist->data,"x_offset"))[0];
      flag_xoffset=TRUE;
    }
    else
      flag_xoffset=FALSE;
    if(ADaPS_exist(plist->data,"y_offset")){
      y_offset    =((double  *)ADaPS_fetch(plist->data,"y_offset"))[0];
      flag_yoffset=TRUE;
    }
    else
      flag_yoffset=FALSE;
    if(ADaPS_exist(plist->data,"z_offset")){
      z_offset    =((double  *)ADaPS_fetch(plist->data,"z_offset"))[0];
      flag_zoffset=TRUE;
    }
    else
      flag_zoffset=FALSE;

    // These values are set when you want to (subsequently) enforce periodic BCs
    if(ADaPS_exist(plist->data,"x_period")){
      x_period    =((double  *)ADaPS_fetch(plist->data,"x_period"))[0];
      flag_xperiod=TRUE;
    }
    else
      flag_xperiod=FALSE;
    if(ADaPS_exist(plist->data,"y_period")){
      y_period    =((double  *)ADaPS_fetch(plist->data,"y_period"))[0];
      flag_yperiod=TRUE;
    }
    else
      flag_yperiod=FALSE;
    if(ADaPS_exist(plist->data,"z_period")){
      z_period    =((double  *)ADaPS_fetch(plist->data,"z_period"))[0];
      flag_zperiod=TRUE;
    }
    else
      flag_zperiod=FALSE;

    pname=plist->species;

    // Header record length
    SID_log("Reading header...",SID_LOG_OPEN);
    SID_fread_all(&record_length_open,4,1,&fp);
    if(record_length_open!=GADGET_HEADER_SIZE)
      SID_log_warning("Problem with GADGET record size (opening size of header is wrong)",ERROR_LOGIC);
    s_load=0;

    // Number of particles for each species in this file
    n_return=SID_fread_all(n_of_type_tmp,sizeof(unsigned int),N_GADGET_TYPE,&fp); 
    s_load  +=n_return*sizeof(unsigned int);
    for(i=0;i<N_GADGET_TYPE;i++)
      n_of_type[i]=(size_t)n_of_type_tmp[i];

    // Particle mass array 
    n_return =SID_fread_all(mass_array,sizeof(double),N_GADGET_TYPE,&fp); 
    s_load  +=n_return*sizeof(double);

    // Expansion factor (or time) 
    n_return=SID_fread_all(&expansion_factor,sizeof(double),1,&fp);
    s_load+=n_return*sizeof(double);
    ADaPS_store(&(plist->data),(void *)(&expansion_factor),"expansion_factor",ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&(plist->data),(void *)(&expansion_factor),"time",            ADaPS_SCALAR_DOUBLE);

    // Redshift
    n_return=SID_fread_all(&d_value,sizeof(double),1,&fp);
    s_load+=n_return*sizeof(double);
    ADaPS_store(&(plist->data),(void *)(&d_value),"redshift",ADaPS_SCALAR_DOUBLE);

    // Some flags 
    n_return=SID_fread_all(&i_value,sizeof(int),1,&fp);
    s_load+=n_return*sizeof(int);
    ADaPS_store(&(plist->data),(void *)(&i_value),"flag_Sfr",ADaPS_SCALAR_INT);

    n_return=SID_fread_all(&i_value,sizeof(int),1,&fp);
    s_load+=n_return*sizeof(int);
    ADaPS_store(&(plist->data),(void *)(&i_value),"flag_feedback",ADaPS_SCALAR_INT);

    // Number of particles for each species in all files
    n_return=SID_fread_all(n_all_tmp,sizeof(int),N_GADGET_TYPE,&fp);
    for(i=0;i<N_GADGET_TYPE;i++)
      n_all[i]=(size_t)n_all_tmp[i];
    s_load+=n_return*sizeof(unsigned int);

    // Another flag
    n_return=SID_fread_all(&i_value,sizeof(int),1,&fp);
    s_load +=n_return*sizeof(int);
    ADaPS_store(&(plist->data),(void *)(&i_value),"flag_cooling",ADaPS_SCALAR_INT);

    // Number of files in this snapshot 
    n_return=SID_fread_all(&n_files_in,sizeof(int),1,&fp);
    s_load+=n_return*sizeof(int);
    ADaPS_store(&(plist->data),(void *)(&n_files),"n_files",ADaPS_SCALAR_INT);
    if(flag_multifile)
      n_files=n_files_in;

    // Box size -- n.b.: we need h_Hubble before we can store box_size ... do so later 
    n_return=SID_fread_all(&box_size,sizeof(double),1,&fp);
    s_load +=n_return*sizeof(double);

    // Cosmology 

    // Omega_o 
    n_return =SID_fread_all(&d_value,sizeof(double),1,&fp);
    s_load  +=n_return*sizeof(double);
    ADaPS_store(&(plist->data),(void *)(&d_value),"Omega_M",ADaPS_SCALAR_DOUBLE);

    // Omega_Lambda 
    n_return =SID_fread_all(&d_value, sizeof(double),1,&fp);
    s_load  +=n_return*sizeof(double);
    ADaPS_store(&(plist->data),(void *)(&d_value),"Omega_Lambda",ADaPS_SCALAR_DOUBLE);

    // Hubble parameter 
    n_return =SID_fread_all(&h_Hubble,sizeof(double),1,&fp);
    s_load  +=n_return*sizeof(double);
    if(h_Hubble<1e-10) h_Hubble=1.;
    if(check_mode_for_flag(mode,READ_GADGET_MODE_NO_HUBBLE))
      h_Hubble=1.;
    box_size*=plist->length_unit/h_Hubble;
    ADaPS_store(&(plist->data),(void *)(&box_size),"box_size",ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&(plist->data),(void *)(&h_Hubble),"h_Hubble",ADaPS_SCALAR_DOUBLE);

    // Last few things in the header
    n_return =SID_fread_all(&flag_metals,sizeof(int),1,&fp);
    s_load  +=n_return*sizeof(int);
    n_return =SID_fread_all(&flag_ages,sizeof(int),1,&fp);
    s_load  +=n_return*sizeof(int);
    n_return =SID_fread_all(&n_all_tmp,sizeof(int),N_GADGET_TYPE,&fp);
    s_load  +=n_return*sizeof(int);
    n_return =SID_fread_all(&flag_entropyICs,sizeof(int),1,&fp);
    s_load  +=n_return*sizeof(int);
    flag_LONGIDS=FALSE;
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_all_tmp[i]>0){
        n_all[i]+=(((size_t)n_all_tmp[i]) << 32);
        flag_LONGIDS=TRUE;
      }
    }

    // Count the total number of particles
    n_particles_all =0;
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

    // Skip unused space
    SID_fread_all(unused,1,GADGET_HEADER_SIZE-s_load,&fp);
    SID_fread_all(&record_length_close,4,1,&fp);
    if(record_length_open!=record_length_close)
      SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);

    // Close file
    SID_fclose(&fp);
    SID_log("Done.",SID_LOG_CLOSE);

    // Store mass array
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_all[i]>0){
        mass_array[i]*=plist->mass_unit/h_Hubble;
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

    // Define the read mode
    read_mode=READ_GADGET_DEFAULT;
    if(ADaPS_exist(plist->data,"flag_read_scatter"))
      read_mode=READ_GADGET_RANDOM;
    if(ADaPS_exist(plist->data,"rank_x_min"))
      read_mode=read_mode|READ_GADGET_X_MIN;
    if(ADaPS_exist(plist->data,"rank_x_max"))
      read_mode=read_mode|READ_GADGET_X_MAX;
    if(ADaPS_exist(plist->data,"rank_y_min"))
      read_mode=read_mode|READ_GADGET_Y_MIN;
    if(ADaPS_exist(plist->data,"rank_y_max"))
      read_mode=read_mode|READ_GADGET_Y_MAX;
    if(ADaPS_exist(plist->data,"rank_z_min"))
      read_mode=read_mode|READ_GADGET_Z_MIN;
    if(ADaPS_exist(plist->data,"rank_z_max"))
      read_mode=read_mode|READ_GADGET_Z_MAX;

    // Which particles are to be read
    //  into this rank if other constraints
    //  are not present?
    if(SID.n_proc>1){
      if(flag_all_read_all){
        id_min=0;
        id_max=n_particles_all-1;
      }
      else{
        if(SID.My_rank==0){
          id_min=0;
          id_max=(n_particles_all/(size_t)SID.n_proc)*(size_t)(SID.My_rank+1)-1;
        }
        else if(SID.My_rank==SID.n_proc-1){
          id_min=(n_particles_all/(size_t)SID.n_proc)*(size_t)(SID.My_rank);
          id_max=n_particles_all-1;
        }
        else{
          id_min=(n_particles_all/(size_t)SID.n_proc)*(size_t)(SID.My_rank);
          id_max=(n_particles_all/(size_t)SID.n_proc)*(size_t)(SID.My_rank+1)-1;
        }
      }
    }
    else{
      id_min=0;
      id_max=n_particles_all-1;
      flag_all_read_all=TRUE;
    }
    ADaPS_store(&(plist->data),(void *)(&id_min),"rank_id_min",ADaPS_SCALAR_SIZE_T);
    ADaPS_store(&(plist->data),(void *)(&id_max),"rank_id_max",ADaPS_SCALAR_SIZE_T);

    // Other constraints
    if(ADaPS_exist(plist->data,"rank_x_min"))
      x_min_read=((double *)ADaPS_fetch(plist->data,"rank_x_min"))[0];
    else
      x_min_read=1e8*M_PER_MPC;
    if(ADaPS_exist(plist->data,"rank_x_max"))
      x_max_read=((double *)ADaPS_fetch(plist->data,"rank_x_max"))[0];
    else
      x_max_read=-1e8*M_PER_MPC;
    if(ADaPS_exist(plist->data,"rank_y_min"))
      y_min_read=((double *)ADaPS_fetch(plist->data,"rank_y_min"))[0];
    else
      y_min_read=1e8*M_PER_MPC;
    if(ADaPS_exist(plist->data,"rank_y_max"))
      y_max_read=((double *)ADaPS_fetch(plist->data,"rank_y_max"))[0];
    else
      y_max_read=-1e8*M_PER_MPC;
    if(ADaPS_exist(plist->data,"rank_z_min"))
      z_min_read=((double *)ADaPS_fetch(plist->data,"rank_z_min"))[0];
    else
      z_min_read=1e8*M_PER_MPC;
    if(ADaPS_exist(plist->data,"rank_z_max"))
      z_max_read=((double *)ADaPS_fetch(plist->data,"rank_z_max"))[0];
    else
      z_max_read=-1e8*M_PER_MPC;

/*
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
      x_min_bcast=x_min_read*h_Hubble/M_PER_MPC;MPI_Bcast(&x_min_bcast,1,MPI_DOUBLE,i_rank,MPI_COMM_WORLD);
      x_max_bcast=x_max_read*h_Hubble/M_PER_MPC;MPI_Bcast(&x_max_bcast,1,MPI_DOUBLE,i_rank,MPI_COMM_WORLD);
      y_min_bcast=y_min_read*h_Hubble/M_PER_MPC;MPI_Bcast(&y_min_bcast,1,MPI_DOUBLE,i_rank,MPI_COMM_WORLD);
      y_max_bcast=y_max_read*h_Hubble/M_PER_MPC;MPI_Bcast(&y_max_bcast,1,MPI_DOUBLE,i_rank,MPI_COMM_WORLD);
      z_min_bcast=z_min_read*h_Hubble/M_PER_MPC;MPI_Bcast(&z_min_bcast,1,MPI_DOUBLE,i_rank,MPI_COMM_WORLD);
      z_max_bcast=z_max_read*h_Hubble/M_PER_MPC;MPI_Bcast(&z_max_bcast,1,MPI_DOUBLE,i_rank,MPI_COMM_WORLD);
      if(SID.I_am_Master)
        fprintf(stderr,"Rank %3d: %le-%le,%le-%le,%le-%le\n",i_rank,x_min_bcast,x_max_bcast,y_min_bcast,y_max_bcast,z_min_bcast,z_max_bcast);
    }
*/

    // PERFORM DOMAIN DECOMPOSITION
    if(SID.n_proc>1)
      SID_log("Performing domain decomposition...",SID_LOG_OPEN|SID_LOG_TIMER);
    else
      SID_log("Initializing read...",SID_LOG_OPEN|SID_LOG_TIMER);
    for(i=0;i<N_GADGET_TYPE;i++){
      list_file_rank[i]     =(size_t **)SID_malloc(sizeof(size_t *)*n_files);
      n_of_type_file_rank[i]=(size_t  *)SID_malloc(sizeof(size_t  )*n_files);
      for(i_file=0;i_file<n_files;i_file++)
        n_of_type_file_rank[i][i_file]=0;
      n_of_type_rank[i]=0;
    }
    n_rank=(size_t *)SID_malloc(sizeof(size_t)*n_files);
    n_particles_mass_rank=0;
    
    // Are we reading pre-marked particles?
    flag_read_marked=FALSE;
    for(i=0;i<N_GADGET_TYPE;i++){
      if(ADaPS_exist(plist->data,"mark_%s",plist->species[i]))
        flag_read_marked=TRUE;
    }

    // Are we reading a catalog of particles?
    flag_read_catalog=FALSE;
    if(ADaPS_exist(plist->data,"read_catalog")){
      if(flag_read_marked)
        SID_trap_error("Can't read a catalog and mark list at the same time (yet).",ERROR_LOGIC);
      read_catalog=(char *)ADaPS_fetch(plist->data,"read_catalog");
      flag_read_catalog=TRUE;
    }
    
    if(flag_read_marked || flag_read_catalog)
      buffer_ids=(void *)SID_malloc(sizeof(long long)*GADGET_BUFFER_SIZE);

    buffer=(void *)SID_malloc(MAX(3*sizeof(float),sizeof(long long))*GADGET_BUFFER_SIZE);
    n_particles_rank=0;
    for(i_file=0,i_p_temp=0;i_file<n_files;i_file++){
      // Open file
      if(flag_multifile){
        sprintf(filename,"%s.%d",filename_in,i_file);
        SID_log("Processing file #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file);
      }
      else
        sprintf(filename,"%s",filename_in);
      SID_fopen(filename,"r",&fp);

      // Read n_of_type array for this file
      //  and skip the rest of the header
      SID_fread_all(&record_length_open,4,1,&fp);
      n_return=SID_fread_all(n_of_type_tmp,sizeof(unsigned int),N_GADGET_TYPE,&fp);
      s_load  =n_return*sizeof(unsigned int);
      n_return=SID_fread_all(unused,sizeof(int),(GADGET_HEADER_SIZE-s_load)/(sizeof(int)),&fp);
      SID_fread_all(&record_length_close,4,1,&fp);
      SID_fread_all(&record_length_open,4,1,&fp);
      for(i=0,n_particles_file=0;i<N_GADGET_TYPE;i++){
        n_of_type_file[i] =(size_t)n_of_type_tmp[i];
        n_particles_file +=n_of_type_file[i];
      }

      // If we are reading marked particles only, then we need another fp to 
      //   read ids for domain decomposition and we need to skip to the ids
      if(flag_read_marked || flag_read_catalog){
        SID_fopen(filename,"r",&fp_ids);
        s_load=4+GADGET_HEADER_SIZE+4;                   // skip header
        s_load+=(4+4);                                   // skip position block size
        for(i=0;i<N_GADGET_TYPE;i++){
          if(n_of_type_file[i]>0)
            s_load+=2*3*n_of_type_file[i]*sizeof(float); // skip positions & velocities
        }
        s_load+=(4+4);                                   // skip velocity block size
        s_load+=4;                                       // skip ids block size
        SID_fseek(&fp_ids,1,s_load,SID_SEEK_CUR);
      }

      // Read positions and count the number needed by this rank from this file
      n_rank[i_file]=0;
      if(ADaPS_exist(plist->data,"flag_read_scatter")){
        RNG=(RNG_info *)SID_malloc(sizeof(RNG_info));
        init_RNG(&seed,RNG,RNG_GLOBAL);
      }
      else
        RNG=NULL;
      for(i=0,i_p=i_p_temp;i<N_GADGET_TYPE;i++){
        if(flag_read_marked){
          id_list  =(size_t *)ADaPS_fetch(plist->data,"mark_%s",plist->species[i]);
          if(ADaPS_exist(plist->data,"n_local_mark_%s",plist->species[i]))
            n_id_list=((size_t *)ADaPS_fetch(plist->data,"n_local_mark_%s",plist->species[i]))[0];
          else
            n_id_list=0;
          merge_sort((void *)id_list,(size_t)n_id_list,&id_list_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
          if(n_id_list==0)
            mark_mode=READ_GADGET_NONE; 
          else
            mark_mode=READ_GADGET_MARKED; 
        }
        else if(flag_read_catalog){
          if(ADaPS_exist(plist->data, "particle_ids_%s",read_catalog))
            id_list  =(size_t *)ADaPS_fetch(plist->data, "particle_ids_%s",read_catalog);
          else
            SID_trap_error("variable particle_ids_%s no present in data structure!",ERROR_LOGIC,read_catalog);
          if(ADaPS_exist(plist->data, "n_particles_%s",read_catalog))
            n_id_list=((size_t *)ADaPS_fetch(plist->data,"n_particles_%s",read_catalog))[0];
          else
            SID_trap_error("variable n_particles_%s no present in data structure!",ERROR_LOGIC,read_catalog);
          merge_sort((void *)id_list,(size_t)n_id_list,&id_list_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
          if(n_id_list==0)
            mark_mode=READ_GADGET_NONE; 
          else
            mark_mode=READ_GADGET_MARKED; 
        }
        else{
          mark_mode=READ_GADGET_DEFAULT;
          id_list      =NULL;
          id_list_index=NULL;
          n_id_list    =0;
        }
        for(j=0,k=0;j<n_of_type_file[i];){
          n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file[i]-j);
          if(flag_read_marked || flag_read_catalog){
            if(flag_LONGIDS)
              SID_fread_all(buffer_ids,sizeof(long long),n_buffer,&fp_ids);
            else
              SID_fread_all(buffer_ids,sizeof(int),      n_buffer,&fp_ids);
          }
          else
            SID_fread_all(buffer,sizeof(float),3*n_buffer,&fp);
          for(jj=0;jj<n_buffer;jj++,j++){
            if(flag_read_catalog){
              if(flag_LONGIDS)
                id_value=(size_t)(((long long *)buffer_ids)[jj]);
              else
                id_value=(size_t)(((int       *)buffer_ids)[jj]);              
              d1_value=0.;
              d2_value=0.;
              d3_value=0.;
            }
            else{
              f1_value=((float *)buffer)[3*jj+0];
              f2_value=((float *)buffer)[3*jj+1];
              f3_value=((float *)buffer)[3*jj+2];
              d1_value=((double)(f1_value))*(plist->length_unit/h_Hubble);
              d2_value=((double)(f2_value))*(plist->length_unit/h_Hubble);
              d3_value=((double)(f3_value))*(plist->length_unit/h_Hubble);
              if(flag_xscale)                        d1_value*=x_scale;
              if(flag_yscale)                        d2_value*=y_scale;
              if(flag_zscale)                        d3_value*=z_scale;
              if(flag_xoffset)                       d1_value+=x_offset;
              if(flag_yoffset)                       d2_value+=y_offset;
              if(flag_zoffset)                       d3_value+=z_offset;
              if(flag_xperiod && d1_value<0.)        d1_value+=x_period;
              if(flag_xperiod && d1_value>=x_period) d1_value-=x_period;
              if(flag_yperiod && d2_value<0.)        d2_value+=y_period;
              if(flag_yperiod && d2_value>=y_period) d2_value-=y_period;
              if(flag_zperiod && d3_value<0.)        d3_value+=z_period;
              if(flag_zperiod && d3_value>=z_period) d3_value-=z_period;
              if(flag_read_marked || flag_read_catalog){
                if(flag_LONGIDS)
                  id_value=(size_t)(((long long *)buffer_ids)[jj]);
                else
                  id_value=(size_t)(((int       *)buffer_ids)[jj]);              
              }
              else
                id_value=0;
            }

            // Judge membership here
            if(assign_particle_to_rank(plist,
                                       d1_value,
                                       d2_value,
                                       d3_value,
                                       id_list,
                                       id_list_index,
                                       n_id_list,
                                       id_value,
                                       i_p,
                                       RNG,
                                       read_mode,mark_mode,
                                       x_min_read,x_max_read,
                                       y_min_read,y_max_read,
                                       z_min_read,z_max_read,
                                       id_min,id_max))
              k++;
            i_p++;
          }
        }
        if(flag_read_marked || flag_read_catalog)
          SID_free((void **)&id_list_index);
        n_of_type_rank[i]+=k;
        n_particles_rank +=k;
        if(mass_array[i]==0.)
          n_particles_mass_rank+=k;
        n_of_type_file_rank[i][i_file]=k;
        n_rank[i_file]+=k;
      } // i

      // Rewind and skip the header -- we need the positions again
      SID_frewind(&fp);
      SID_fread_all(&record_length_close,4,1,&fp);
      SID_fread_all(unused,sizeof(int),GADGET_HEADER_SIZE/(sizeof(int)),&fp);
      SID_fread_all(&record_length_open,4,1,&fp);
      SID_fread_all(&record_length_close,4,1,&fp);

      if(flag_multifile)
        SID_log("membership counted...",SID_LOG_CONTINUE);

      // Skip to the ids for the second fp (if needed)
      if(flag_read_marked || flag_read_catalog)
        SID_fseek(&fp_ids,1,s_load,SID_SEEK_CUR);
      
      // Generate array listing particles needed by this rank from this file
      if(ADaPS_exist(plist->data,"flag_read_scatter"))
        init_RNG(&seed,RNG,RNG_GLOBAL);
      else
        RNG=NULL;
      for(i=0,i_p=i_p_temp;i<N_GADGET_TYPE;i++){
        if(flag_read_marked){
          id_list  =(size_t *)ADaPS_fetch(plist->data,"mark_%s",plist->species[i]);
          if(ADaPS_exist(plist->data,"n_local_mark_%s",plist->species[i]))
            n_id_list=((size_t *)ADaPS_fetch(plist->data,"n_local_mark_%s",plist->species[i]))[0];
          else
            n_id_list=0;
          merge_sort((void *)id_list,(size_t)n_id_list,&id_list_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
          if(n_id_list==0)
            mark_mode=READ_GADGET_NONE; 
          else
            mark_mode=READ_GADGET_MARKED; 
        }
        else if(flag_read_catalog){
          id_list  =(size_t *)ADaPS_fetch(plist->data, "particle_ids_%s",read_catalog);
          n_id_list=((size_t *)ADaPS_fetch(plist->data,"n_particles_%s",read_catalog))[0];
          merge_sort((void *)id_list,(size_t)n_id_list,&id_list_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
          if(n_id_list==0)
            mark_mode=READ_GADGET_NONE; 
          else
            mark_mode=READ_GADGET_MARKED; 
        }
        else{
          mark_mode=READ_GADGET_DEFAULT;
          id_list  =NULL;
          n_id_list=0;
        }
        list_file_rank[i][i_file]=
          (size_t *)SID_malloc(sizeof(size_t)*MAX(1,n_of_type_file_rank[i][i_file]));
        for(j=0,k=0;j<n_of_type_file[i];){
          n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file[i]-j);
          if(flag_read_marked || flag_read_catalog){
            if(flag_LONGIDS)
              SID_fread_all(buffer_ids,sizeof(long long),n_buffer,&fp_ids);
            else
              SID_fread_all(buffer_ids,sizeof(int),      n_buffer,&fp_ids);
          }
          else
            SID_fread_all(buffer,sizeof(float),3*n_buffer,&fp);
          for(jj=0;jj<n_buffer;jj++,j++){
            if(flag_read_marked || flag_read_catalog){
              if(flag_LONGIDS)
                id_value=(size_t)(((long long *)buffer_ids)[jj]);
              else
                id_value=(size_t)(((int       *)buffer_ids)[jj]);              
              d1_value=0.;
              d2_value=0.;
              d3_value=0.;
            }
            else{
              f1_value=((float *)buffer)[3*jj+0];
              f2_value=((float *)buffer)[3*jj+1];
              f3_value=((float *)buffer)[3*jj+2];
              d1_value=((double)(f1_value))*(plist->length_unit/h_Hubble);
              d2_value=((double)(f2_value))*(plist->length_unit/h_Hubble);
              d3_value=((double)(f3_value))*(plist->length_unit/h_Hubble);
              if(flag_xscale)                        d1_value*=x_scale;
              if(flag_yscale)                        d2_value*=y_scale;
              if(flag_zscale)                        d3_value*=z_scale;
              if(flag_xoffset)                       d1_value+=x_offset;
              if(flag_yoffset)                       d2_value+=y_offset;
              if(flag_zoffset)                       d3_value+=z_offset;
              if(flag_xperiod && d1_value<0.)        d1_value+=x_period;
              if(flag_xperiod && d1_value>=x_period) d1_value-=x_period;
              if(flag_yperiod && d2_value<0.)        d2_value+=y_period;
              if(flag_yperiod && d2_value>=y_period) d2_value-=y_period;
              if(flag_zperiod && d3_value<0.)        d3_value+=z_period;
              if(flag_zperiod && d3_value>=z_period) d3_value-=z_period;
              if(flag_read_marked || flag_read_catalog){
                if(flag_LONGIDS)
                  id_value=(size_t)(((long long *)buffer_ids)[jj]);
                else
                  id_value=(size_t)(((int       *)buffer_ids)[jj]);              
              }
              else
                id_value=0;
            }
          
            // Judge membership here
            if(assign_particle_to_rank(plist,
                                       d1_value,
                                       d2_value,
                                       d3_value,
                                       id_list,
                                       id_list_index,
                                       n_id_list,
                                       id_value,
                                       i_p,
                                       RNG,
                                       read_mode,mark_mode,
                                       x_min_read,x_max_read,
                                       y_min_read,y_max_read,
                                       z_min_read,z_max_read,
                                       id_min,id_max)){
              list_file_rank[i][i_file][k]=j;
              k++;
            }
            i_p++;
          }
        }
        if(flag_read_marked || flag_read_catalog)
          SID_free((void **)&id_list_index);
      }
      i_p_temp=i_p;
      // Close file
      SID_fclose(&fp);
      if(flag_read_marked || flag_read_catalog)
        SID_fclose(&fp_ids);
      if(ADaPS_exist(plist->data,"flag_read_scatter"))
        SID_free((void **)&RNG);
      if(flag_multifile)
        SID_log("Done.",SID_LOG_CLOSE);
    }
    if(ADaPS_exist(plist->data,"flag_read_scatter"))
      SID_free((void **)&RNG);
    if(n_of_type_rank[GADGET_TYPE_GAS]==0)
      flag_gas=FALSE;
    n_particles_all=0;
#ifdef USE_MPI
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
      n_particles_bcast=n_particles_rank;
      MPI_Bcast(&n_particles_bcast,1,MPI_SIZE_T,i_rank,MPI_COMM_WORLD);
      SID_log("Number of particles on rank #%d:\t%lld",SID_LOG_COMMENT,i_rank,n_particles_bcast);
      n_particles_all+=n_particles_bcast;
    }
    SID_log("Done.",SID_LOG_CLOSE);
#else
    n_particles_all=n_particles_rank;
    SID_log("%lld particles will be read...Done.",SID_LOG_CLOSE,n_particles_rank);
#endif

    // Allocate data arrays
    x_array=(REAL *)SID_malloc(sizeof(REAL)*(size_t)(MAX(1,n_particles_rank)));
    y_array=(REAL *)SID_malloc(sizeof(REAL)*(size_t)(MAX(1,n_particles_rank)));
    z_array=(REAL *)SID_malloc(sizeof(REAL)*(size_t)(MAX(1,n_particles_rank)));
    if(!flag_initpositions){
      if(!flag_no_velocities){
        vx_array=(REAL   *)SID_malloc(sizeof(REAL  )*(size_t)(MAX(1,n_particles_rank)));
        vy_array=(REAL   *)SID_malloc(sizeof(REAL  )*(size_t)(MAX(1,n_particles_rank)));
        vz_array=(REAL   *)SID_malloc(sizeof(REAL  )*(size_t)(MAX(1,n_particles_rank)));
      }
      if(flag_keep_IDs)
        id_array=(size_t *)SID_malloc(sizeof(size_t)*(size_t)(MAX(1,n_particles_rank)));
      if(flag_multimass)
        M_array  =(double *)SID_malloc(sizeof(double)*(size_t)(MAX(1,n_particles_mass_rank)));
      if(flag_gas){
        u_array     =(REAL *)SID_malloc(sizeof(REAL)*(size_t)(MAX(1,n_of_type_rank[GADGET_TYPE_GAS])));
        rho_array   =(REAL *)SID_malloc(sizeof(REAL)*(size_t)(MAX(1,n_of_type_rank[GADGET_TYPE_GAS])));
        T_array     =(REAL *)SID_malloc(sizeof(REAL)*(size_t)(MAX(1,n_of_type_rank[GADGET_TYPE_GAS])));
        smooth_array=(REAL *)SID_malloc(sizeof(REAL)*(size_t)(MAX(1,n_of_type_rank[GADGET_TYPE_GAS])));
      }
    }

    // Set offsets
    for(i=0,j=0,k=0;i<N_GADGET_TYPE;i++){
      k_offset[i]     =k;
      k_offset_mass[i]=j;
      k+=n_of_type_rank[i];
      if(mass_array[i]==0.)
        j+=n_of_type_rank[i];
    }
    // Read data
    for(i_file=0;i_file<n_files;i_file++){
      if(n_files>1)
        SID_log("Reading file %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file+1,n_files);

      if(flag_multifile)
        sprintf(filename,"%s.%d",filename_in,i_file);
      else
        sprintf(filename,"%s",filename_in);
      SID_fopen(filename,"r",&fp);

      // Read n_of_type array for this file
      //  and skip the rest of the header
      SID_fread_all(&record_length_open,4,1,&fp);
      n_return=SID_fread_all(n_of_type_tmp,sizeof(unsigned int),N_GADGET_TYPE,&fp);
      for(i=0;i<N_GADGET_TYPE;i++)
        n_of_type_file[i]=(size_t)n_of_type_tmp[i];
      s_load  =n_return*sizeof(unsigned int);
      SID_fread_all(unused,sizeof(int),(GADGET_HEADER_SIZE-s_load)/(sizeof(int)),&fp);
      SID_fread_all(&record_length_close,4,1,&fp);
        
      // Read positions
      SID_log("Reading positions...",SID_LOG_OPEN|SID_LOG_TIMER);
      SID_fread_all(&record_length_open,4,1,&fp);
      for(i=0;i<N_GADGET_TYPE;i++){
        for(j=0,k=0;j<n_of_type_file[i];){
          n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file[i]-j);
          SID_fread_all(buffer,sizeof(float),3*n_buffer,&fp);
          for(jj=0;jj<n_buffer;jj++,j++){
            f1_value=((float *)buffer)[3*jj+0];
            f2_value=((float *)buffer)[3*jj+1];
            f3_value=((float *)buffer)[3*jj+2];
            if(n_of_type_file_rank[i][i_file]>0 && k<n_of_type_file_rank[i][i_file]){
              if(j==list_file_rank[i][i_file][k]){
                x_array[k+k_offset[i]]=((REAL)(f1_value))*(REAL)(plist->length_unit/h_Hubble);
                y_array[k+k_offset[i]]=((REAL)(f2_value))*(REAL)(plist->length_unit/h_Hubble);
                z_array[k+k_offset[i]]=((REAL)(f3_value))*(REAL)(plist->length_unit/h_Hubble);
                if(flag_xscale)  x_array[k+k_offset[i]]*=(REAL)x_scale;
                if(flag_yscale)  y_array[k+k_offset[i]]*=(REAL)y_scale;
                if(flag_zscale)  z_array[k+k_offset[i]]*=(REAL)z_scale;
                if(flag_xoffset) x_array[k+k_offset[i]]+=(REAL)x_offset;
                if(flag_yoffset) y_array[k+k_offset[i]]+=(REAL)y_offset;
                if(flag_zoffset) z_array[k+k_offset[i]]+=(REAL)z_offset;
                if(flag_xperiod && x_array[k+k_offset[i]]<0.)
                  x_array[k+k_offset[i]]+=(REAL)x_period;
                if(flag_xperiod && x_array[k+k_offset[i]]>=(REAL)x_period)
                  x_array[k+k_offset[i]]-=(REAL)x_period;
                if(flag_yperiod && y_array[k+k_offset[i]]<0.)
                  y_array[k+k_offset[i]]+=(REAL)y_period;
                if(flag_yperiod && y_array[k+k_offset[i]]>=(REAL)y_period) 
                  y_array[k+k_offset[i]]-=(REAL)y_period;
                if(flag_zperiod && z_array[k+k_offset[i]]<0.)
                  z_array[k+k_offset[i]]+=(REAL)z_period;
                if(flag_zperiod && z_array[k+k_offset[i]]>=(REAL)z_period)
                  z_array[k+k_offset[i]]-=(REAL)z_period;
                k++;
              }
            }
          }
        }
      }
      SID_fread_all(&record_length_close,4,1,&fp);
      if(record_length_open!=record_length_close)
        SID_log_warning("Problem with GADGET record size (close of positions)",ERROR_LOGIC);
      SID_log("Done.",SID_LOG_CLOSE);

      // Read the rest only if flag_initpositions!=TRUE
      if(!flag_initpositions){
        if(!flag_no_velocities){
          // Read velocities
          SID_log("Reading velocities...",SID_LOG_OPEN|SID_LOG_TIMER);
          SID_fread_all(&record_length_open,4,1,&fp);
          for(i=0;i<N_GADGET_TYPE;i++){
            for(j=0,k=0;j<n_of_type_file[i];){
              n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file[i]-j);
              SID_fread_all(buffer,sizeof(float),3*n_buffer,&fp);
              for(jj=0;jj<n_buffer;jj++,j++){
                f1_value=((float *)buffer)[3*jj+0];
                f2_value=((float *)buffer)[3*jj+1];
                f3_value=((float *)buffer)[3*jj+2];
                if(n_of_type_file_rank[i][i_file]>0 && k<n_of_type_file_rank[i][i_file]){
                  if(j==list_file_rank[i][i_file][k]){
                    vx_array[k+k_offset[i]]=((REAL)(f1_value))*(REAL)(plist->velocity_unit*sqrt(expansion_factor));
                    vy_array[k+k_offset[i]]=((REAL)(f2_value))*(REAL)(plist->velocity_unit*sqrt(expansion_factor));
                    vz_array[k+k_offset[i]]=((REAL)(f3_value))*(REAL)(plist->velocity_unit*sqrt(expansion_factor));
                    k++;
                  }
                }
              }
            }
          }
          SID_fread_all(&record_length_close,4,1,&fp);
          if(record_length_open!=record_length_close)
            SID_log_warning("Problem with GADGET record size (close of velocities)",ERROR_LOGIC);
          SID_log("Done.",SID_LOG_CLOSE);
        }
        else{
          SID_fread_all(&record_length_open, 4,1,&fp);
          SID_fskip(record_length_open,1,&fp); 
          SID_fread_all(&record_length_close,4,1,&fp);
        }
        
        // Read ids
        if(flag_LONGIDS)        
          SID_log("Reading IDs (long long)...",SID_LOG_OPEN|SID_LOG_TIMER);
        else
          SID_log("Reading IDs (int)...",SID_LOG_OPEN|SID_LOG_TIMER);
        if(!flag_keep_IDs)
          SID_log("(will be discarded)...",SID_LOG_CONTINUE);
        SID_fread_all(&record_length_open,4,1,&fp);
        for(i=0;i<N_GADGET_TYPE;i++){
          for(j=0,k=0;j<n_of_type_file[i];){
            n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file[i]-j);
            if(flag_LONGIDS)
              SID_fread_all(buffer,sizeof(long long),n_buffer,&fp);
            else
              SID_fread_all(buffer,sizeof(int),n_buffer,&fp);
            for(jj=0;jj<n_buffer;jj++,j++){
              if(flag_LONGIDS)
                id_value_LL=((long long *)buffer)[jj];
              else
                id_value_I =((int *)buffer)[jj];
              if(n_of_type_file_rank[i][i_file]>0 && k<n_of_type_file_rank[i][i_file]){
                if(j==list_file_rank[i][i_file][k]){
                  if(flag_keep_IDs){
                    if(flag_LONGIDS)
                      id_array[k+k_offset[i]]=(size_t)id_value_LL;
                    else
                      id_array[k+k_offset[i]]=(size_t)id_value_I;
                  }
                  k++;
                }
              }
            }
          }
        }
        SID_fread_all(&record_length_close,4,1,&fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of ids)",ERROR_LOGIC);
        SID_log("Done.",SID_LOG_CLOSE);

        // Read masses (if needed)
        if(flag_multimass){
#ifdef USE_MPI
          MPI_Allreduce(&n_particles_mass_rank,&n_particles_mass,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
#else
          n_particles_mass=n_particles_mass_rank;
#endif
          SID_log("Reading mass array (%d multimass particles)...",SID_LOG_OPEN|SID_LOG_TIMER,n_particles_mass);
          SID_fread_all(&record_length_open,4,1,&fp);
          for(i=0;i<N_GADGET_TYPE;i++){
            if(mass_array[i]==0.){
              for(j=0,k=0;j<n_of_type_file[i];){
                n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file[i]-j);
                SID_fread_all(buffer,sizeof(float),n_buffer,&fp);
                for(jj=0;jj<n_buffer;jj++,j++){
                  f1_value=((float *)buffer)[jj];
                  if(n_of_type_file_rank[i][i_file]>0 && k<n_of_type_file_rank[i][i_file]){
                    if(j==list_file_rank[i][i_file][k]){
                      M_array[k+k_offset_mass[i]]=(double)(f1_value)*(double)(plist->mass_unit/h_Hubble);
                      k++;
                    }
                  }
                }
              }
            }
          }
          SID_fread_all(&record_length_close,4,1,&fp);
          if(record_length_open!=record_length_close)
            SID_log_warning("Problem with GADGET record size (close of masses)",ERROR_LOGIC);
          SID_log("Done.",SID_LOG_CLOSE);
        }

        // Read internal energies/densities and compute temperatures for gas
        if(flag_gas && n_of_type_file_rank[GADGET_TYPE_GAS][i_file]>0){
          i=GADGET_TYPE_GAS;
          // Read energies
          SID_log("Reading gas properties...",SID_LOG_OPEN|SID_LOG_TIMER);
          SID_fread_all(&record_length_open,4,1,&fp);
          for(j=0,k=0;j<n_of_type_file[i];){
            n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file[i]-j);
            SID_fread_all(buffer,sizeof(float),n_buffer,&fp);
            for(jj=0;jj<n_buffer;jj++,j++){
              f1_value=((float *)buffer)[jj];
              if(k<n_of_type_file_rank[i][i_file]){
                if(j==list_file_rank[i][i_file][k]){
                  u_array[k+k_offset[i]]=(REAL)(f1_value)*(REAL)(plist->velocity_unit*plist->velocity_unit);
                  k++;
                }
              }
            }
          }
          SID_fread_all(&record_length_close,4,1,&fp);
          if(record_length_open!=record_length_close)
            SID_log_warning("Problem with GADGET record size (close of internal energies)",ERROR_LOGIC);

          // Read densities
          SID_fread_all(&record_length_open,4,1,&fp);
          for(j=0,k=0;j<n_of_type_file[i];){
            n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file[i]-j);
            SID_fread_all(buffer,sizeof(float),n_buffer,&fp);
            for(jj=0;jj<n_buffer;jj++,j++){
              f1_value=((float *)buffer)[jj];
              if(k<n_of_type_file_rank[i][i_file]){
                if(j==list_file_rank[i][i_file][k]){
                  rho_array[k+k_offset[i]]=
                    (REAL)(f1_value)*(REAL)(h_Hubble*h_Hubble*plist->mass_unit/pow(plist->length_unit,3.));
                  k++;
                }
              }
            }
          }
          SID_fread_all(&record_length_close,4,1,&fp);
          if(record_length_open!=record_length_close)
            SID_log_warning("Problem with GADGET record size (close of densities)",ERROR_LOGIC);

          // Read smoothing lengths
          SID_fread_all(&record_length_open,4,1,&fp);
          for(j=0,k=0;j<n_of_type_file[i];){
            n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file[i]-j);
            SID_fread_all(buffer,sizeof(float),n_buffer,&fp);
            for(jj=0;jj<n_buffer;jj++,j++){
              f1_value=((float *)buffer)[jj];
              if(k<n_of_type_file_rank[i][i_file]){
                if(j==list_file_rank[i][i_file][k]){
                  smooth_array[k+k_offset[i]]=
                    ((REAL)(f1_value))*(REAL)(plist->length_unit/h_Hubble);
                  k++;
                }
              }
            }
          }
          SID_fread_all(&record_length_close,4,1,&fp);
          if(record_length_open!=record_length_close)
            SID_log_warning("Problem with GADGET record size (close of smoothing lengths)",ERROR_LOGIC);

          // Compute temperatures for gas
          for(j=0,k=0;j<n_of_type_file[GADGET_TYPE_GAS];j++){
            if(rho_array[k]>0.)
              T_array[k+k_offset[i]]=u_array[k]*MU_MMW*M_PROTON/(1.5*K_BOLTZMANN);    
            else
              T_array[k+k_offset[i]]=0.;
            k++;
          }
          SID_log("Done.",SID_LOG_CLOSE);
        }
      }
      for(i=0,j=0,k=0;i<N_GADGET_TYPE;i++){
        k_offset[i]+=n_of_type_file_rank[i][i_file];
        if(mass_array[i]==0.)
          k_offset_mass[i]+=n_of_type_file_rank[i][i_file];
      }
      // Close file
      SID_fclose(&fp);
      if(n_files>1)
        SID_log("Done.",SID_LOG_CLOSE);
    }
    SID_free((void **)&buffer);
    if(flag_read_marked || flag_read_catalog)
      SID_free((void **)&buffer_ids);
    SID_Barrier(SID.COMM_WORLD);

    // Store everything in the data structure...

    //   ... particle counts ...
    for(i=0,k=0,j=0;i<N_GADGET_TYPE;i++){
      k_offset[i]     =k;
      k_offset_mass[i]=j;
      k+=n_of_type_rank[i];
      if(mass_array[i]==0.)
        j+=n_of_type_rank[i];
#ifdef USE_MPI
      MPI_Allreduce(&(n_of_type_rank[i]),&(n_of_type[i]),1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
#else
      n_of_type[i]=n_of_type_rank[i];
#endif
    }
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type[i]>0){
        if(flag_initpositions){
          ADaPS_store(&(plist->data),(void *)(&(n_of_type_rank[i])),"n_%s",    ADaPS_SCALAR_SIZE_T,name_initpositions);
          ADaPS_store(&(plist->data),(void *)(&(n_of_type[i])),     "n_all_%s",ADaPS_SCALAR_SIZE_T,name_initpositions);
        }
        else{
          ADaPS_store(&(plist->data),(void *)(&(n_of_type_rank[i])),"n_%s",    ADaPS_SCALAR_SIZE_T,pname[i]);
          ADaPS_store(&(plist->data),(void *)(&(n_of_type[i])),     "n_all_%s",ADaPS_SCALAR_SIZE_T,pname[i]);
        }
      }
    }
    ADaPS_store(&(plist->data),(void *)(&n_particles_all),"n_particles_all",ADaPS_SCALAR_SIZE_T);

    //   ... positions ...
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type_rank[i]>0){
        if(flag_initpositions){
          ADaPS_store(&(plist->data),(void *)(&(x_array[k_offset[i]])),"x_%s_init",ADaPS_DEFAULT,name_initpositions);
          ADaPS_store(&(plist->data),(void *)(&(y_array[k_offset[i]])),"y_%s_init",ADaPS_DEFAULT,name_initpositions);
          ADaPS_store(&(plist->data),(void *)(&(z_array[k_offset[i]])),"z_%s_init",ADaPS_DEFAULT,name_initpositions);
        }
        else{
          ADaPS_store(&(plist->data),(void *)(&(x_array[k_offset[i]])),"x_%s",     ADaPS_DEFAULT,pname[i]);
          ADaPS_store(&(plist->data),(void *)(&(y_array[k_offset[i]])),"y_%s",     ADaPS_DEFAULT,pname[i]);
          ADaPS_store(&(plist->data),(void *)(&(z_array[k_offset[i]])),"z_%s",     ADaPS_DEFAULT,pname[i]);
        }
      }
    }

    // Don't store the rest of these quantities if
    //  we are just reading positions
    if(!flag_initpositions){

      //  ... velocities ...
      if(!flag_no_velocities){
        for(i=0;i<N_GADGET_TYPE;i++){
          if(n_of_type_rank[i]>0){
            ADaPS_store(&(plist->data),(void *)(&(vx_array[k_offset[i]])),"vx_%s",ADaPS_DEFAULT,pname[i]);
            ADaPS_store(&(plist->data),(void *)(&(vy_array[k_offset[i]])),"vy_%s",ADaPS_DEFAULT,pname[i]);
            ADaPS_store(&(plist->data),(void *)(&(vz_array[k_offset[i]])),"vz_%s",ADaPS_DEFAULT,pname[i]);
          }
        }
      }

      //  ... masses ...
      for(i=0;i<N_GADGET_TYPE;i++){
        if(n_all[i]>0)
          ADaPS_store(&(plist->data),(void *)(&(mass_array[i])),"mass_array_%s",ADaPS_SCALAR_DOUBLE,pname[i]);
        if(flag_multimass && mass_array[i]==0.){
          ADaPS_store(&(plist->data),(void *)(&(M_array[k_offset_mass[i]])),"M_%s",ADaPS_DEFAULT,pname[i]);
        }
      }

      //  ... ids ...
      if(flag_keep_IDs){
        for(i=0;i<N_GADGET_TYPE;i++){
          if(n_of_type_rank[i]>0)
            ADaPS_store(&(plist->data),(void *)(&(id_array[k_offset[i]])),"id_%s",ADaPS_DEFAULT,pname[i]);
        }
      }

      //  ... gas properties ...
      if(flag_gas && n_of_type[GADGET_TYPE_GAS]>0){
        ADaPS_store(&(plist->data),u_array,  "u_%s",       ADaPS_DEFAULT,pname[GADGET_TYPE_GAS]);
        ADaPS_store(&(plist->data),rho_array,"rho_%s",     ADaPS_DEFAULT,pname[GADGET_TYPE_GAS]);
        ADaPS_store(&(plist->data),T_array,  "T_%s",       ADaPS_DEFAULT,pname[GADGET_TYPE_GAS]);
        ADaPS_store(&(plist->data),u_array,  "r_smooth_%s",ADaPS_DEFAULT,pname[GADGET_TYPE_GAS]);
      }
    }

    // Clean-up
    SID_free((void **)&n_rank);
    for(i=0;i<N_GADGET_TYPE;i++){
      SID_free((void **)&n_of_type_file_rank[i]);
      for(i_file=0;i_file<n_files;i_file++)
        SID_free((void **)&list_file_rank[i][i_file]);
      SID_free((void **)&list_file_rank[i]);
    }
    SID_Barrier(SID.COMM_WORLD);
    SID_log("Done.",SID_LOG_CLOSE);
  }
}

