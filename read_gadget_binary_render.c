#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void read_gadget_binary_render(char       *filename_root_in,
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
  size_t    n_of_type[N_GADGET_TYPE];
  size_t    i_particle;
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
  GBPREAL     *x_array[N_GADGET_TYPE];
  GBPREAL     *y_array[N_GADGET_TYPE];
  GBPREAL     *z_array[N_GADGET_TYPE];
  GBPREAL     *vx_array[N_GADGET_TYPE];
  GBPREAL     *vy_array[N_GADGET_TYPE];
  GBPREAL     *vz_array[N_GADGET_TYPE];
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
  RNG_info *RNG;
 
  // Determine file format
  flag_filefound=init_gadget_read(filename_root_in,snapshot_number,&flag_multifile,&flag_file_type,&header);

  // A file was found ... 
  if(flag_filefound){
    SID_log("Reading GADGET binary file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_root_in);

    pname=plist->species;

    // Number of particles for each species in this file
    for(i=0;i<N_GADGET_TYPE;i++)
      n_of_type[i]=(size_t)header.n_file[i];

    // Expansion factor (or time) 
    ADaPS_store(&(plist->data),(void *)(&(header.time)),"expansion_factor",ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&(plist->data),(void *)(&(header.time)),"time",            ADaPS_SCALAR_DOUBLE);

    // Redshift
    d_value=(double)header.redshift;
    ADaPS_store(&(plist->data),(void *)(&d_value),"redshift",ADaPS_SCALAR_DOUBLE);

    // Number of particles for each species in all files
    for(i=0;i<N_GADGET_TYPE;i++){
      n_all[i]     =(size_t)header.n_all[i];
      mass_array[i]=header.mass_array[i];
    }

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
    h_Hubble=(double)header.h_Hubble; 
    if(h_Hubble<1e-10) h_Hubble=1.;
    box_size=header.box_size*plist->length_unit/h_Hubble;
    ADaPS_store(&(plist->data),(void *)(&box_size),"box_size",ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&(plist->data),(void *)(&h_Hubble),"h_Hubble",ADaPS_SCALAR_DOUBLE);

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
    SID_log("%zd",SID_LOG_CONTINUE,n_particles_all);
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

    for(i=0;i<N_GADGET_TYPE;i++)
      n_of_type_rank[i]=0;

    // Count the number of particles that will be scattered to each rank
    pcounter_info pcounter;
    size_t        k_particle;
    SID_init_pcounter(&pcounter,n_particles_all,10);
    SID_log("Counting the number of particles that will be scattered to each rank...",SID_LOG_OPEN|SID_LOG_TIMER);
    RNG=(RNG_info *)SID_malloc(sizeof(RNG_info));
    init_RNG(&seed,RNG,RNG_GLOBAL);
    for(i_file=0,k_particle=0;i_file<n_files;i_file++){
      read_rank=i_file%SID.n_proc;
      read_rank=0;

      set_gadget_filename(filename_root_in,snapshot_number,i_file,flag_multifile,flag_file_type,filename);

      // Open file and read header
      if(SID.My_rank==read_rank){
        fp=fopen(filename,"r");
        fread(&record_length_open,4,1,fp);
        fread(&header,sizeof(gadget_header_info),1,fp);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
        fclose(fp);
      }
      SID_Bcast(&header,(int)sizeof(gadget_header_info),read_rank,SID.COMM_WORLD);
      for(i=0;i<N_GADGET_TYPE;i++){
        for(i_particle=0;i_particle<header.n_file[i];i_particle++,k_particle++){
          scatter_rank=(int)(random_number(RNG)*(GBPREAL)SID.n_proc);
          if(scatter_rank<0)
            scatter_rank=0;
          else if(scatter_rank>=SID.n_proc)
            scatter_rank=SID.n_proc-1;
          if(scatter_rank==SID.My_rank)
            n_of_type_rank[i]++;
          SID_check_pcounter(&pcounter,k_particle);
        }
      }
    }
    for(i=0,n_particles_rank=0;i<N_GADGET_TYPE;i++)
      n_particles_rank+=n_of_type_rank[i];
    init_RNG(&seed2,RNG,RNG_GLOBAL);
    SID_log("Done.",SID_LOG_CLOSE);

    // Allocate data arrays
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type_rank[i]>0){
        x_array[i] =(GBPREAL   *)SID_malloc(sizeof(GBPREAL)  *(size_t)n_of_type_rank[i]);
        y_array[i] =(GBPREAL   *)SID_malloc(sizeof(GBPREAL)  *(size_t)n_of_type_rank[i]);
        z_array[i] =(GBPREAL   *)SID_malloc(sizeof(GBPREAL)  *(size_t)n_of_type_rank[i]);
/*
        vx_array[i]=(GBPREAL   *)SID_malloc(sizeof(GBPREAL)  *(size_t)n_of_type_rank[i]);
        vy_array[i]=(GBPREAL   *)SID_malloc(sizeof(GBPREAL)  *(size_t)n_of_type_rank[i]);
        vz_array[i]=(GBPREAL   *)SID_malloc(sizeof(GBPREAL)  *(size_t)n_of_type_rank[i]);
*/
        id_array[i]=(size_t *)SID_malloc(sizeof(size_t)*(size_t)n_of_type_rank[i]);
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

      set_gadget_filename(filename_root_in,snapshot_number,i_file,flag_multifile,flag_file_type,filename);
        
      // Open file and read header
      if(SID.My_rank==read_rank){
        fp=fopen(filename,"r");
        fread(&record_length_open,4,1,fp);
        fread(&header,sizeof(gadget_header_info),1,fp);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
      }
      SID_Bcast(&header,(int)sizeof(gadget_header_info),read_rank,SID.COMM_WORLD);
      for(i=0,n_particles_file=0;i<N_GADGET_TYPE;i++)
        n_particles_file+=(size_t)header.n_file[i];

      // Initialize buffer
      SID_log("Initializing buffer...",SID_LOG_OPEN|SID_LOG_TIMER);
      if(SID.My_rank==read_rank)
        fread(&record_length_open,4,1,fp);
      SID_Bcast(&record_length_open,(int)sizeof(int),read_rank,SID.COMM_WORLD);
      buffer=SID_malloc((size_t)record_length_open);
      keep  =(char *)SID_malloc(sizeof(char)*n_particles_file);
      for(i=0,jj=0;i<N_GADGET_TYPE;i++){
        n_keep[i]=0;
        for(j=0,k=0;j<header.n_file[i];j++,jj++){
          scatter_rank=(int)(random_number(RNG)*(GBPREAL)SID.n_proc);
          if(scatter_rank<0)
            scatter_rank=0;
          else if(scatter_rank>=SID.n_proc)
            scatter_rank=SID.n_proc-1;
          if(scatter_rank==SID.My_rank){
            n_keep[i]++;
            n_particles_kept++;
            keep[jj]=TRUE;
          }
          else
            keep[jj]=FALSE;
        }
      }
      SID_log("Done.",SID_LOG_CLOSE);

      // Read positions
      SID_log("Reading positions...",SID_LOG_OPEN|SID_LOG_TIMER);
      if(SID.My_rank==read_rank){
        fread(buffer,(size_t)record_length_open,1,fp);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of positions)",ERROR_LOGIC);
      }
      SID_Barrier(SID.COMM_WORLD);
      SID_Bcast(buffer,(int)record_length_open,read_rank,SID.COMM_WORLD);
      for(i=0,jj=0;i<N_GADGET_TYPE;i++){
        for(j=0,k=0;j<header.n_file[i];j++,jj++){
          if(keep[jj]){
            f1_value=((float *)buffer)[3*jj+0];
            f2_value=((float *)buffer)[3*jj+1];
            f3_value=((float *)buffer)[3*jj+2];
            x_array[i][k+k_offset[i]]=((GBPREAL)(f1_value))*(GBPREAL)(plist->length_unit/h_Hubble);
            y_array[i][k+k_offset[i]]=((GBPREAL)(f2_value))*(GBPREAL)(plist->length_unit/h_Hubble);
            z_array[i][k+k_offset[i]]=((GBPREAL)(f3_value))*(GBPREAL)(plist->length_unit/h_Hubble);
            k++;
          }
        }
        if(k!=n_keep[i])
          SID_trap_error("Particle count mismatch (ie. %d!=%d) during positions read",ERROR_LOGIC,k,n_keep[i]);
      }
      SID_log("Done.",SID_LOG_CLOSE);

      // Read velocities
      /*
        SID_log("Reading velocities...",SID_LOG_OPEN|SID_LOG_TIMER);
        if(SID.My_rank==read_rank){
        fread(&record_length_open,4,1,fp);
        fread(buffer,record_length_open,1,fp);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
        SID_log_warning("Problem with GADGET record size (close of velocities)",ERROR_LOGIC);
        }
        SID_Bcast(&record_length_open,sizeof(int),read_rank,SID.COMM_WORLD);
        SID_Bcast(buffer,record_length_open,read_rank,SID.COMM_WORLD);
        for(i=0,jj=0;i<N_GADGET_TYPE;i++){
        for(j=0,k=0;j<header.n_file[i];j++,jj++){
        if(keep[jj]){
        f1_value=((float *)buffer)[3*jj+0];
        f2_value=((float *)buffer)[3*jj+1];
        f3_value=((float *)buffer)[3*jj+2];
        vx_array[i][k+k_offset[i]]=((GBPREAL)(f1_value))*(GBPREAL)(plist->velocity_unit*sqrt(expansion_factor));
        vy_array[i][k+k_offset[i]]=((GBPREAL)(f2_value))*(GBPREAL)(plist->velocity_unit*sqrt(expansion_factor));
        vz_array[i][k+k_offset[i]]=((GBPREAL)(f3_value))*(GBPREAL)(plist->velocity_unit*sqrt(expansion_factor));
        k++;
        }
        }
        if(k!=n_keep[i])
        SID_trap_error("Particle count mismatch during velocity read",ERROR_LOGIC);
        }
        SID_log("Done.",SID_LOG_CLOSE);
      */
      
      SID_log("Skipping velocities...",SID_LOG_OPEN|SID_LOG_TIMER);
      if(SID.My_rank==read_rank){
        fread(&record_length_open,4,1,fp);
        fseeko(fp,(size_t)record_length_open,SEEK_CUR);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of velocities)",ERROR_LOGIC);
      }
      SID_log("Done.",SID_LOG_CLOSE);
      
      // Read IDs
      SID_log("Reading IDs...",SID_LOG_OPEN|SID_LOG_TIMER);
      if(SID.My_rank==read_rank){
        fread(&record_length_open,4,1,fp);
        fread(buffer,(size_t)record_length_open,1,fp);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of IDs)",ERROR_LOGIC);
      }
      SID_Bcast(&record_length_open,(int)sizeof(int),       read_rank,SID.COMM_WORLD);
      SID_Barrier(SID.COMM_WORLD);
      SID_Bcast(buffer,             (int)record_length_open,read_rank,SID.COMM_WORLD);

      // Decide what kind of IDs we have
      if(record_length_open/(int)n_particles_file==sizeof(long long)){
        SID_log("(long long)...",SID_LOG_CONTINUE);
        flag_LONGIDs=TRUE;
      }
      else{
        SID_log("(int)...",SID_LOG_CONTINUE);
        flag_LONGIDs=FALSE;
      }

      // Process IDs
      for(i=0,jj=0;i<N_GADGET_TYPE;i++){
        for(j=0,k=0;j<header.n_file[i];j++,jj++){
          if(keep[jj]){
            switch(flag_LONGIDs){
            case TRUE:
              id_test=(size_t)((long long *)buffer)[jj];
              break;
            case FALSE:
              id_test=(size_t)((int *)buffer)[jj];
              break;
            }
            id_array[i][k+k_offset[i]]=id_test;
            k++;
          }
        }
        if(k!=n_keep[i])
          SID_trap_error("Particle count mismatch (ie. %d!=%d) during IDs read",ERROR_LOGIC,k,n_keep[i]);
      }
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
    SID_free((void **)&RNG);
    
    // Store everything in the data structure...

    //   ... particle counts ...
    for(i=0,n_particles_kept_all=0;i<N_GADGET_TYPE;i++){
      SID_Allreduce(&(n_of_type_rank[i]),&(n_of_type[i]),1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
      if(n_of_type[i]>0){
        n_particles_kept_all+=n_of_type[i];
        ADaPS_store(&(plist->data),(void *)(&(n_of_type_rank[i])),"n_%s",    ADaPS_SCALAR_SIZE_T,pname[i]);
        ADaPS_store(&(plist->data),(void *)(&(n_of_type[i])),     "n_all_%s",ADaPS_SCALAR_SIZE_T,pname[i]);
      }
    }
    ADaPS_store(&(plist->data),(void *)(&n_particles_all),"n_particles_all",ADaPS_SCALAR_SIZE_T);

    // Check that the right number of particles have been read
    if(n_particles_kept!=n_particles_rank)
      SID_log_warning("Rank %d did not receive the right number of particles (ie. %d!=%d)",ERROR_LOGIC,n_particles_kept,n_particles_rank);
    if(n_particles_kept_all!=n_particles_all)
      SID_log_warning("The right number of particles were not read (ie. %d!=%d)",ERROR_LOGIC,n_particles_kept_all,n_particles_all);

    //   ... positions ...
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type_rank[i]>0){
        ADaPS_store(&(plist->data),(void *)x_array[i],"x_%s",ADaPS_DEFAULT,pname[i]);
        ADaPS_store(&(plist->data),(void *)y_array[i],"y_%s",ADaPS_DEFAULT,pname[i]);
        ADaPS_store(&(plist->data),(void *)z_array[i],"z_%s",ADaPS_DEFAULT,pname[i]);
      }
    }
    /*
    //  ... velocities ...
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type_rank[i]>0){
        ADaPS_store(&(plist->data),(void *)vx_array[i],"vx_%s",ADaPS_DEFAULT,pname[i]);
        ADaPS_store(&(plist->data),(void *)vy_array[i],"vy_%s",ADaPS_DEFAULT,pname[i]);
        ADaPS_store(&(plist->data),(void *)vz_array[i],"vz_%s",ADaPS_DEFAULT,pname[i]);
      }
    }
    */
    //  ... ids ...
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type_rank[i]>0)
        ADaPS_store(&(plist->data),(void *)id_array[i],"id_%s",ADaPS_DEFAULT,pname[i]);
    }
    if(flag_LONGIDs)
      ADaPS_store(&(plist->data),(void *)&flag_LONGIDs,"flag_LONGIDs",ADaPS_SCALAR_INT);

    SID_log("Done.",SID_LOG_CLOSE);
  }
  else
    SID_trap_error("Could not find file with root {%s}",ERROR_IO_OPEN,filename_root_in);
}

