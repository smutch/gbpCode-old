#define  _MAIN
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpCosmo.h>
#include <gbpClustering.h>

#define _FILE_OFFSET_BITS 64

#define READ_BUFFER_SIZE_LOCAL    (1024*1024)
#define READ_BUFFER_ALLOC_LOCAL 3*(1024*1024)

void read_gadget_binary_local(char       *filename_root_in,
                              int         snapshot_number,
                              int         i_coord,
                              slab_info  *slab,
                              cosmo_info *cosmo,
                              plist_info *plist);
void read_gadget_binary_local(char       *filename_root_in,
                              int         snapshot_number,
                              int         i_coord,
                              slab_info  *slab,
                              cosmo_info *cosmo,
                              plist_info *plist){
  size_t     n_of_type_local[N_GADGET_TYPE];
  size_t     n_of_type[N_GADGET_TYPE];
  size_t     type_counter[N_GADGET_TYPE];
  GBPREAL   *x_array[N_GADGET_TYPE];
  GBPREAL   *y_array[N_GADGET_TYPE];
  GBPREAL   *z_array[N_GADGET_TYPE];
  int        i_type;

  // Determine file format and read the header
  gadget_read_info   fp_gadget;
  int                flag_filefound=init_gadget_read(filename_root_in,snapshot_number,&fp_gadget);
  int                flag_multifile=fp_gadget.flag_multifile;
  int                flag_file_type=fp_gadget.flag_file_type;
  gadget_header_info header        =fp_gadget.header;

  // A file was found ... 
  if(flag_filefound){
    char **pname;
    SID_log("Reading GADGET binary file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_root_in);

    pname=plist->species;

    // Expansion factor (or time) 
    ADaPS_store(&(plist->data),(void *)(&(header.time)),"expansion_factor",ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&(plist->data),(void *)(&(header.time)),"time",            ADaPS_SCALAR_DOUBLE);

    // Redshift
    double d_value;
    d_value=(double)header.redshift;
    ADaPS_store(&(plist->data),(void *)(&d_value),"redshift",ADaPS_SCALAR_DOUBLE);

    // Number of particles for each species in all files
    size_t n_all[N_GADGET_TYPE];
    for(i_type=0;i_type<N_GADGET_TYPE;i_type++)
       n_all[i_type]=(size_t)header.n_all_lo_word[i_type]+((size_t)header.n_all_hi_word[i_type])<<32;

    // Number of files in this snapshot 
    int n_files;
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
    double h_Hubble;
    double redshift;
    h_Hubble=(double)header.h_Hubble; 
    if(h_Hubble<1e-10) h_Hubble=1.;
    ADaPS_store(&(plist->data),(void *)(&h_Hubble),"h_Hubble",ADaPS_SCALAR_DOUBLE);
    redshift=header.redshift;
    ADaPS_store(&(plist->data),(void *)(&redshift),"redshift",ADaPS_SCALAR_DOUBLE);

    // Count and report the total number of particles
    size_t n_particles_all;
    int    n_non_zero;
    n_particles_all=0;
    for(i_type=0,n_non_zero=0;i_type<N_GADGET_TYPE;i_type++){
      if(n_all[i_type]>0){
        n_particles_all+=n_all[i_type];
        n_non_zero++;
      }
    }
    SID_log("%zd",SID_LOG_CONTINUE,n_particles_all);
    if(n_non_zero>0)
      SID_log(" (",SID_LOG_CONTINUE,n_particles_all);
    for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
      if(n_all[i_type]>0){
        if(i_type==n_non_zero-1){
          if(n_non_zero>1)
            SID_log("and %lld %s",SID_LOG_CONTINUE,n_all[i_type],pname[i_type]);
          else
            SID_log("%lld %s",SID_LOG_CONTINUE,n_all[i_type],pname[i_type]);
        }
        else{
          if(n_non_zero>1)
            SID_log("%lld %s, ",SID_LOG_CONTINUE,n_all[i_type],pname[i_type]);        
          else
            SID_log("%lld %s",SID_LOG_CONTINUE,n_all[i_type],pname[i_type]);        
        }
      }
    }
    if(n_non_zero>0)
      SID_log(") particles...",SID_LOG_CONTINUE);
    else
      SID_log(" particles...",SID_LOG_CONTINUE);

    // Count the number of particles that will be scattered to each rank
    char     filename[MAX_FILENAME_LENGTH];
    size_t   k_particle;
    int      i_file;
    int      record_length_open;
    int      record_length_close;
    size_t   i_particle;
    size_t   i_buffer;
    size_t   i_step;
    int      i_type;
    size_t   index;
    GBPREAL *pos_buffer;
    GBPREAL *vel_buffer;
    double   pos_test;

    // Initialize some arrays
    pos_buffer=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*READ_BUFFER_ALLOC_LOCAL);
    vel_buffer=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*READ_BUFFER_ALLOC_LOCAL);
    for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
       n_of_type_local[i_type]=0;
       n_of_type[i_type]      =0;
       type_counter[i_type]   =0;
    }

    // Determine how many particles of each type will end-up on each core
    SID_log("Performing domain decomposition...",SID_LOG_OPEN|SID_LOG_TIMER);
    for(i_file=0;i_file<n_files;i_file++){

      set_gadget_filename(&fp_gadget,i_file,filename);
      if(n_files>1)
         SID_log("Processing file #%d of %d...",SID_LOG_OPEN,i_file+1,n_files);

      // Read header and move to the positions
      FILE *fp_pos;
      FILE *fp_vel;
      fp_pos=fopen(filename,"r");
      fread(&record_length_open,4,1,fp_pos);
      fread(&header,sizeof(gadget_header_info),1,fp_pos);
      fread(&record_length_close,4,1,fp_pos);
      if(record_length_open!=record_length_close)
        SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
      fread(&record_length_open,4,1,fp_pos);

      // Create a file pointer to the velocities
      fp_vel=fopen(filename,"r");
      fread(&record_length_open,4,1,fp_vel);
      fseeko(fp_vel,(off_t)(record_length_open),SEEK_CUR);
      fread(&record_length_close,4,1,fp_vel);
      fread(&record_length_open,4,1,fp_vel);
      fseeko(fp_vel,(off_t)(record_length_open),SEEK_CUR);
      fread(&record_length_close,4,1,fp_vel);
      if(record_length_open!=record_length_close)
        SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
      fread(&record_length_open,4,1,fp_vel);

      // We only have to worry about z-space effects for domain decomposition in this one case.
      if(i_coord==1){
         for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
            for(i_particle=0;i_particle<header.n_file[i_type];i_particle+=i_step){
               i_step=MIN(READ_BUFFER_SIZE_LOCAL,header.n_file[i_type]-i_particle);
               if(SID.I_am_Master){
                  fread(pos_buffer,sizeof(GBPREAL),3*i_step,fp_pos);
                  fread(vel_buffer,sizeof(GBPREAL),3*i_step,fp_vel);
               }
               SID_Bcast(pos_buffer,sizeof(GBPREAL)*3*i_step,MASTER_RANK,SID.COMM_WORLD);
               SID_Bcast(vel_buffer,sizeof(GBPREAL)*3*i_step,MASTER_RANK,SID.COMM_WORLD);
               for(i_buffer=0;i_buffer<i_step;i_buffer++){
                  index=3*i_buffer;
                  pos_test =(double)(pos_buffer[index]);
                  pos_test+=(double)(1e3*h_Hubble*((double)vel_buffer[index])/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
                  if(pos_test<0)                pos_test+=header.box_size;
                  if(pos_test>=header.box_size) pos_test-=header.box_size;
                  if(pos_test>=slab->x_min_local && pos_test<slab->x_max_local)
                     n_of_type_local[i_type]++;
               }
            }
         }
      }
      else{
         for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
            for(i_particle=0;i_particle<header.n_file[i_type];i_particle+=i_step){
               i_step=MIN(READ_BUFFER_SIZE_LOCAL,header.n_file[i_type]-i_particle);
               if(SID.I_am_Master)
                  fread(pos_buffer,sizeof(GBPREAL),3*i_step,fp_pos);
               SID_Bcast(pos_buffer,sizeof(GBPREAL)*3*i_step,MASTER_RANK,SID.COMM_WORLD);
               for(i_buffer=0;i_buffer<i_step;i_buffer++){
                  pos_test=pos_buffer[3*i_buffer];
                  if(pos_test<0)                pos_test+=header.box_size;
                  if(pos_test>=header.box_size) pos_test-=header.box_size;
                  if(pos_test>=slab->x_min_local && pos_test<slab->x_max_local)
                    n_of_type_local[i_type]++;
               }
            }
         }
      }
      if(n_files>1)
         SID_log("Done.",SID_LOG_CLOSE);
      fclose(fp_pos);
      fclose(fp_vel);
    }
    size_t n_local;
    for(i_type=0,n_local=0;i_type<N_GADGET_TYPE;i_type++) 
       n_local+=n_of_type_local[i_type];
    if(SID.n_proc>1){
       SID_log("Results:",SID_LOG_OPEN);
       int i_rank;
       for(i_rank=0;i_rank<SID.n_proc;i_rank++){
          size_t n_report;
          n_report=n_local;
          SID_Bcast(&n_report,sizeof(size_t),i_rank,SID.COMM_WORLD);
          SID_log("rank #%04d: n_particles=%zd",SID_LOG_COMMENT,i_rank,n_report);
       }
       SID_log("Done.",SID_LOG_SILENT_CLOSE);
    }
    SID_log("Done.",SID_LOG_CLOSE);

    // Allocate arrays
    for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
       if(n_all[i_type]>0){
          x_array[i_type]=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_of_type_local[i_type]);
          y_array[i_type]=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_of_type_local[i_type]);
          z_array[i_type]=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_of_type_local[i_type]);
       }
    }

    // Perform read
    SID_log("Performing read...",SID_LOG_OPEN|SID_LOG_TIMER);
    for(i_file=0;i_file<n_files;i_file++){

      set_gadget_filename(&fp_gadget,i_file,filename);
      if(n_files>1)
         SID_log("Processing file #%d of %d...",SID_LOG_OPEN,i_file+1,n_files);

      // Read header and move to the positions
      FILE *fp_pos;
      FILE *fp_vel;
      fp_pos=fopen(filename,"r");
      fread(&record_length_open,4,1,fp_pos);
      fread(&header,sizeof(gadget_header_info),1,fp_pos);
      fread(&record_length_close,4,1,fp_pos);
      if(record_length_open!=record_length_close)
        SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
      fread(&record_length_open,4,1,fp_pos);

      // Create a file pointer to the velocities
      fp_vel=fopen(filename,"r");
      fread(&record_length_open,4,1,fp_vel);
      fseeko(fp_vel,(off_t)(record_length_open),SEEK_CUR);
      fread(&record_length_close,4,1,fp_vel);
      fread(&record_length_open,4,1,fp_vel);
      fseeko(fp_vel,(off_t)(record_length_open),SEEK_CUR);
      fread(&record_length_close,4,1,fp_vel);
      if(record_length_open!=record_length_close)
        SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
      fread(&record_length_open,4,1,fp_vel);

      // Perform the read and populate the local position arrays
      size_t   i_particle;
      size_t   i_step;
      int      i_type;
      for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
         for(i_particle=0;i_particle<header.n_file[i_type];i_particle+=i_step){
            i_step=MIN(READ_BUFFER_SIZE_LOCAL,header.n_file[i_type]-i_particle);
            if(SID.I_am_Master){
               fread(pos_buffer,sizeof(GBPREAL),3*i_step,fp_pos);
               if(i_coord>0)
                  fread(vel_buffer,sizeof(GBPREAL),3*i_step,fp_vel);
            }
            SID_Bcast(pos_buffer,sizeof(GBPREAL)*3*i_step,MASTER_RANK,SID.COMM_WORLD);
            if(i_coord>0)
               SID_Bcast(vel_buffer,sizeof(GBPREAL)*3*i_step,MASTER_RANK,SID.COMM_WORLD);
            for(i_buffer=0;i_buffer<i_step;i_buffer++){
               double x_test;
               double y_test;
               double z_test;
               index=3*i_buffer;
               x_test=pos_buffer[index+0];
               y_test=pos_buffer[index+1];
               z_test=pos_buffer[index+2];
               switch(i_coord){
                  case 1:
                     x_test+=(1e3*h_Hubble*((double)vel_buffer[index+0])/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
                     break;
                  case 2:
                     y_test+=(1e3*h_Hubble*((double)vel_buffer[index+1])/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
                     break;
                  case 3:
                     z_test+=(1e3*h_Hubble*((double)vel_buffer[index+2])/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
                     break;
               }
               if(x_test<0)                x_test+=header.box_size;
               if(x_test>=header.box_size) x_test-=header.box_size;
               if(y_test<0)                y_test+=header.box_size;
               if(y_test>=header.box_size) y_test-=header.box_size;
               if(z_test<0)                z_test+=header.box_size;
               if(z_test>=header.box_size) z_test-=header.box_size;
               if(x_test>=slab->x_min_local && x_test<slab->x_max_local){
                  x_array[i_type][type_counter[i_type]]=x_test;
                  y_array[i_type][type_counter[i_type]]=y_test;
                  z_array[i_type][type_counter[i_type]]=z_test;
                  type_counter[i_type]++;
               }
            }
         }
      }

      // Close file pointers
      fclose(fp_pos);
      fclose(fp_vel);

      if(n_files>1)
         SID_log("Done.",SID_LOG_CLOSE);
    }
    SID_free(SID_FARG pos_buffer);
    SID_free(SID_FARG vel_buffer);
    SID_log("Done.",SID_LOG_CLOSE);

    // Sanity checks
    size_t n_particles_local;
    size_t n_particles_read;
    size_t n_particles_test;
    for(i_type=0,n_particles_local=0,n_particles_test=0;i_type<N_GADGET_TYPE;i_type++){
      n_particles_local+=n_of_type_local[i_type];
      n_particles_test +=n_all[i_type];
    }
    SID_Allreduce(&n_particles_local,&n_particles_read,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
    if(n_particles_read!=n_particles_test)
       SID_trap_error("Total particle counts don't make sense after read_gadget (ie. %zd!=%zd).",ERROR_LOGIC,n_particles_read,n_particles_test);
    for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
       SID_Allreduce(&(n_of_type_local[i_type]),&(n_of_type[i_type]),1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
       if(n_of_type[i_type]!=n_all[i_type])
          SID_trap_error("Particle counts don't make sense after read_gadget (ie. %zd!=%zd).",ERROR_LOGIC,n_of_type[i_type],n_all[i_type]);
    }

    // Store results
    for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
      if(n_of_type[i_type]>0){
        ADaPS_store(&(plist->data),(void *)(&(n_of_type_local[i_type])),"n_%s",    ADaPS_SCALAR_SIZE_T,pname[i_type]);
        ADaPS_store(&(plist->data),(void *)(&(n_of_type[i_type])),      "n_all_%s",ADaPS_SCALAR_SIZE_T,pname[i_type]);
      }
    }
    ADaPS_store(&(plist->data),(void *)(&n_particles_all),"n_particles_all",ADaPS_SCALAR_SIZE_T);
    for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
      if(n_of_type_local[i_type]>0){
        ADaPS_store(&(plist->data),(void *)x_array[i_type],"x_%s",ADaPS_DEFAULT,pname[i_type]);
        ADaPS_store(&(plist->data),(void *)y_array[i_type],"y_%s",ADaPS_DEFAULT,pname[i_type]);
        ADaPS_store(&(plist->data),(void *)z_array[i_type],"z_%s",ADaPS_DEFAULT,pname[i_type]);
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
}

int main(int argc, char *argv[]){
  int     n_species;
  char    species_name[256];
  int     i_species;
  char    n_string[64];
  int             n[3];
  double          L[3];
  FILE           *fp_1D;
  FILE           *fp_2D;
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

  // Parse arguments
  int grid_size;
  char filename_in_root[MAX_FILENAME_LENGTH];
  char filename_out_root[MAX_FILENAME_LENGTH];
  char cosmo_name[16];
  strcpy(filename_in_root,  argv[1]);
  snapshot_number=(int)atoi(argv[2]);
  strcpy(filename_out_root, argv[3]);
  grid_size      =(int)atoi(argv[4]);
  strcpy(cosmo_name,        argv[5]);
  SID_log("Processing the power spectra of {%s}, snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in_root,snapshot_number);

  // Initialization -- fetch header info
  SID_log("Reading Gadget header...",SID_LOG_OPEN);
  gadget_read_info   fp_gadget;
  int                flag_filefound=init_gadget_read(filename_in_root,snapshot_number,&fp_gadget);
  int                flag_multifile=fp_gadget.flag_multifile;
  int                flag_file_type=fp_gadget.flag_file_type;
  gadget_header_info header        =fp_gadget.header;

  // Set the cosmology
  cosmo_info *cosmo=NULL;
  double box_size        =header.box_size;
  double h_Hubble        =header.h_Hubble;
  double redshift        =header.redshift;
  double expansion_factor=header.time;
/*
  double Omega_M         =header.Omega_M;
  double Omega_Lambda    =header.Omega_Lambda;
  double Omega_k         =1.-Omega_Lambda-Omega_M;
  double Omega_b         =0.; // not needed, so doesn't matter
  double f_gas           =Omega_b/Omega_M;
  double sigma_8         =0.; // not needed, so doesn't matter
  double n_spec          =0.; // not needed, so doesn't matter
  init_cosmo(&cosmo,
             cosmo_name,
             Omega_Lambda,
             Omega_M,
             Omega_k,
             Omega_b,
             f_gas,
             h_Hubble,
             sigma_8,
             n_spec);
*/
  read_gbpCosmo_file(&cosmo,cosmo_name);
  SID_log("Done.",SID_LOG_CLOSE);

  // Set the k ranges
  double k_Nyq;
  k_Nyq   =(TWO_PI*(double)grid_size/box_size)/2.;
  k_min_1D=0.02;
  dk_1D   =0.02;
  k_max_1D=dk_1D*(float)((int)(k_Nyq/dk_1D));
  k_min_2D=0.0;
  dk_2D   =0.02;
  k_max_2D=dk_2D*(float)((int)(k_Nyq/dk_2D));

  // Initialize the power spectrum
  pspec_info *pspec;
  pspec=(pspec_info *)SID_malloc(sizeof(pspec_info)*N_GADGET_TYPE);
  size_t  n_total;
  size_t  n_all[N_GADGET_TYPE];
  for(i_species=0,n_total=0;i_species<N_GADGET_TYPE;i_species++){
     n_all[i_species]=(size_t)header.n_all_lo_word[i_species]+((size_t)header.n_all_hi_word[i_species])<<32;
     if(i_species>0) SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
     init_pspec(&(pspec[i_species]),NULL,cosmo,
                MAP2GRID_DIST_DWT20,
                redshift,box_size,grid_size,
                k_min_1D,k_max_1D,dk_1D,
                k_min_2D,k_max_2D,dk_2D);
     if(i_species>0) SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
     n_total+=n_all[i_species];
  }

  // Only process a species if there are >0 particles present
  if(n_total>0){

     // Loop over ithe real-space and 3 redshift-space frames
     int i_run;
     for(i_run=0;i_run<4;i_run++){

        // Read catalog
        switch(i_run){
        case 0:
           SID_log("Processing real-space ...",SID_LOG_OPEN|SID_LOG_TIMER);
           break;
        case 1:
           SID_log("Processing v_x redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
           break;
        case 2:
           SID_log("Processing v_y redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
           break;
        case 3:
           SID_log("Processing v_z redsift space...",SID_LOG_OPEN|SID_LOG_TIMER);
           break;
        }

        // Initialization -- data structure which holds all   
        //                   the (local) particle information  
        init_plist(&plist,&(pspec[0].FFT.slab),GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);

        // Initialization -- read gadget file
        char filename_root[MAX_FILENAME_LENGTH];
        read_gadget_binary_local(filename_in_root,
                                 snapshot_number,
                                 i_run,
                                 &(pspec[0].FFT.slab),
                                 cosmo,
                                 &plist);

        // Generate power spectra
        for(i_species=0;i_species<plist.n_species;i_species++){

           // Compute power spectrum and write results
           if(n_all[i_species]>0){
              compute_pspec(&plist,plist.species[i_species],&(pspec[i_species]),i_run);
              char filename_out_species[MAX_FILENAME_LENGTH];
              SID_log("Writing results for the %s particles...",SID_LOG_OPEN,plist.species[i_species]);
              sprintf(filename_out_species,"%s_%s",filename_out_root,plist.species[i_species]);
              write_pspec(&(pspec[i_species]),filename_out_species,&plist,plist.species[i_species]);
              SID_log("Done.",SID_LOG_CLOSE);
           }
        }

        // Clean-up
        free_plist(&plist);

        SID_log("Done.",SID_LOG_CLOSE);
     }
  } 

  // Clean-up
  SID_log("Clean-up...",SID_LOG_OPEN);
  free_cosmo(&cosmo);
  for(i_species=0;i_species<N_GADGET_TYPE;i_species++)
     free_pspec(&(pspec[i_species]));
  SID_free(SID_FARG pspec);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}


