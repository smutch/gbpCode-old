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
                              int         i_load,
                              int         n_load,
                              slab_info  *slab,
                              cosmo_info *cosmo,
                              plist_info *plist);
void read_gadget_binary_local(char       *filename_root_in,
                              int         snapshot_number,
                              int         i_coord,
                              int         i_load,
                              int         n_load,
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
  gadget_header_info header;
  int                flag_filefound;
  int                flag_multifile;
  int                flag_file_type;
  flag_filefound=init_gadget_read(filename_root_in,snapshot_number,&flag_multifile,&flag_file_type,&header);

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
      n_all[i_type]     =(size_t)header.n_all[i_type];

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
    int n_read;
    if(n_load<n_files)
       n_read=n_files;
    else
       n_read=1;
    for(i_file=i_load;i_file<(i_load+n_read);i_file++){

      set_gadget_filename(filename_root_in,snapshot_number,i_file,flag_multifile,flag_file_type,filename);
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
                  if(pos_test>=slab->x_min_local && pos_test<slab->x_max_local)
                    n_of_type_local[i_type]++;
               }
            }
            i_step=MIN(READ_BUFFER_SIZE_LOCAL,header.n_file[i_type]-i_particle);
         }
      }
      if(n_files>1)
         SID_log("Done.",SID_LOG_CLOSE);
      fclose(fp_pos);
      fclose(fp_vel);
    }
    SID_log("Done.",SID_LOG_CLOSE);

    // Allocate arrays
    for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
       if(header.n_all[i_type]>0){
          x_array[i_type]=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_of_type_local[i_type]);
          y_array[i_type]=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_of_type_local[i_type]);
          z_array[i_type]=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_of_type_local[i_type]);
       }
    }

    // Perform read
    SID_log("Performing read...",SID_LOG_OPEN|SID_LOG_TIMER);
    for(i_file=i_load;i_file<(i_load+n_read);i_file++){

      set_gadget_filename(filename_root_in,snapshot_number,i_file,flag_multifile,flag_file_type,filename);
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
                     if(x_test<0)                x_test+=header.box_size;
                     if(x_test>=header.box_size) x_test-=header.box_size;
                     break;
                  case 2:
                     y_test+=(1e3*h_Hubble*((double)vel_buffer[index+1])/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
                     if(y_test<0)                y_test+=header.box_size;
                     if(y_test>=header.box_size) y_test-=header.box_size;
                     break;
                  case 3:
                     z_test+=(1e3*h_Hubble*((double)vel_buffer[index+2])/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
                     if(z_test<0)                z_test+=header.box_size;
                     if(z_test>=header.box_size) z_test-=header.box_size;
                     break;
               }
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
      n_particles_test +=header.n_all[i_type];
    }
    SID_Allreduce(&n_particles_local,&n_particles_read,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
    if(n_particles_read!=n_particles_test && n_load==1)
       SID_trap_error("Total particle counts don't make sense after read_gadget (ie. %zd!=%zd).",ERROR_LOGIC,n_particles_read,n_particles_test);
    for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
       SID_Allreduce(&(n_of_type_local[i_type]),&(n_of_type[i_type]),1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
       if(n_of_type[i_type]!=header.n_all[i_type] && n_load==1)
          SID_trap_error("Particle counts don't make sense after read_gadget (ie. %zd!=%zd).",ERROR_LOGIC,n_of_type[i_type],header.n_all[i_type]);
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
  int     n_load;
  int     n_used;
  int     flag_used[N_GADGET_TYPE];
  char    species_name[256];
  double  h_Hubble;
  double  n_spec;
  double  redshift;
  int     i_species;
  char    n_string[64];
  int             n[3];
  double          L[3];
  size_t          n_all;
  FILE           *fp_1D;
  FILE           *fp_2D;
  cosmo_info     *cosmo;
  field_info     *field[N_GADGET_TYPE];
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
  int     snapshot_number;
  int     i_compute;
  int     distribution_scheme;
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
  strcpy(filename_in_root,  argv[1]);
  snapshot_number=(int)atoi(argv[2]);
  strcpy(filename_out_root, argv[3]);
  grid_size      =(int)atoi(argv[4]);
  if(!strcmp(argv[5],"ngp") || !strcmp(argv[5],"NGP"))
     distribution_scheme=MAP2GRID_DIST_NGP;
  else if(!strcmp(argv[5],"cic") || !strcmp(argv[5],"CIC"))
     distribution_scheme=MAP2GRID_DIST_CIC;
  else if(!strcmp(argv[5],"tsc") || !strcmp(argv[5],"TSC"))
     distribution_scheme=MAP2GRID_DIST_TSC;
  else if(!strcmp(argv[5],"d12") || !strcmp(argv[5],"D12"))
     distribution_scheme=MAP2GRID_DIST_DWT12;
  else if(!strcmp(argv[5],"d20") || !strcmp(argv[5],"D20"))
     distribution_scheme=MAP2GRID_DIST_DWT20;
  else
     SID_trap_error("Invalid distribution scheme {%s} specified.",ERROR_SYNTAX,argv[5]);

  SID_log("Smoothing Gadget file {%s;snapshot=#%d} to a %dx%dx%d grid with %s kernel...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in_root,snapshot_number,grid_size,grid_size,grid_size,argv[5]);

  // Initialization -- default cosmology
  init_cosmo_std(&cosmo);

  // Initialization -- fetch header info
  double box_size;
  gadget_header_info header;
  int    flag_filefound;
  int    flag_multifile;
  int    flag_file_type;
  SID_log("Reading Gadget header...",SID_LOG_OPEN);
  flag_filefound=init_gadget_read(filename_in_root,snapshot_number,&flag_multifile,&flag_file_type,&header);
  if(flag_filefound){
     if(SID.I_am_Master){
        FILE *fp_in;
        char  filename[MAX_FILENAME_LENGTH];
        int   block_length_open;
        int   block_length_close;
        set_gadget_filename(filename_in_root,snapshot_number,0,flag_multifile,flag_file_type,filename);
        fp_in=fopen(filename,"r");
        fread(&block_length_open, sizeof(int),1,fp_in);
        fread(&header,            sizeof(gadget_header_info),1,fp_in);
        fread(&block_length_close,sizeof(int),1,fp_in);
        fclose(fp_in);
        if(block_length_open!=block_length_close)
           SID_trap_error("Block lengths don't match (ie. %d!=%d).",ERROR_LOGIC,block_length_open,block_length_close);
     }
     SID_Bcast(&header,sizeof(gadget_header_info),MASTER_RANK,SID.COMM_WORLD);
     redshift=header.redshift;
     h_Hubble=header.h_Hubble;
     box_size=header.box_size;
     if(SID.n_proc>1)
        n_load=1;
     else
        n_load=header.n_files;
     for(i_species=0,n_all=0,n_used=0;i_species<N_GADGET_TYPE;i_species++){
        n_all+=header.n_all[i_species];
        if(header.n_all[i_species]>0){
           n_used++;
           flag_used[i_species]=TRUE;
        }
        else
           flag_used[i_species]=FALSE;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Only process if there are >0 particles present
  if(n_used>0){

     // Loop over ithe real-space and 3 redshift-space frames
     int i_run;
     int n_run=4;
     n_run=1; // For now, hard-wire real-space calculation only
     for(i_run=0;i_run<n_run;i_run++){

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

        // Initialize the field that will hold the grid
        int        n[]={grid_size,grid_size,grid_size};
        double     L[]={box_size, box_size, box_size};
        int        i_init;
        for(i_species=0;i_species<N_GADGET_TYPE;i_species++){
           if(flag_used[i_species]){
              field[i_species]=(field_info *)SID_malloc(sizeof(field_info));
              init_field(3,n,L,field[i_species]);
              i_init=i_species;
           }
           else
              field[i_species]=NULL;
        }

        int i_load;
        for(i_load=0;i_load<n_load;i_load++){
           // Initialization -- read gadget file
           init_plist(&plist,&((field[i_init])->slab),GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
           char filename_root[MAX_FILENAME_LENGTH];
           read_gadget_binary_local(filename_in_root,
                                    snapshot_number,
                                    i_run,
                                    i_load,
                                    n_load,
                                    &(field[i_init]->slab),
                                    cosmo,
                                    &plist);

           // Generate power spectra
           for(i_species=0;i_species<plist.n_species;i_species++){

              // Determine how many particles of species i_species there are
              if(ADaPS_exist(plist.data,"n_all_%s",plist.species[i_species]))
                 n_all=((size_t *)ADaPS_fetch(plist.data,"n_all_%s",plist.species[i_species]))[0];
              else
                 n_all=0;
              // Compute power spectrum and write results
              if(n_all>0){
                 // Fetch the needed information
                 size_t   n_particles;
                 size_t   n_particles_local;
                 GBPREAL *x_particles_local;
                 GBPREAL *y_particles_local;
                 GBPREAL *z_particles_local;
                 GBPREAL *m_particles_local;
                 n_particles      =((size_t  *)ADaPS_fetch(plist.data,"n_all_%s",plist.species[i_species]))[0];
                 n_particles_local=((size_t  *)ADaPS_fetch(plist.data,"n_%s",    plist.species[i_species]))[0];
                 x_particles_local= (GBPREAL *)ADaPS_fetch(plist.data,"x_%s",    plist.species[i_species]);
                 y_particles_local= (GBPREAL *)ADaPS_fetch(plist.data,"y_%s",    plist.species[i_species]);
                 z_particles_local= (GBPREAL *)ADaPS_fetch(plist.data,"z_%s",    plist.species[i_species]);
                 if(ADaPS_exist(plist.data,"M_%s",plist.species[i_species]))
                    m_particles_local=(GBPREAL *)ADaPS_fetch(plist.data,"M_%s",plist.species[i_species]);
                 else
                    m_particles_local=NULL;

                 // Generate mass-field
                 if(i_load==0)
                    map_to_grid(n_particles_local,
                                x_particles_local,
                                y_particles_local,
                                z_particles_local,
                                m_particles_local,
                                cosmo,
                                redshift,
                                distribution_scheme,
                                1.,
                                field[i_species],
                                MAP2GRID_MODE_DEFAULT|MAP2GRID_MODE_NONORM);
                 else if(i_load==(n_load-1))
                    map_to_grid(n_particles_local,
                                x_particles_local,
                                y_particles_local,
                                z_particles_local,
                                m_particles_local,
                                cosmo,
                                redshift,
                                distribution_scheme,
                                1.,
                                field[i_species],
                                MAP2GRID_MODE_NOCLEAN);
                 else
                    map_to_grid(n_particles_local,
                                x_particles_local,
                                y_particles_local,
                                z_particles_local,
                                m_particles_local,
                                cosmo,
                                redshift,
                                distribution_scheme,
                                1.,
                                field[i_species],
                                MAP2GRID_MODE_NOCLEAN|MAP2GRID_MODE_NONORM);
              }
           }

           // Clean-up
           free_plist(&plist);
        }
        
        // Write results to disk
        char filename_out_species[MAX_FILENAME_LENGTH];
        init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
        for(i_species=0;i_species<plist.n_species;i_species++){
           if(flag_used[i_species]){
              SID_log("Writing results for the %s particles...",SID_LOG_OPEN,plist.species[i_species]);
              sprintf(filename_out_species,"%s_%s",filename_out_root,plist.species[i_species]);
              write_grid(field[i_species],filename_out_species,i_run,n_run,distribution_scheme,header.box_size);
              free_field(field[i_species]);
              SID_free(SID_FARG field[i_species]);
              SID_log("Done.",SID_LOG_CLOSE);
           }
        }
        free_plist(&plist);

        SID_log("Done.",SID_LOG_CLOSE);
     }
  } 

  // Clean-up
  SID_log("Clean-up...",SID_LOG_OPEN);
  free_cosmo(&cosmo);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}


