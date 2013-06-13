#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHalos.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>

#define GADGET_BUFFER_SIZE_LOCAL  128*SIZE_OF_MEGABYTE

#define _FILE_OFFSET_BITS 64

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
  size_t    i_list_offset[N_GADGET_TYPE];
  size_t    n_particles_file;
  size_t    n_particles_kept;
  size_t    n_of_type_rank[N_GADGET_TYPE];
  size_t    n_of_type[N_GADGET_TYPE];
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
  unsigned int       record_length_open;
  unsigned int       record_length_close;
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
  FILE     *fp;
  void     *buffer;
  size_t   *buffer_index=NULL;
  void     *buffer_i;
  size_t    id_search;
  size_t    id_value;
  size_t    n_buffer;
  int       n_non_zero;
  int       seed=1073743;
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
  
  // Determine file format
  n_files=1;
  for(i_file=0;i_file<3 && !flag_filefound;i_file++){  
    if(i_file==0)
      sprintf(filename,"%s/snapshot_%03d/snapshot_%03d",filename_root_in,snapshot_number,snapshot_number);
    else if(i_file==1)
      sprintf(filename,"%s/snapshot_%03d",filename_root_in,snapshot_number);
    else if(i_file==2)
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
  SID_log("Reading GADGET binary file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_root_in);
  if(flag_filefound){

    pname=plist->species;

    // Header record length
    SID_log("Reading header...",SID_LOG_OPEN);
    fread(&record_length_open,4,1,fp);
    if(record_length_open!=GADGET_HEADER_SIZE)
      SID_log_warning("Problem with GADGET record size (opening size of header is wrong)",ERROR_LOGIC);

    // Read header
    fread(&header,sizeof(gadget_header_info),1,fp); 

    // Number of particles for each species in this file
    for(i=0;i<N_GADGET_TYPE;i++)
      n_of_type[i]=(size_t)header.n_file[i];

    // Expansion factor (or time) 
    expansion_factor=header.time;
    ADaPS_store(&(plist->data),(void *)(&expansion_factor),"expansion_factor",ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&(plist->data),(void *)(&(header.time)),   "time",            ADaPS_SCALAR_DOUBLE);

    // Redshift
    d_value=(double)header.redshift;
    ADaPS_store(&(plist->data),(void *)(&d_value),"redshift",ADaPS_SCALAR_DOUBLE);

    // Number of particles for each species in all files
    for(i=0;i<N_GADGET_TYPE;i++)
      n_all[i]=(size_t)header.n_all[i];

    // Number of files in this snapshot 
    ADaPS_store(&(plist->data),(void *)(&(header.n_files)),"n_files",ADaPS_SCALAR_INT);
    if(flag_multifile)
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
      mass_array[i]=header.mass_array[i];
      if(header.n_all_high_word[i]>0)
        n_all[i]+=(((size_t)(header.n_all_high_word[i])) << 32);
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

    // Fetch (and sort by id) local group particles
    SID_log("Sorting group particle ids...",SID_LOG_OPEN);
    n_particles_rank         =((size_t *)ADaPS_fetch(plist->data,"n_particles_%03d",    snapshot_number))[0];
    n_particles_all_in_groups=((size_t *)ADaPS_fetch(plist->data,"n_particles_all_%03d",snapshot_number))[0];
    id_list                  = (size_t *)ADaPS_fetch(plist->data,"particle_ids_%03d",   snapshot_number);
    merge_sort((void *)id_list,(size_t)n_particles_rank,&id_list_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
    for(i=0;i<N_GADGET_TYPE;i++)
      n_of_type_rank[i]=0;
    n_of_type[GADGET_TYPE_DARK]     =n_particles_all_in_groups;
    n_of_type_rank[GADGET_TYPE_DARK]=n_particles_rank;
    SID_log("Done.",SID_LOG_CLOSE);

    // Check integretry of particle ids
    SID_log("Checking integrity of group particle ids...",SID_LOG_OPEN);  
    for(j=1,n_warning=0;j<n_particles_rank && n_warning<10;j++){
      if(id_list[id_list_index[j]]==id_list[id_list_index[j-1]] && n_warning<100){
        fprintf(stderr,"WARNING: group particle id=%zd is duplicated on rank=%d!\n",id_list[id_list_index[j]],SID.My_rank);
        n_warning++;
      }
    }
    if(n_warning>0)
      SID_trap_error("Error in group particle ids!",ERROR_LOGIC);
    SID_log("Done.",SID_LOG_CLOSE);

    // Allocate data arrays
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type_rank[i]>0){
        id_array[i]=(size_t    *)SID_malloc(sizeof(size_t) *(size_t)n_of_type_rank[i]);
        x_array[i] =(GBPREAL   *)SID_malloc(sizeof(GBPREAL)*(size_t)n_of_type_rank[i]);
        y_array[i] =(GBPREAL   *)SID_malloc(sizeof(GBPREAL)*(size_t)n_of_type_rank[i]);
        z_array[i] =(GBPREAL   *)SID_malloc(sizeof(GBPREAL)*(size_t)n_of_type_rank[i]);
        vx_array[i]=(GBPREAL   *)SID_malloc(sizeof(GBPREAL)*(size_t)n_of_type_rank[i]);
        vy_array[i]=(GBPREAL   *)SID_malloc(sizeof(GBPREAL)*(size_t)n_of_type_rank[i]);
        vz_array[i]=(GBPREAL   *)SID_malloc(sizeof(GBPREAL)*(size_t)n_of_type_rank[i]);
      }
    }
    
    // Set offsets
    for(i=0;i<N_GADGET_TYPE;i++)
      i_list_offset[i]=0;

    // Read data
    for(i_file=0,n_particles_kept=0;i_file<n_files;i_file++){
      if(n_files>1)
        SID_log("Reading file %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file+1,n_files);
      else
        SID_log("Performing read...",SID_LOG_OPEN|SID_LOG_TIMER);

      if(flag_file_type==0)
        sprintf(filename,"%s/snapshot_%03d/snapshot_%03d",filename_root_in,snapshot_number,snapshot_number);
      else if(flag_file_type==1)
        sprintf(filename,"%s/snapshot_%03d",filename_root_in,snapshot_number);
      else if(flag_file_type==2)
        sprintf(filename,"%s_%03d",filename_root_in,snapshot_number);
      if(flag_multifile)
        sprintf(filename,"%s.%d",filename,i_file);
        
      // Read header and set positions pointer
      FILE *fp_positions;
      if(SID.I_am_Master){
        fp_positions=fopen(filename,"r");
        fread(&record_length_open, sizeof(int),1,               fp_positions);
        fread(&header,             sizeof(gadget_header_info),1,fp_positions);
        fread(&record_length_close,sizeof(int),1,               fp_positions);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
        fread(&record_length_open,sizeof(int),1,fp_positions);
      }
      SID_Bcast(&header,sizeof(gadget_header_info),MASTER_RANK,SID.COMM_WORLD);
      for(i=0,n_particles_file=0;i<N_GADGET_TYPE;i++)
        n_particles_file+=(size_t)header.n_file[i];

      // Set velocities pointer
      FILE *fp_velocities;
      if(SID.I_am_Master){
        fp_velocities=fopen(filename,"r");
        fread(&record_length_open,sizeof(int),1,fp_velocities);
        fseeko(fp_velocities,(off_t)record_length_open,SEEK_CUR);
        fread(&record_length_close,sizeof(int),1,fp_velocities);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
        fread(&record_length_open,sizeof(int),1,fp_velocities);
        fseeko(fp_velocities,(off_t)record_length_open,SEEK_CUR);
        fread(&record_length_close,sizeof(int),1,fp_velocities);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of positions)",ERROR_LOGIC);
        fread(&record_length_open,sizeof(int),1,fp_velocities);
      }

      // Set IDs pointer
      FILE *fp_IDs;
      if(SID.I_am_Master){
         fp_IDs=fopen(filename,"r");
         fread(&record_length_open,sizeof(int),1,fp_IDs);
         fseeko(fp_IDs,(off_t)record_length_open,SEEK_CUR);
         fread(&record_length_close,sizeof(int),1,fp_IDs);
         if(record_length_open!=record_length_close)
           SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
         fread(&record_length_open,sizeof(int),1,fp_IDs);
         fseeko(fp_IDs,(off_t)record_length_open,SEEK_CUR);
         fread(&record_length_close,sizeof(int),1,fp_IDs);
         if(record_length_open!=record_length_close)
           SID_log_warning("Problem with GADGET record size (close of positions)",ERROR_LOGIC);
         fread(&record_length_open,sizeof(int),1,fp_IDs);
         fseeko(fp_IDs,(off_t)record_length_open,SEEK_CUR);
         fread(&record_length_close,sizeof(int),1,fp_IDs);
         if(record_length_open!=record_length_close)
           SID_log_warning("Problem with GADGET record size (close of velocities)",ERROR_LOGIC);
         fread(&record_length_open,sizeof(int),1,fp_IDs);
         if((size_t)record_length_open/n_particles_file==sizeof(long long)){
            SID_log("(long long IDs)...",SID_LOG_CONTINUE);
            flag_LONGIDS=TRUE;
         }
         else if((size_t)record_length_open/n_particles_file==sizeof(int)){
            SID_log("(int IDs)...",SID_LOG_CONTINUE);
            flag_LONGIDS=FALSE;
         }
         else
            SID_trap_error("IDs record length (%d) does not set a sensible id byte size for the number of particles given in the header (%s)",ERROR_LOGIC,record_length_open,n_particles_file);
      }
      SID_Bcast(&flag_LONGIDS,sizeof(int),MASTER_RANK,SID.COMM_WORLD);

      // Allocate buffers
      int      n_buffer_max=MIN(n_particles_file,4*1024*1024);
      GBPREAL *buffer_positions ;
      GBPREAL *buffer_velocities;
      size_t  *buffer_IDs       ;
      buffer_positions =(GBPREAL *)SID_malloc(3*n_buffer_max*sizeof(GBPREAL));
      buffer_velocities=(GBPREAL *)SID_malloc(3*n_buffer_max*sizeof(GBPREAL));
      buffer_IDs       =(size_t  *)SID_malloc(  n_buffer_max*sizeof(size_t));

      // Perform read
      int     i_buffer;
      int     n_buffer_left;
      int     n_buffer;
      int     i_species;
      int     i_particle;
      int     j_particle;
      int     i_list;
      size_t *buffer_index=NULL;
      for(i_species=0,i_buffer=n_buffer_max,i_list=0,n_buffer_left=n_particles_file;i_species<N_GADGET_TYPE;i_species++){
         n_keep[i_species]=0;
         if(header.n_file[i_species]>0){
            for(i_particle=0;i_particle<header.n_file[i_species];i_particle++){
               // Perform a buffered read
               if(i_buffer>=n_buffer_max){
                  n_buffer=MIN(n_buffer_max,n_buffer_left);
                  if(SID.I_am_Master){
                     fread(buffer_positions, sizeof(GBPREAL),3*n_buffer,fp_positions);
                     fread(buffer_velocities,sizeof(GBPREAL),3*n_buffer,fp_velocities);
                     if(!flag_LONGIDS){
                        int *buffer_IDs_int=(int *)buffer_IDs;
                        fread(buffer_IDs_int,sizeof(int),n_buffer,fp_IDs);
                        for(j_particle=n_buffer-1;j_particle>=0;j_particle--)
                           buffer_IDs[j_particle]=(size_t)buffer_IDs_int[j_particle];
                     }
                     else
                        fread(buffer_IDs,sizeof(size_t),n_buffer,fp_IDs);
                  }
                  SID_Bcast(buffer_positions, 3*n_buffer*sizeof(GBPREAL),MASTER_RANK,SID.COMM_WORLD);
                  SID_Bcast(buffer_velocities,3*n_buffer*sizeof(GBPREAL),MASTER_RANK,SID.COMM_WORLD);
                  SID_Bcast(buffer_IDs,         n_buffer*sizeof(size_t), MASTER_RANK,SID.COMM_WORLD);
                  SID_free(SID_FARG buffer_index);
                  merge_sort(buffer_IDs,(size_t)n_buffer,&buffer_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
                  n_buffer_left-=n_buffer;
                  i_buffer=0;
                  i_list  =0;
               }

               if(i_list<n_particles_rank){
                  // Move to the next particle
                  size_t id_test;
                  size_t j_buffer=buffer_index[i_buffer];
                  id_test=buffer_IDs[j_buffer];
                  while(id_list[id_list_index[i_list]]<id_test && i_list<(n_particles_rank-1))
                    i_list++;

                  // If we've found a particle, add it to the arrays
                  if(id_list[id_list_index[i_list]]==id_test){
                    GBPREAL f1_value;
                    GBPREAL f2_value;
                    GBPREAL f3_value;
                    GBPREAL f4_value;
                    GBPREAL f5_value;
                    GBPREAL f6_value;
                    f1_value=((GBPREAL *)buffer_positions)[3*j_buffer+0];
                    f2_value=((GBPREAL *)buffer_positions)[3*j_buffer+1];
                    f3_value=((GBPREAL *)buffer_positions)[3*j_buffer+2];
                    f4_value=((GBPREAL *)buffer_velocities)[3*j_buffer+0];
                    f5_value=((GBPREAL *)buffer_velocities)[3*j_buffer+1];
                    f6_value=((GBPREAL *)buffer_velocities)[3*j_buffer+2];
                    id_array[i_species][n_keep[i_species]+i_list_offset[i_species]]=id_test;
                    x_array[i_species][n_keep[i_species]+i_list_offset[i_species]] =f1_value*(GBPREAL)(plist->length_unit/h_Hubble);
                    y_array[i_species][n_keep[i_species]+i_list_offset[i_species]] =f2_value*(GBPREAL)(plist->length_unit/h_Hubble);
                    z_array[i_species][n_keep[i_species]+i_list_offset[i_species]] =f3_value*(GBPREAL)(plist->length_unit/h_Hubble);
                    vx_array[i_species][n_keep[i_species]+i_list_offset[i_species]]=f4_value*(GBPREAL)(plist->velocity_unit*sqrt(expansion_factor));
                    vy_array[i_species][n_keep[i_species]+i_list_offset[i_species]]=f5_value*(GBPREAL)(plist->velocity_unit*sqrt(expansion_factor));
                    vz_array[i_species][n_keep[i_species]+i_list_offset[i_species]]=f6_value*(GBPREAL)(plist->velocity_unit*sqrt(expansion_factor));
                    n_keep[i_species]++;
                    n_particles_kept++;
                    i_list++;
                  }
               }
               i_buffer++;
            } // i_particle
         } // if n>0
      } // i_species

      // Update offsets
      for(i_species=0;i_species<N_GADGET_TYPE;i_species++)
        i_list_offset[i_species]+=n_keep[i_species];

      // Clean-up
      SID_free(SID_FARG buffer_positions);
      SID_free(SID_FARG buffer_velocities);
      SID_free(SID_FARG buffer_IDs);
      SID_free(SID_FARG buffer_index);
      if(SID.I_am_Master){
         fclose(fp_positions);
         fclose(fp_velocities);
         fclose(fp_IDs);
      }
      SID_log("Done.",SID_LOG_CLOSE);
    } // i_file
    
    // Check that the right number of particles have been read
    if(n_particles_kept!=n_particles_rank)
      SID_log_warning("Rank %d did not receive the right number of particles (ie. %d!=%d)",ERROR_LOGIC,SID.My_rank,n_particles_kept,n_particles_rank);

    // Store everything in the data structure...

    //   ... particle counts ...
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type[i]>0){
        SID_Allreduce(&(n_of_type_rank[i]),&(n_of_type[i]),1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
        ADaPS_store(&(plist->data),(void *)(&(n_of_type_rank[i])),"n_%s",    ADaPS_SCALAR_SIZE_T,pname[i]);
        ADaPS_store(&(plist->data),(void *)(&(n_of_type[i])),     "n_all_%s",ADaPS_SCALAR_SIZE_T,pname[i]);
      }
    }
    ADaPS_store(&(plist->data),(void *)(&n_particles_all),"n_particles_all",ADaPS_SCALAR_SIZE_T);

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

    //  ... ids ...
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type_rank[i]>0)
        ADaPS_store(&(plist->data),(void *)id_array[i],"id_%s",ADaPS_DEFAULT,pname[i]);
    }

    SID_free(SID_FARG id_list_index);

    SID_log("Done.",SID_LOG_CLOSE);
  }
  else
    SID_trap_error("Could not find file with root {%s}",ERROR_IO_OPEN,filename_root_in);
}


int main(int argc, char *argv[]){
  plist_info  plist;
  char        filename_groups_root[256];
  char        filename_snapshot_root[256];
  char        filename_snapshot[256];
  char       *filename_number;
  char        filename_output_properties_dir[256];
  char        filename_output_properties[256];
  char        filename_output_profiles_dir[256];
  char        filename_output_profiles[256];
  char        filename_output_indices[256];
  char        filename_output_indices_dir[256];
  char        filename_output_properties_temp[256];
  char        filename_output_profiles_temp[256];
  char        filename_output_indices_temp[256];
  char        group_text_prefix[4];
  int         n_groups_process;
  int         n_groups;
  int         n_groups_all;
  int         i_rank;
  int         i_group;
  int         i_file_lo;
  int         i_file_hi;
  int         i_file;
  int         i_file_skip;
  int         i_particle;
  int         j_particle;
  int         i_process;
  int         n_particles;
  int         n_particles_max;
  GBPREAL       *x_array;
  GBPREAL       *y_array;
  GBPREAL       *z_array;
  GBPREAL       *vx_array;
  GBPREAL       *vy_array;
  GBPREAL       *vz_array;
  int        *n_particles_groups_process;
  int        *n_particles_groups;
  int        *n_particles_subgroups;
  size_t     *group_offset;
  size_t      n_particles_in_groups;
  size_t     *ids_snapshot;
  size_t     *ids_groups;
  size_t     *ids_sort_index;
  size_t     *ids_snapshot_sort_index;
  size_t     *ids_groups_sort_index;
  size_t      n_particles_snapshot;
  int         i_value;
  double      h_Hubble;
  double      Omega_M;
  double      Omega_b;
  double      Omega_Lambda;
  double      f_gas;
  double      Omega_k;
  double      sigma_8;
  double      n_spec;
  double      redshift;
  double      expansion_factor;
  double      box_size;
  double      particle_mass;
  
  int         r_val;
  struct      stat file_stats;
  size_t      n_bytes;
  size_t      n_bytes_buffer;
  void       *buffer;
  
  FILE       *fp_properties;
  FILE       *fp_profiles;
  FILE       *fp_indices;
  cosmo_info *cosmo;
  halo_properties_info  properties;
  halo_profile_info     profile;
  int                   n_temp;
  int                   n_truncated_local;
  int                   n_truncated;
  int                   largest_truncated;
  int                   largest_truncated_local;
  int flag_write_properties=TRUE;
  int flag_write_profiles  =TRUE;
  int flag_write_indices   =TRUE;
  int flag_manual_centre   =TRUE;
  //flag_manual_centre=FALSE;

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  strcpy(filename_snapshot_root,argv[1]);
  strcpy(filename_groups_root,  argv[2]);
  i_file_lo  =atoi(argv[3]);
  i_file_hi  =atoi(argv[4]);
  i_file_skip=atoi(argv[5]);

  SID_log("Processing group/subgroup statistics for files #%d->#%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file_lo,i_file_hi);
  SID_log("Properties structure size=%d bytes",SID_LOG_COMMENT,sizeof(halo_properties_info));
  SID_log("Profiles   structure size=%d bytes",SID_LOG_COMMENT,sizeof(halo_profile_bin_info));

  for(i_file=i_file_lo;i_file<=i_file_hi;i_file+=i_file_skip){
    filename_number=(char *)SID_malloc(sizeof(char)*10);
    sprintf(filename_number,"%03d", i_file);
    SID_log("Processing file #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file);

    // Read group and particle info
    init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
    //sprintf(filename_number,"read_catalog");
    ADaPS_store(&(plist.data),(void *)filename_number,"read_catalog",ADaPS_DEFAULT);
    read_groups(filename_groups_root,i_file,READ_GROUPS_ALL,&plist,filename_number);
    n_particles_in_groups=((size_t *)ADaPS_fetch(plist.data,"n_particles_%s", filename_number))[0];
    n_groups_all         =((int    *)ADaPS_fetch(plist.data,"n_groups_all_%s",filename_number))[0];
    if(n_groups_all>0){
      read_gadget_binary_local(filename_snapshot_root,i_file,&plist);
      //read_gadget_binary(filename_snapshot_root,i_file,&plist,READ_GADGET_DEFAULT);
      if(ADaPS_exist(plist.data,"id_dark")){
         n_particles_snapshot=((size_t *)ADaPS_fetch(plist.data,"n_dark"))[0];
         ids_snapshot        = (size_t *)ADaPS_fetch(plist.data,"id_dark");
         ids_groups          = (size_t *)ADaPS_fetch(plist.data,"particle_ids_%s",filename_number);
         particle_mass       =((double *)ADaPS_fetch(plist.data,"mass_array_dark"))[0];
         x_array             =(GBPREAL *)ADaPS_fetch(plist.data,"x_dark");
         y_array             =(GBPREAL *)ADaPS_fetch(plist.data,"y_dark");
         z_array             =(GBPREAL *)ADaPS_fetch(plist.data,"z_dark");
         vx_array            =(GBPREAL *)ADaPS_fetch(plist.data,"vx_dark");
         vy_array            =(GBPREAL *)ADaPS_fetch(plist.data,"vy_dark");
         vz_array            =(GBPREAL *)ADaPS_fetch(plist.data,"vz_dark");
      }
      else{
         n_particles_snapshot=0;
      }

      // Initialize cosmology
      box_size        =((double *)ADaPS_fetch(plist.data,"box_size"))[0];
      h_Hubble        =((double *)ADaPS_fetch(plist.data,"h_Hubble"))[0];
      redshift        =((double *)ADaPS_fetch(plist.data,"redshift"))[0];
      expansion_factor=((double *)ADaPS_fetch(plist.data,"expansion_factor"))[0];
      Omega_M         =((double *)ADaPS_fetch(plist.data,"Omega_M"))[0];
      Omega_Lambda    =((double *)ADaPS_fetch(plist.data,"Omega_Lambda"))[0];
      f_gas  =Omega_b/Omega_M;
      Omega_k=1.-Omega_Lambda-Omega_M;
      Omega_b=0.;
      sigma_8=0.;
      n_spec =0.;
      init_cosmo(&cosmo,
                 Omega_Lambda,
                 Omega_M,
                 Omega_k,
                 Omega_b,
                 f_gas,
                 h_Hubble,
                 sigma_8,
                 n_spec);


      // Compute sort indices for ids in snapshot and group catalog
      SID_log("Sorting IDs...",SID_LOG_OPEN|SID_LOG_TIMER);
      if(n_particles_snapshot>0){
         merge_sort((void *)ids_snapshot,n_particles_snapshot, &ids_snapshot_sort_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
         merge_sort((void *)ids_groups,  n_particles_in_groups,&ids_groups_sort_index,  SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);

         // Create snapshot-indices for each particle in the group catalog
         ids_sort_index=(size_t *)SID_malloc(sizeof(size_t)*n_particles_in_groups);
         for(i_particle=0,j_particle=0;i_particle<n_particles_in_groups;i_particle++){
           while(ids_snapshot[ids_snapshot_sort_index[j_particle]]<ids_groups[ids_groups_sort_index[i_particle]]){
             j_particle++;
             if(j_particle>=n_particles_snapshot)
               SID_trap_error("There's a particle id in the group catalog that's not in the snapshot!",ERROR_LOGIC);
           }
           ids_sort_index[ids_groups_sort_index[i_particle]]=ids_snapshot_sort_index[j_particle];
         }
      }
      else{
         ids_snapshot_sort_index=NULL;
         ids_groups_sort_index  =NULL;
         ids_sort_index         =NULL;
      }
      SID_log("Done.",SID_LOG_CLOSE);

      // Print stats; groups first and then subgroups
      for(i_process=0;i_process<2;i_process++){
      
        // Initialize a bunch of stuff depending on whether we are
        //   computing group or subgroup properties
        switch(i_process){
        case 0:
          sprintf(group_text_prefix,"");
          break;
        case 1:
          sprintf(group_text_prefix,"sub");
          break;
        }
        n_groups    =((int *)ADaPS_fetch(plist.data,"n_%sgroups_%s",    group_text_prefix,filename_number))[0];
        n_groups_all=((int *)ADaPS_fetch(plist.data,"n_%sgroups_all_%s",group_text_prefix,filename_number))[0];

        // Fetch some stuff
        n_particles_groups = (int    *)ADaPS_fetch(plist.data,"n_particles_%sgroup_%s"    ,group_text_prefix,filename_number);
        group_offset       = (size_t *)ADaPS_fetch(plist.data,"particle_offset_%sgroup_%s",group_text_prefix,filename_number);

        // Create filenames, directories, etc
        sprintf(filename_output_properties_temp,"%s_%s.catalog_%sgroups_properties",filename_groups_root,filename_number,group_text_prefix);
        sprintf(filename_output_profiles_temp,  "%s_%s.catalog_%sgroups_profiles",  filename_groups_root,filename_number,group_text_prefix);
        sprintf(filename_output_indices_temp,   "%s_%s.catalog_%sgroups_indices",   filename_groups_root,filename_number,group_text_prefix);
        if(SID.n_proc>1){
           // Create property filenames
           if(flag_write_properties){
              char properties_root[MAX_FILENAME_LENGTH];
              strcpy(properties_root,filename_output_properties_temp);
              strip_path(properties_root);
              strcpy(filename_output_properties_dir,filename_output_properties_temp);
              if(SID.I_am_Master)
                mkdir(filename_output_properties_dir,02755);
              SID_Barrier(SID.COMM_WORLD);
              sprintf(filename_output_properties,"%s/%s.%d",filename_output_properties_dir,properties_root,SID.My_rank);
           }
           // Create profile filenames
           if(flag_write_profiles){
              char profiles_root[MAX_FILENAME_LENGTH];
              strcpy(profiles_root,filename_output_profiles_temp);
              strip_path(profiles_root);
              strcpy(filename_output_profiles_dir,filename_output_profiles_temp);
              if(SID.I_am_Master)
                mkdir(filename_output_profiles_dir,02755);
              SID_Barrier(SID.COMM_WORLD);
              sprintf(filename_output_profiles,"%s/%s.%d",filename_output_profiles_dir,profiles_root,SID.My_rank);
           }
           // Create sort indices filenames
           if(flag_write_indices){
              char indices_root[MAX_FILENAME_LENGTH];
              strcpy(indices_root,filename_output_indices_temp);
              strip_path(indices_root);
              strcpy(filename_output_indices_dir,filename_output_indices_temp);
              if(SID.I_am_Master)
                mkdir(filename_output_indices_dir,02755);
              SID_Barrier(SID.COMM_WORLD);
              sprintf(filename_output_indices,"%s/%s.%d",filename_output_indices_dir,indices_root,SID.My_rank);
           }
        }
        else{
           strcpy(filename_output_properties,filename_output_properties_temp);
           strcpy(filename_output_profiles,  filename_output_profiles_temp);
           strcpy(filename_output_indices,   filename_output_indices_temp);
        }
        SID_Barrier(SID.COMM_WORLD); // This makes sure that any created directories are there before proceeding

        // Open files
        SID_log("Processing %sgroups...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix);
        if(flag_write_properties)
           fp_properties=fopen(filename_output_properties,"w");
        else
           fp_properties=NULL;
        if(flag_write_profiles)
           fp_profiles=fopen(filename_output_profiles,"w");
        else
           fp_profiles=NULL;
        if(flag_write_indices)
           fp_indices=fopen(filename_output_indices,"w");
        else
           fp_indices=NULL;

        // Write header
        if(fp_properties!=NULL){
           fwrite(&(SID.My_rank), sizeof(int),1,fp_properties);
           fwrite(&(SID.n_proc),  sizeof(int),1,fp_properties);
           fwrite(&n_groups,      sizeof(int),1,fp_properties);
           fwrite(&n_groups_all,  sizeof(int),1,fp_properties);
        }
        if(fp_profiles!=NULL){
           fwrite(&(SID.My_rank), sizeof(int),1,fp_profiles);
           fwrite(&(SID.n_proc),  sizeof(int),1,fp_profiles);
           fwrite(&n_groups,      sizeof(int),1,fp_profiles);
           fwrite(&n_groups_all,  sizeof(int),1,fp_profiles);
        }
        if(fp_indices!=NULL){
           fwrite(&(SID.My_rank), sizeof(int),1,fp_indices);
           fwrite(&(SID.n_proc),  sizeof(int),1,fp_indices);
           fwrite(&n_groups,      sizeof(int),1,fp_indices);
           fwrite(&n_groups_all,  sizeof(int),1,fp_indices);
        }

        // Turn off the gsl error handler
        gsl_error_handler_t *original_handler;
        original_handler=gsl_set_error_handler_off();

        // Create and write the properties and profiles of each group/subgroup in turn

        // Allocate some temporary arrays for particle positions/velocities
        int     n_particles_alloc;
        double *x;
        double *y;
        double *z;
        double *vx;
        double *vy;
        double *vz;
        double *R;
        size_t *R_index=NULL;
        calc_max(n_particles_groups,&n_particles_alloc,n_groups,SID_INT,CALC_MODE_DEFAULT);
        x                =(double *)SID_malloc(sizeof(double)*n_particles_alloc);
        y                =(double *)SID_malloc(sizeof(double)*n_particles_alloc);
        z                =(double *)SID_malloc(sizeof(double)*n_particles_alloc);
        vx               =(double *)SID_malloc(sizeof(double)*n_particles_alloc);
        vy               =(double *)SID_malloc(sizeof(double)*n_particles_alloc);
        vz               =(double *)SID_malloc(sizeof(double)*n_particles_alloc);
        R                =(double *)SID_malloc(sizeof(double)*n_particles_alloc);

        for(i_group=0,n_truncated_local=0,largest_truncated_local=0;i_group<n_groups;i_group++){
          if(compute_group_analysis(&properties,
                                    &profile,
                                    ids_snapshot,
                                    x_array,
                                    y_array,
                                    z_array,
                                    vx_array,
                                    vy_array,
                                    vz_array,
                                    &(ids_sort_index[group_offset[i_group]]),
                                    box_size,
                                    h_Hubble,
                                    Omega_M,
                                    particle_mass,
                                    n_particles_groups[i_group],
                                    redshift,
                                    expansion_factor,
                                    i_group,
                                    x,y,z,vx,vy,vz,R,&R_index,
                                    flag_manual_centre,
                                    cosmo)!=TRUE){
             n_truncated_local++;
             largest_truncated_local=MAX(largest_truncated_local,n_particles_groups[i_group]);
          }
          write_group_analysis(fp_properties,
                               fp_profiles,
                               fp_indices,
                               &properties,
                               &profile,
                               R_index,
                               n_particles_groups[i_group]);

          SID_free(SID_FARG R_index);
        }

        // Free arrays
        SID_free(SID_FARG x);
        SID_free(SID_FARG y);
        SID_free(SID_FARG z);
        SID_free(SID_FARG vx);
        SID_free(SID_FARG vy);
        SID_free(SID_FARG vz);
        SID_free(SID_FARG R);

        if(fp_properties!=NULL)
          fclose(fp_properties);
        if(fp_profiles!=NULL)
          fclose(fp_profiles);
        if(fp_indices!=NULL)
          fclose(fp_indices);
        calc_sum_global(&n_truncated_local,      &n_truncated,      1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
        calc_max_global(&largest_truncated_local,&largest_truncated,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD); 

        // Restore the original handler
        gsl_set_error_handler(original_handler);

        SID_Barrier(SID.COMM_WORLD);
        SID_log("Done. (f_truncated=%6.2lf%% largest=%d)",SID_LOG_CLOSE,100.*(double)n_truncated/(double)n_groups_all,largest_truncated);
      }

      // Clean-up
      free_cosmo(&cosmo);
      SID_free(SID_FARG ids_snapshot_sort_index);
      SID_free(SID_FARG ids_groups_sort_index);
      SID_free(SID_FARG ids_sort_index);
    }
    // If the group catalog or snapshot is empty, create an empty file
    else{
      SID_log("Creating empty analysis files...",SID_LOG_OPEN);
      if(SID.I_am_Master){
         for(i_process=0;i_process<2;i_process++){
            switch(i_process){
            case 0:
              sprintf(group_text_prefix,"");
              break;
            case 1:
              sprintf(group_text_prefix,"sub");
              break;
            }
            n_groups_all=0;
            n_temp      =1;
            if(flag_write_properties){
               sprintf(filename_output_properties,"%s_%s.catalog_%sgroups_properties",filename_groups_root,filename_number,group_text_prefix);
               fp_properties=fopen(filename_output_properties,"w");
               fwrite(&(SID.My_rank), sizeof(int),1,fp_properties);
               fwrite(&n_temp,        sizeof(int),1,fp_properties);
               fwrite(&n_groups_all,  sizeof(int),1,fp_properties);
               fwrite(&n_groups_all,  sizeof(int),1,fp_properties);
               fclose(fp_properties);
            }
            if(flag_write_profiles){
               sprintf(filename_output_profiles,  "%s_%s.catalog_%sgroups_profiles",filename_groups_root,filename_number,group_text_prefix);
               fp_profiles  =fopen(filename_output_profiles,  "w");
               fwrite(&(SID.My_rank), sizeof(int),1,fp_profiles);
               fwrite(&n_temp,        sizeof(int),1,fp_profiles);
               fwrite(&n_groups_all,  sizeof(int),1,fp_profiles);
               fwrite(&n_groups_all,  sizeof(int),1,fp_profiles);
               fclose(fp_profiles);
            }
         }
      }
      SID_log("Done.",SID_LOG_CLOSE);
    }
    free_plist(&plist);
    SID_log("Done.",SID_LOG_CLOSE);
  }

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

