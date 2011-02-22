#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpSPH.h>

void write_gadget_binary(char        *filename_in,
                         plist_info  *plist){
  char    **pname;
  int       i,j,k,jj; 
  char      filename[256];
  size_t    n_of_type[N_GADGET_TYPE];
  size_t    n_of_type_rank[N_GADGET_TYPE];
  size_t  **n_of_type_file;
  size_t  **n_of_type_file_rank;
  size_t   *i_p_hi_file;
  size_t   *n_particles_file;
  size_t   *n_gas_file;
  size_t   *n_particles_multimass_file;
  size_t    n_particles;
  size_t    n_particles_rank;
  uint      n_temp_uint[N_GADGET_TYPE];
  int       n_temp_int[N_GADGET_TYPE];
  int       n_files;
  double    mass_array[N_GADGET_TYPE];
  int       unused[256];
  double    h_Hubble;
  double    d_value;
  double   *d1_array;
  REAL     *R1_array;
  REAL     *R2_array;
  REAL     *R3_array;
  float     f_temp;
  float     f1_temp;
  float     f2_temp;
  float     f3_temp;
  int       i_value;
  size_t   *id_array;
  int       record_length;
  int       n_return;
  size_t    s_load;
  int       i_file;
  int       j_file;
  int       i_rank;
  size_t    i_p;
  long long i_p_LL;
  int       i_p_I;
  int       flag_multifile;
  int       flag_multimass;
  int       flag_multimass_species[N_GADGET_TYPE];
  int       flag_noTemps;
  int       flag_recordwritten;
  int       flag_LONGIDS;
  int       n_x_replicate;
  int       n_y_replicate;
  int       n_z_replicate;
  int       i_x_replicate;
  int       i_y_replicate;
  int       i_z_replicate;
  size_t    i_replicate;
  size_t    n_replicate;
  double    box_size;
  double    box_size_offset;
  REAL      x_offset;
  REAL      y_offset;
  REAL      z_offset;
  FILE     *fp;
  void     *buffer;
  size_t    n_buffer;
  size_t    n_write;
  int       n_non_zero;
  int       flag_recordlength_written;

  SID_log("Writing GADGET bindary file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in);      

  pname=plist->species;

  // Determine if we are supposed to replicate this box in x &/or y &/or z
  if(ADaPS_exist(plist->data,"n_x_replicate"))
    n_x_replicate=((int *)ADaPS_fetch(plist->data,"n_x_replicate"))[0];
  else
    n_x_replicate=1;
  if(ADaPS_exist(plist->data,"n_y_replicate"))
    n_y_replicate=((int *)ADaPS_fetch(plist->data,"n_y_replicate"))[0];
  else
    n_y_replicate=1;
  if(ADaPS_exist(plist->data,"n_z_replicate"))
    n_z_replicate=((int *)ADaPS_fetch(plist->data,"n_z_replicate"))[0];
  else
    n_z_replicate=1;
  n_replicate=(size_t)(n_x_replicate*n_y_replicate*n_z_replicate);

  // We're going to need h_Hubble for a bunch of stuff
  if(ADaPS_exist(plist->data,"h_Hubble"))
    h_Hubble=((double *)ADaPS_fetch(plist->data,"h_Hubble"))[0];
  else
    h_Hubble=1.;

  // Fetch the number of files
  if(ADaPS_exist(plist->data,"n_files"))
    n_files=((int *)ADaPS_fetch(plist->data,"n_files"))[0];
  else
    n_files=1;
  if(n_files>1)
    flag_multifile=TRUE;
  else
    flag_multifile=FALSE;

  // Set various particle counts
  n_particles     =0;
  n_particles_rank=0;
  for(i=0,n_non_zero=0,j=0;i<N_GADGET_TYPE;i++){
    if(i<plist->n_species){
      // Set numbers of particles (total across all ranks)
      if(ADaPS_exist(plist->data,"n_all_%s",pname[i]))
        n_of_type[i]=((size_t *)ADaPS_fetch(plist->data,"n_all_%s",pname[i]))[0];
      else
        n_of_type[i]=0;
      n_particles+=n_of_type[i];
      // Set numbers of particles (local)
      if(ADaPS_exist(plist->data,"n_%s",pname[i]))
        n_of_type_rank[i]=((size_t *)ADaPS_fetch(plist->data,"n_%s",pname[i]))[0];
      else
        n_of_type_rank[i]=0;
      n_particles_rank+=n_of_type_rank[i];
    }
    else{
      n_of_type[i]     =0;
      n_of_type_rank[i]=0;
    }
    if(n_of_type[i]>0)
      n_non_zero++;
  }

  // Determine multimass info
  flag_multimass=FALSE;
  for(i=0;i<N_GADGET_TYPE;i++){
    flag_multimass_species[i]=FALSE;
    if(n_of_type[i]>0){
      if(ADaPS_exist(plist->data,"mass_array_%s",pname[i]))
        mass_array[i]=
          (((double *)ADaPS_fetch(plist->data,"mass_array_%s",pname[i]))[0])/
          (plist->mass_unit/h_Hubble);
      else
        mass_array[i] =0.;
      if(mass_array[i]==0.){
        flag_multimass_species[i]=TRUE;
        flag_multimass           =TRUE;
      }
    }
    else
      mass_array[i] =0.;
  }

  // Decide how each rank's particles
  //  will be divided between files
  i_p_hi_file=(size_t *)SID_malloc(sizeof(size_t)*n_files);
  for(i_file=0;i_file<n_files;i_file++){
    if(i_file==n_files-1)
      i_p_hi_file[i_file]=n_particles-1;
    else
      i_p_hi_file[i_file]=n_particles*(i_file+1)/n_files-1;
  }

  // Figure-out how the particles will be split across files
  n_of_type_file            =(size_t **)SID_malloc(sizeof(size_t *)*n_files);
  n_of_type_file_rank       =(size_t **)SID_malloc(sizeof(size_t *)*n_files);
  n_particles_file          =(size_t  *)SID_malloc(sizeof(size_t  )*n_files);
  n_gas_file                =(size_t  *)SID_malloc(sizeof(size_t  )*n_files);
  n_particles_multimass_file=(size_t  *)SID_malloc(sizeof(size_t  )*n_files);
  for(i_file=0;i_file<n_files;i_file++){
    n_of_type_file[i_file]     =(size_t *)SID_malloc(sizeof(size_t)*N_GADGET_TYPE);
    n_of_type_file_rank[i_file]=(size_t *)SID_malloc(sizeof(size_t)*N_GADGET_TYPE);
    for(i=0;i<N_GADGET_TYPE;i++){
      n_of_type_file[i_file][i]       =0;
      n_of_type_file_rank[i_file][i]  =0;
    }
    n_gas_file[i_file]                =0;
    n_particles_file[i_file]          =0;
    n_particles_multimass_file[i_file]=0;
  }
  for(i=0,i_p=0,i_file=0;i<N_GADGET_TYPE;i++) {
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
      if(i_rank==SID.My_rank){
        if(n_of_type_rank[i]>0){
          for(j=0;j<n_of_type_rank[i];j++){
            n_of_type_file[i_file][i]++;
            n_of_type_file_rank[i_file][i]++;
            n_particles_file[i_file]++;
            if(flag_multimass_species[i])
              n_particles_multimass_file[i_file]++;
            if(i==GADGET_TYPE_GAS)
              n_gas_file[i_file]++;
            i_p++;
            if(i_p>i_p_hi_file[i_file])
              i_file++;
          }
        }
      }
#ifdef USE_MPI
      MPI_Bcast(&i_p,   1,MPI_SIZE_T, i_rank,MPI_COMM_WORLD);
      MPI_Bcast(&i_file,1,MPI_INTEGER,i_rank,MPI_COMM_WORLD);
      for(j=0;j<n_files;j++){
        MPI_Bcast(  n_of_type_file[j],             N_GADGET_TYPE,MPI_SIZE_T,i_rank,MPI_COMM_WORLD);
        MPI_Bcast(&(n_particles_file[j]),          1,            MPI_SIZE_T,i_rank,MPI_COMM_WORLD);
        MPI_Bcast(&(n_particles_multimass_file[j]),1,            MPI_SIZE_T,i_rank,MPI_COMM_WORLD);
        MPI_Bcast(&(n_gas_file[j]),                1,            MPI_SIZE_T,i_rank,MPI_COMM_WORLD);
      }
#endif
    }
  }
 /* 
  for(i_file=0;i_file<n_files;i_file++){
    for(i=0;i<N_GADGET_TYPE;i++) {
      for(i_rank=0;i_rank<SID.n_proc;i_rank++){
        if(i_rank==SID.My_rank){
          fprintf(stderr,"rank=%d of %d type=%d of %d n_of_type_file=%d n_particles_file=%d\n",SID.My_rank,n_files,i,N_GADGET_TYPE,n_of_type_file[i_file][i],n_particles_file[i_file]);
        }
        SID_Barrier(SID.COMM_WORLD);
      }
    }
  }
*/

  // Write the headers of all the files
  if(n_particles>0){
    SID_log("Writing headers...",SID_LOG_OPEN);
    SID_log("%lld",SID_LOG_CONTINUE,n_particles);
    if(n_non_zero>0)
      SID_log(" (",SID_LOG_CONTINUE,n_particles);
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type[i]>0){
        if(i==n_non_zero-1){
          if(n_non_zero>1)
            SID_log("and %lld %s",SID_LOG_CONTINUE,n_of_type[i],pname[i]);
          else
            SID_log("%lld %s",SID_LOG_CONTINUE,n_of_type[i],pname[i]);
        }
        else{
          if(n_non_zero>1)
            SID_log("%lld %s, ",SID_LOG_CONTINUE,n_of_type[i],pname[i]);
          else
            SID_log("%lld %s",SID_LOG_CONTINUE,n_of_type[i],pname[i]);
        }
      }
    }
    if(n_non_zero>0)
      SID_log(") particles...",SID_LOG_CONTINUE);
    else
      SID_log(" particles...",SID_LOG_CONTINUE);
    if(n_replicate>1)
      SID_log("%dx replication...",SID_LOG_CONTINUE,n_replicate);
      
    for(i_file=0;i_file<n_files;i_file++){
      if(SID.I_am_Master){
        if(n_files>1)
          sprintf(filename,"%s.%d",filename_in,i_file);
        else
          sprintf(filename,"%s",filename_in);
        if((fp=fopen(filename,"w"))!=NULL){
          // Header record length
          record_length=GADGET_HEADER_SIZE;
          fwrite(&record_length,4,1,fp);
          // No. of particles
          for(i=0;i<N_GADGET_TYPE;i++)
            n_temp_int[i]=(int)(((size_t)(n_replicate*n_of_type_file[i_file][i]) <<32)>>32);
          n_return=fwrite(n_temp_int,sizeof(int),N_GADGET_TYPE,fp); 
          s_load  =n_return*sizeof(int);
          // Particle mass array
          n_return =fwrite(mass_array,sizeof(double),N_GADGET_TYPE,fp); 
          s_load  +=n_return*sizeof(double);
          // Expansion factor (or time)
          if(ADaPS_exist(plist->data,"expansion_factor"))
            d_value=((double *)ADaPS_fetch(plist->data,"expansion_factor"))[0];
          else if(ADaPS_exist(plist->data,"time"))
            d_value=((double *)ADaPS_fetch(plist->data,"time"))[0];
          else
            d_value=1.;
          n_return=fwrite(&d_value,sizeof(double),1,fp); 
          s_load+=n_return*sizeof(double);
          // Redshift 
          if(ADaPS_exist(plist->data,"redshift"))
            d_value=((double *)ADaPS_fetch(plist->data,"redshift"))[0];
          else
            d_value=0.;
          n_return=fwrite(&d_value,sizeof(double),1,fp); 
          s_load+=n_return*sizeof(double);
          // Some flags 
          if(ADaPS_exist(plist->data,"flag_Sfr"))
            i_value=((int *)ADaPS_fetch(plist->data,"flag_Sfr"))[0];
          else
            i_value=0;
          n_return =fwrite(&i_value,sizeof(int),1,fp); 
          s_load  +=n_return*sizeof(int);
          if(ADaPS_exist(plist->data,"flag_feedback"))
            i_value=((int *)ADaPS_fetch(plist->data,"flag_feedback"))[0];
          else
            i_value=0;
          n_return=fwrite(&i_value,sizeof(int),1,fp); 
          s_load+=n_return*sizeof(int);
          // Number of particles in all files 
          for(i=0;i<N_GADGET_TYPE;i++)
            n_temp_uint[i]=(unsigned int)(((size_t)(n_replicate*n_of_type[i]) <<32)>>32);
          n_return=fwrite(n_temp_uint,sizeof(unsigned int),N_GADGET_TYPE,fp); 
          s_load  +=n_return*sizeof(unsigned int);
          // Another flag 
          if(ADaPS_exist(plist->data,"flag_cooling"))
            i_value=((int *)ADaPS_fetch(plist->data,"flag_cooling"))[0];
          else
            i_value=0;
          n_return=fwrite(&i_value,sizeof(int),1,fp); 
          s_load  +=n_return*sizeof(int);
          // Number of files per snapshot 
          n_return=fwrite(&n_files,sizeof(int),1,fp); 
          s_load  +=n_return*sizeof(int);
          // Box size 
          if(ADaPS_exist(plist->data,"box_size"))
            box_size=((double *)ADaPS_fetch(plist->data,"box_size"))[0];
          else
            box_size=0.;
          box_size_offset=box_size;
          if((REAL)box_size<=0.)
            SID_trap_error("Box size is <=0 in write_gadget_binary.\n",ERROR_LOGIC);
          if(n_x_replicate!=n_y_replicate || n_x_replicate!=n_z_replicate)
            SID_log("WARNING: replication is not cubical ... using n_x_replicate to scale box size!",SID_LOG_OPEN);
          if(n_x_replicate!=1)
            box_size*=(double)n_x_replicate;
          d_value =box_size/(plist->length_unit/h_Hubble);
          n_return=fwrite(&d_value,sizeof(double),1,fp); 
          s_load  +=n_return*sizeof(double);
        
          // Cosmology 
          // Omega_o 
          if(ADaPS_exist(plist->data,"Omega_M"))
            d_value=((double *)ADaPS_fetch(plist->data,"Omega_M"))[0];
          else
            d_value=0.;
          n_return=fwrite(&d_value,sizeof(double),1,fp); 
          s_load+=n_return*sizeof(double);
          // Omega_Lambda 
          if(ADaPS_exist(plist->data,"Omega_Lambda"))
            d_value=((double *)ADaPS_fetch(plist->data,"Omega_Lambda"))[0];
          else
            d_value=0.;
          n_return=fwrite(&d_value,sizeof(double),1,fp); 
          s_load+=n_return*sizeof(double);
          // Hubble parameter 
          n_return =fwrite(&h_Hubble,sizeof(double),1,fp); 
          s_load  +=n_return*sizeof(double);

          // Misc. stuff
          // Another flag 
          if(ADaPS_exist(plist->data,"flag_metals"))
            i_value=((int *)ADaPS_fetch(plist->data,"flag_metals"))[0];
          else
            i_value=0;
          n_return =fwrite(&i_value,sizeof(int),1,fp); 
          s_load  +=n_return*sizeof(int);
          // Another flag 
          if(ADaPS_exist(plist->data,"flag_ages"))
            i_value=((int *)ADaPS_fetch(plist->data,"flag_ages"))[0];
          else
            i_value=0;
          n_return =fwrite(&i_value,sizeof(int),1,fp); 
          s_load  +=n_return*sizeof(int);
          
          // High words for VERY LARGE total particle numbers
          flag_LONGIDS=FALSE;
          for(i=0;i<N_GADGET_TYPE;i++){
            n_temp_uint[i]=(unsigned int)((size_t)(n_replicate*n_of_type[i]) >> 32);
            if(n_temp_uint[i]>0)
              flag_LONGIDS=TRUE;
          }
          n_return =fwrite(n_temp_uint,sizeof(unsigned int),N_GADGET_TYPE,fp); 
          s_load  +=n_return*sizeof(unsigned int);

          // Another flag 
          if(ADaPS_exist(plist->data,"flag_entropyICs"))
            i_value=((int *)ADaPS_fetch(plist->data,"flag_cooling"))[0];
          else
            i_value=0;
          n_return =fwrite(&i_value,sizeof(int),1,fp); 
          s_load  +=n_return*sizeof(int);
          // Unused space 
          n_return=fwrite(unused,sizeof(int),(GADGET_HEADER_SIZE-s_load)/sizeof(int),fp);
          n_return=fwrite(&record_length,4,1,fp);
          fclose(fp);
        }
        else
          SID_log("Error writing headers!",SID_LOG_COMMENT);
      }
    }
#ifdef USE_MPI
    MPI_Bcast(&flag_LONGIDS,   1,MPI_INTEGER,MASTER_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&box_size,       1,MPI_DOUBLE, MASTER_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&box_size_offset,1,MPI_DOUBLE, MASTER_RANK,MPI_COMM_WORLD);
#endif
    if(flag_LONGIDS)
      SID_log("(using LONG IDS)...",SID_LOG_CONTINUE);      
    SID_log("Done.",SID_LOG_CLOSE);      

    // Allocate a buffer to make writes faster
    buffer=(void *)SID_malloc(MAX(MAX(3*sizeof(float),sizeof(long long)),sizeof(double))*GADGET_BUFFER_SIZE);

    // Write positions
    SID_log("Writing positions...",SID_LOG_OPEN);
    for(i_file=0,i_p=0;i_file<n_files;i_file++){
      if(n_files>1)
        sprintf(filename,"%s.%d",filename_in,i_file);
      else
        sprintf(filename,"%s",filename_in);
      record_length=3*n_replicate*n_particles_file[i_file]*sizeof(float);
      for(i=0,flag_recordlength_written=FALSE;i<N_GADGET_TYPE;i++){
        for(i_rank=0;i_rank<SID.n_proc;i_rank++){
          if(SID.My_rank==i_rank){
            fp=fopen(filename,"a");
            if(SID.I_am_Master && !flag_recordlength_written){
              fwrite(&record_length,4,1,fp);
              flag_recordlength_written=TRUE;
            }
            R1_array=(REAL *)ADaPS_fetch(plist->data,"x_%s",pname[i]);
            R2_array=(REAL *)ADaPS_fetch(plist->data,"y_%s",pname[i]);
            R3_array=(REAL *)ADaPS_fetch(plist->data,"z_%s",pname[i]);
            for(k=0,j_file=0;j_file<i_file;j_file++) k+=n_of_type_file_rank[j_file][i];
            for(j=0;j<n_of_type_file_rank[i_file][i];j+=n_buffer,k+=n_buffer){
              n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file_rank[i_file][i]-j);
              for(i_x_replicate=0;i_x_replicate<n_x_replicate;i_x_replicate++){
                x_offset=((REAL)i_x_replicate)*(REAL)box_size_offset;
                for(jj=0;jj<n_buffer;jj++)
                  ((float *)buffer)[3*jj+0]=((float)(R1_array[k+jj]+x_offset))/((REAL)(plist->length_unit/h_Hubble));
                for(i_y_replicate=0;i_y_replicate<n_y_replicate;i_y_replicate++){
                  y_offset=((REAL)i_y_replicate)*(REAL)box_size_offset;
                  for(jj=0;jj<n_buffer;jj++)
                    ((float *)buffer)[3*jj+1]=((float)(R2_array[k+jj]+y_offset))/((REAL)(plist->length_unit/h_Hubble));
                  for(i_z_replicate=0;i_z_replicate<n_z_replicate;i_z_replicate++){
                    z_offset=((REAL)i_z_replicate)*(REAL)box_size_offset;
                    for(jj=0;jj<n_buffer;jj++)
                      ((float *)buffer)[3*jj+2]=((float)(R3_array[k+jj]+z_offset))/((REAL)(plist->length_unit/h_Hubble));
                    fwrite(buffer,sizeof(float),3*n_buffer,fp);
                    i_p+=n_buffer;
                  }
                }
              }
            }
            fclose(fp);
          }
          SID_Barrier(SID.COMM_WORLD);
        }
      }
      if(SID.I_am_Master && flag_recordlength_written){
        fp=fopen(filename,"a");
        fwrite(&record_length,4,1,fp);
        fclose(fp);
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);      
  
    // Write velocities
    SID_log("Writing velocities...",SID_LOG_OPEN);
    for(i_file=0,i_p=0;i_file<n_files;i_file++){
      if(n_files>1)
        sprintf(filename,"%s.%d",filename_in,i_file);
      else
        sprintf(filename,"%s",filename_in);
      record_length=n_replicate*3*n_particles_file[i_file]*sizeof(float);
      for(i=0,flag_recordlength_written=FALSE;i<N_GADGET_TYPE;i++){
        for(i_rank=0;i_rank<SID.n_proc;i_rank++){
          if(SID.My_rank==i_rank){
            fp=fopen(filename,"a");
            if(SID.I_am_Master && !flag_recordlength_written){
              fwrite(&record_length,4,1,fp);
              flag_recordlength_written=TRUE;
            }
            R1_array=(REAL *)ADaPS_fetch(plist->data,"vx_%s",pname[i]);
            R2_array=(REAL *)ADaPS_fetch(plist->data,"vy_%s",pname[i]);
            R3_array=(REAL *)ADaPS_fetch(plist->data,"vz_%s",pname[i]);
            for(k=0,j_file=0;j_file<i_file;j_file++) k+=n_of_type_file_rank[j_file][i];
            for(j=0;j<n_of_type_file_rank[i_file][i];j+=n_buffer,k+=n_buffer){
              n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file_rank[i_file][i]-j);
              for(jj=0;jj<n_buffer;jj++){
                ((float *)buffer)[3*jj+0]=((float)(R1_array[k+jj]))/((REAL)(plist->velocity_unit));
                ((float *)buffer)[3*jj+1]=((float)(R2_array[k+jj]))/((REAL)(plist->velocity_unit));
                ((float *)buffer)[3*jj+2]=((float)(R3_array[k+jj]))/((REAL)(plist->velocity_unit));
              }
              for(i_x_replicate=0;i_x_replicate<n_x_replicate;i_x_replicate++){
                for(i_y_replicate=0;i_y_replicate<n_y_replicate;i_y_replicate++){
                  for(i_z_replicate=0;i_z_replicate<n_z_replicate;i_z_replicate++){
                    fwrite(buffer,sizeof(float),3*n_buffer,fp);
                  }
                }
              }
            }
            fclose(fp);
          }
          SID_Barrier(SID.COMM_WORLD);
        }
      }
      if(SID.I_am_Master && flag_recordlength_written){
        fp=fopen(filename,"a");
        fwrite(&record_length,4,1,fp);
        fclose(fp);
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);      

    // Write ids
    SID_log("Writing ids...",SID_LOG_OPEN);
    for(i_file=0,i_p=0;i_file<n_files;i_file++){
      if(n_files>1)
        sprintf(filename,"%s.%d",filename_in,i_file);
      else
        sprintf(filename,"%s",filename_in);
      if(flag_LONGIDS)
        record_length=n_replicate*n_particles_file[i_file]*sizeof(long long);
      else
        record_length=n_replicate*n_particles_file[i_file]*sizeof(int);
      for(i=0,flag_recordlength_written=FALSE;i<N_GADGET_TYPE;i++){
        for(i_rank=0;i_rank<SID.n_proc;i_rank++){
          if(SID.My_rank==i_rank){
            fp=fopen(filename,"a");
            if(SID.I_am_Master && !flag_recordlength_written){
              fwrite(&record_length,4,1,fp);
              flag_recordlength_written=TRUE;
            }
            // If ids are defined, then write them ...
            if(ADaPS_exist(plist->data,"id_%s",pname[i])){
              id_array=(size_t *)ADaPS_fetch(plist->data,"id_%s",pname[i]);
              for(k=0,j_file=0;j_file<i_file;j_file++) k+=n_of_type_file_rank[j_file][i];
              for(j=0;j<n_of_type_file_rank[i_file][i];j+=n_buffer,k+=n_buffer){
                n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file_rank[i_file][i]-j);
                for(i_x_replicate=0,i_replicate=0;i_x_replicate<n_x_replicate;i_x_replicate++){
                  for(i_y_replicate=0;i_y_replicate<n_y_replicate;i_y_replicate++){
                    for(i_z_replicate=0;i_z_replicate<n_z_replicate;i_z_replicate++,i_replicate++){
                      if(flag_LONGIDS){
                        for(jj=0;jj<n_buffer;jj++)
                          ((long long *)buffer)[jj]=(long long)(id_array[k+jj]+i_replicate*n_particles);
                        fwrite(buffer,sizeof(long long),n_buffer,fp);
                      }
                      else{
                        for(jj=0;jj<n_buffer;jj++)
                          ((int *)buffer)[jj]=(int)(id_array[k+jj]+i_replicate*n_particles);
                        fwrite(buffer,sizeof(int),n_buffer,fp);
                      }
                    }
                  }
                }
              }
            }
            // ... else generate them now
            else{
              for(j=0;j<n_of_type_file_rank[i_file][i];j+=n_buffer,k+=n_buffer){
                n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file_rank[i_file][i]-j);
                for(i_x_replicate=0,i_replicate=0;i_x_replicate<n_x_replicate;i_x_replicate++){
                  for(i_y_replicate=0;i_y_replicate<n_y_replicate;i_y_replicate++){
                    for(i_z_replicate=0;i_z_replicate<n_z_replicate;i_z_replicate++,i_replicate++){
                      if(flag_LONGIDS){
                        for(jj=0;jj<n_buffer;jj++)
                          ((long long *)buffer)[jj]=(long long)(i_p++);
                        fwrite(buffer,sizeof(long long),n_buffer,fp);
                      }
                      else{
                        for(jj=0;jj<n_buffer;jj++)
                          ((int *)buffer)[jj]=(int)(i_p++);
                        fwrite(buffer,sizeof(int),n_buffer,fp);
                      }
                    }
                  }
                }
              }
            }
            fclose(fp);
          }
          SID_Bcast(&i_p,sizeof(size_t),i_rank,SID.COMM_WORLD);
        }
      }
      if(SID.I_am_Master && flag_recordlength_written){
        fp=fopen(filename,"a");
        fwrite(&record_length,4,1,fp);
        fclose(fp);
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);      
    
    // Write masses (if needed)
    if(flag_multimass){
      SID_log("Writing masses...",SID_LOG_OPEN);      
      for(i_file=0,i_p=0;i_file<n_files;i_file++){
        if(n_files>1)
          sprintf(filename,"%s.%d",filename_in,i_file);
        else
          sprintf(filename,"%s",filename_in);
        record_length=n_replicate*n_particles_file[i_file]*sizeof(float);
        for(i=0,flag_recordlength_written=FALSE;i<N_GADGET_TYPE;i++){
          if(flag_multimass_species[i]){
            for(i_rank=0;i_rank<SID.n_proc;i_rank++){
              if(SID.My_rank==i_rank){
                fp=fopen(filename,"a");
                if(SID.I_am_Master && !flag_recordlength_written){
                  fwrite(&record_length,4,1,fp);
                  flag_recordlength_written=TRUE;
                }
                d1_array=(double *)ADaPS_fetch(plist->data,"M_%s",pname[i]);
                for(k=0,j_file=0;j_file<i_file;j_file++) k+=n_of_type_file_rank[j_file][i];
                for(j=0;j<n_of_type_file_rank[i_file][i];j+=n_buffer,k+=n_buffer){
                  n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file_rank[i_file][i]-j);
                  for(jj=0;jj<n_buffer;jj++)
                    ((float *)buffer)[jj]=(float)(d1_array[k+jj]/(plist->mass_unit/h_Hubble));
                  for(i_x_replicate=0;i_x_replicate<n_x_replicate;i_x_replicate++){
                    for(i_y_replicate=0;i_y_replicate<n_y_replicate;i_y_replicate++){
                      for(i_z_replicate=0;i_z_replicate<n_z_replicate;i_z_replicate++){
                        fwrite(buffer,sizeof(float),n_buffer,fp);
                        i_p+=n_buffer;
                      }
                    }
                  }
                }
                fclose(fp);
              }
              SID_Barrier(SID.COMM_WORLD);
            }
          }
        }
        if(SID.I_am_Master && flag_recordlength_written){
          fp=fopen(filename,"a");
          fwrite(&record_length,4,1,fp);
          fclose(fp);
        }
      }
      SID_log("Done.",SID_LOG_CLOSE);      
    }
      
    // Write gas properties ...
    if(n_of_type[GADGET_TYPE_GAS]>0){
      SID_log("Writing gas properties...",SID_LOG_OPEN);
      i=GADGET_TYPE_GAS;

      // Internal energies
      for(i_file=0,i_p=0;i_file<n_files;i_file++){
        if(n_files>1)
          sprintf(filename,"%s.%d",filename_in,i_file);
        else
          sprintf(filename,"%s",filename_in);
        record_length=n_replicate*n_particles_file[i_file]*sizeof(float);
        for(i_rank=0,flag_recordlength_written=FALSE;i_rank<SID.n_proc;i_rank++){
          if(SID.My_rank==i_rank){
            fp=fopen(filename,"a");
            if(SID.I_am_Master && !flag_recordlength_written){
              fwrite(&record_length,4,1,fp);
              flag_recordlength_written=TRUE;
            }
            // If internal energies are defined, then write them ...
            if(ADaPS_exist(plist->data,"u_%s",pname[i])){
              R1_array=(REAL *)ADaPS_fetch(plist->data,"u_%s",pname[i]);
              for(k=0,j_file=0;j_file<i_file;j_file++) k+=n_of_type_file_rank[j_file][i];
              for(j=0;j<n_of_type_file_rank[i_file][i];j+=n_buffer,k+=n_buffer){
                n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file_rank[i_file][i]-j);
                for(jj=0;jj<n_buffer;jj++)
                  ((float *)buffer)[jj]=(float)(R1_array[k+jj])/((REAL)(plist->velocity_unit*plist->velocity_unit));
                for(i_x_replicate=0;i_x_replicate<n_x_replicate;i_x_replicate++){
                  for(i_y_replicate=0;i_y_replicate<n_y_replicate;i_y_replicate++){
                    for(i_z_replicate=0;i_z_replicate<n_z_replicate;i_z_replicate++){
                      fwrite(buffer,sizeof(float),n_buffer,fp);
                      i_p+=n_buffer;
                    }
                  }
                }
              }
            }
            // ... else use temperatures (else write zeros)
            else{
              if(ADaPS_exist(plist->data,"T_%s",pname[i])){
                flag_noTemps=FALSE;
                R1_array=(REAL *)ADaPS_fetch(plist->data,"T_%s",pname[i]);
              }
              else
                flag_noTemps=TRUE;
              for(k=0,j_file=0;j_file<i_file;j_file++) k+=n_of_type_file_rank[j_file][i];
              for(j=0;j<n_of_type_file_rank[i_file][i];j+=n_buffer,k+=n_buffer){
                n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file_rank[i_file][i]-j);
                if(flag_noTemps){
                  for(jj=0;jj<n_buffer;jj++)
                    ((float *)buffer)[jj]=0.;
                }
                else{
                  for(jj=0;jj<n_buffer;jj++)
                    ((float *)buffer)[jj]=(float)(1.5*K_BOLTZMANN*R1_array[k+jj])/
                      (float)(MU_MMW*M_PROTON*plist->velocity_unit*plist->velocity_unit);
                }
                for(i_x_replicate=0;i_x_replicate<n_x_replicate;i_x_replicate++){
                  for(i_y_replicate=0;i_y_replicate<n_y_replicate;i_y_replicate++){
                    for(i_z_replicate=0;i_z_replicate<n_z_replicate;i_z_replicate++){
                      fwrite(buffer,sizeof(float),n_buffer,fp);
                      i_p+=n_buffer;
                    }
                  }
                }
              }
            }
            fclose(fp);
          }
          SID_Barrier(SID.COMM_WORLD);
        }
        if(SID.I_am_Master && flag_recordlength_written){
          fp=fopen(filename,"a");
          fwrite(&record_length,4,1,fp);
          fclose(fp);
        }
      }

      // Densities
      for(i_file=0,i_p=0;i_file<n_files;i_file++){
        if(n_files>1)
          sprintf(filename,"%s.%d",filename_in,i_file);
        else
          sprintf(filename,"%s",filename_in);
        record_length=n_replicate*n_particles_file[i_file]*sizeof(float);
        for(i_rank=0,flag_recordlength_written=FALSE;i_rank<SID.n_proc;i_rank++){
          if(SID.My_rank==i_rank){
            fp=fopen(filename,"a");
            if(SID.I_am_Master && !flag_recordlength_written){
              fwrite(&record_length,4,1,fp);
              flag_recordlength_written=TRUE;
            }
            // If densities are defined, then write them ...
            if(ADaPS_exist(plist->data,"rho_%s",pname[i])){
              R1_array=(REAL *)ADaPS_fetch(plist->data,"rho_%s",pname[i]);
              for(k=0,j_file=0;j_file<i_file;j_file++) k+=n_of_type_file_rank[j_file][i];
              for(j=0;j<n_of_type_file_rank[i_file][i];j+=n_buffer,k+=n_buffer){
                n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file_rank[i_file][i]-j);
                for(jj=0;jj<n_buffer;jj++)
                  ((float *)buffer)[jj]=
                    (float)(R1_array[k+jj])/
                    (float)(h_Hubble*h_Hubble*plist->mass_unit/pow(plist->length_unit,3.));
                for(i_x_replicate=0;i_x_replicate<n_x_replicate;i_x_replicate++){
                  for(i_y_replicate=0;i_y_replicate<n_y_replicate;i_y_replicate++){
                    for(i_z_replicate=0;i_z_replicate<n_z_replicate;i_z_replicate++){
                      fwrite(buffer,sizeof(float),n_buffer,fp);
                      i_p+=n_buffer;
                    }
                  }
                }
              }
            }
            // ... else write zeros
            else{
              for(j=0;j<n_of_type_file_rank[i_file][i];j+=n_buffer){
                n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file_rank[i_file][i]-j);
                for(jj=0;jj<n_buffer;jj++)
                  ((float *)buffer)[jj]=0.;
                for(i_x_replicate=0;i_x_replicate<n_x_replicate;i_x_replicate++){
                  for(i_y_replicate=0;i_y_replicate<n_y_replicate;i_y_replicate++){
                    for(i_z_replicate=0;i_z_replicate<n_z_replicate;i_z_replicate++){
                      fwrite(buffer,sizeof(float),n_buffer,fp);
                      i_p+=n_buffer;
                    }
                  }
                }
              }
            }
            fclose(fp);
          }
          SID_Barrier(SID.COMM_WORLD);
        }
        if(SID.I_am_Master && flag_recordlength_written){
          fp=fopen(filename,"a");
          fwrite(&record_length,4,1,fp);
          fclose(fp);
        }
      }

      // Smoothing lengths
      for(i_file=0,i_p=0;i_file<n_files;i_file++){
        if(n_files>1)
          sprintf(filename,"%s.%d",filename_in,i_file);
        else
          sprintf(filename,"%s",filename_in);
        record_length=n_replicate*n_particles_file[i_file]*sizeof(float);
        for(i_rank=0,flag_recordlength_written=FALSE;i_rank<SID.n_proc;i_rank++){
          if(SID.My_rank==i_rank){
            fp=fopen(filename,"a");
            if(SID.I_am_Master && !flag_recordlength_written){
              fwrite(&record_length,4,1,fp);
              flag_recordlength_written=TRUE;
            }
            // If smoothing lengths are defined, then write them ...
            if(ADaPS_exist(plist->data,"r_smooth_%s",pname[i])){
              R1_array=(REAL *)ADaPS_fetch(plist->data,"r_smooth_%s",pname[i]);
              for(k=0,j_file=0;j_file<i_file;j_file++) k+=n_of_type_file_rank[j_file][i];
              for(j=0;j<n_of_type_file_rank[i_file][i];j+=n_buffer,k+=n_buffer){
                n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file_rank[i_file][i]-j);
                for(jj=0;jj<n_buffer;jj++)
                  ((float *)buffer)[jj]=
                    (float)(R1_array[k+jj])/
                    (float)(h_Hubble*h_Hubble*plist->mass_unit/pow(plist->length_unit,3.));
                for(i_x_replicate=0;i_x_replicate<n_x_replicate;i_x_replicate++){
                  for(i_y_replicate=0;i_y_replicate<n_y_replicate;i_y_replicate++){
                    for(i_z_replicate=0;i_z_replicate<n_z_replicate;i_z_replicate++){
                      fwrite(buffer,sizeof(float),n_buffer,fp);
                      i_p+=n_buffer;
                    }
                  }
                }
              }
            }
            // ... else write zeros
            else{
              for(j=0;j<n_of_type_file_rank[i_file][i];j+=n_buffer){
                n_buffer=MIN(GADGET_BUFFER_SIZE,n_of_type_file_rank[i_file][i]-j);
                for(jj=0;jj<n_buffer;jj++)
                  ((float *)buffer)[jj]=0.;
                for(i_x_replicate=0;i_x_replicate<n_x_replicate;i_x_replicate++){
                  for(i_y_replicate=0;i_y_replicate<n_y_replicate;i_y_replicate++){
                    for(i_z_replicate=0;i_z_replicate<n_z_replicate;i_z_replicate++){
                      fwrite(buffer,sizeof(float),n_buffer,fp);
                      i_p+=n_buffer;
                    }
                  }
                }
              }
            }
            fclose(fp);
          }
          SID_Barrier(SID.COMM_WORLD);
        }
        if(SID.I_am_Master && flag_recordlength_written){
          fp=fopen(filename,"a");
          fwrite(&record_length,4,1,fp);
          fclose(fp);
        }
      }
      SID_log("Done.",SID_LOG_CLOSE);      
    }
    SID_free((void **)&buffer);
  }

  // Clean-up
  SID_log("Cleaning-up...",SID_LOG_OPEN);
  SID_free((void **)&i_p_hi_file);
  for(i_file=0;i_file<n_files;i_file++){
    SID_free((void **)&n_of_type_file[i_file]);
    SID_free((void **)&n_of_type_file_rank[i_file]);
  }
  SID_free((void **)&n_of_type_file);
  SID_free((void **)&n_of_type_file_rank);
  SID_free((void **)&n_particles_file);
  SID_free((void **)&n_gas_file);
  SID_free((void **)&n_particles_multimass_file);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_log("Done.",SID_LOG_CLOSE);      
}

