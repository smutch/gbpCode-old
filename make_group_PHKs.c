#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpHalos.h>

#define N_BITS_MIN 1

#define GADGET_BUFFER_SIZE_LOCAL  128*SIZE_OF_MEGABYTE

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
  size_t   *buffer_index;
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
  int                read_rank;
  
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

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

    // Number of particles for each species in this file
    fread(&header,sizeof(gadget_header_info),1,fp); 
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
    for(i=0;i<N_GADGET_TYPE;i++){
      n_of_type[i]     =0;
      n_of_type_rank[i]=0;
    }
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
        x_array[i] =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*(size_t)n_of_type_rank[i]);
        y_array[i] =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*(size_t)n_of_type_rank[i]);
        z_array[i] =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*(size_t)n_of_type_rank[i]);
        id_array[i]=(size_t  *)SID_malloc(sizeof(size_t) *(size_t)n_of_type_rank[i]);
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
      SID_Bcast(&header,sizeof(gadget_header_info),read_rank,SID.COMM_WORLD);
      for(i=0,n_particles_file=0;i<N_GADGET_TYPE;i++)
        n_particles_file+=(size_t)header.n_file[i];

      // Skip to IDs
      if(SID.My_rank==read_rank){
        // Skip positions
        fread(&record_length_open,4,1,fp);
        fseeko(fp,record_length_open,SEEK_CUR);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of positions)",ERROR_LOGIC);
        // Skip velocities
        fread(&record_length_open,4,1,fp);
        fseeko(fp,record_length_open,SEEK_CUR);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of velocities)",ERROR_LOGIC);
      }

      // The largest block will always be the positions/velocities.  Use
      //   it's size here to allocate the buffer.
      SID_Bcast(&record_length_open,sizeof(int),read_rank,SID.COMM_WORLD);

      // Allocate buffers
      buffer=SID_malloc(record_length_open);
      keep  =(char *)SID_malloc(sizeof(char)*n_particles_file);
      for(jj=0;jj<n_particles_file;jj++)
        keep[jj]=FALSE;
      
      // Read IDs
      SID_log("Reading IDs...",SID_LOG_OPEN|SID_LOG_TIMER);
      if(SID.My_rank==read_rank){
        fread(&record_length_open,4,1,fp);
        fread(buffer,record_length_open,1,fp);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of ids)",ERROR_LOGIC);
      }

      // Decide what kind of IDs we have
      SID_Bcast(buffer,record_length_open,read_rank,SID.COMM_WORLD);
      if(record_length_open/n_particles_file==sizeof(long long)){
        SID_log("(long long)...",SID_LOG_CONTINUE);
        flag_LONGIDS=TRUE;
      }
      else{
        SID_log("(int)...",SID_LOG_CONTINUE);
        flag_LONGIDS=FALSE;
      }

      // Process IDs
      for(i=0,jj=0;i<N_GADGET_TYPE;i++){
        n_keep[i]=0;
        if(header.n_file[i]>0){
          if(flag_LONGIDS){
            buffer_i=&(((long long *)buffer)[jj]);
            merge_sort(buffer_i,(size_t)header.n_file[i],&buffer_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
          }
          else{
            buffer_i=&(((int *)buffer)[jj]);
            merge_sort(buffer_i,(size_t)header.n_file[i],&buffer_index,SID_INT,   SORT_COMPUTE_INDEX,FALSE);
          }
          for(j=0,k=0,l=jj;j<header.n_file[i] && k<n_particles_rank;j++,jj++){
            switch(flag_LONGIDS){
            case TRUE:
              id_test=(size_t)(((long long *)buffer_i)[buffer_index[j]]);
              break;
            case FALSE:
              id_test=(size_t)(((int       *)buffer_i)[buffer_index[j]]);
              break;
            }
            while(id_list[id_list_index[k]]<id_test && k<(n_particles_rank-1))
              k++;
            if(id_list[id_list_index[k]]==id_test){
              keep[l+buffer_index[j]]=TRUE;
              n_keep[i]++;
              k++;
            }
          }
          SID_free(SID_FARG buffer_index);
        }
        n_particles_kept+=n_keep[i];
      }
      if(flag_LONGIDS){
        for(i=0,jj=0;i<N_GADGET_TYPE;i++){
          for(j=0,k=0;j<header.n_file[i];j++,jj++){
            if(keep[jj]){
              id_test=(size_t)((long long *)buffer)[jj];
              id_array[i][k+k_offset[i]]=id_test;
              k++;
            }
          }
          if(k!=n_keep[i])
            SID_trap_error("Particle count mismatch (ie. %d!=%d) during IDs read",ERROR_LOGIC,k,n_keep[i]);
        }
      }
      else{
        for(i=0,jj=0;i<N_GADGET_TYPE;i++){
          for(j=0,k=0;j<header.n_file[i];j++,jj++){
            if(keep[jj]){
              id_test=(size_t)((int *)buffer)[jj];
              id_array[i][k+k_offset[i]]=id_test;
              k++;
            }
          }
          if(k!=n_keep[i])
            SID_trap_error("Particle count mismatch (ie. %d!=%d) during IDs read",ERROR_LOGIC,k,n_keep[i]);
        }
      }
      SID_log("Done.",SID_LOG_CLOSE);

      // Skip back to positions
      if(SID.My_rank==read_rank){
        rewind(fp);
        // Skip header
        fread(&record_length_open,4,1,fp);
        fseeko(fp,record_length_open,SEEK_CUR);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (skip of header)",ERROR_LOGIC);
      }

      // Read positions
      SID_log("Reading positions...",SID_LOG_OPEN|SID_LOG_TIMER);
      if(SID.My_rank==read_rank){
        fread(&record_length_open,4,1,fp);
        fread(buffer,record_length_open,1,fp);
        fread(&record_length_close,4,1,fp);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of positions) (ie. %d!=%d)",ERROR_LOGIC,record_length_open,record_length_close);
      }
      SID_Bcast(&record_length_open,sizeof(int),read_rank,SID.COMM_WORLD);
      SID_Bcast(buffer,record_length_open,      read_rank,SID.COMM_WORLD);
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
    SID_free(SID_FARG id_list_index);
    
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

    //  ... ids ...
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_of_type_rank[i]>0)
        ADaPS_store(&(plist->data),(void *)id_array[i],"id_%s",ADaPS_DEFAULT,pname[i]);
    }

    SID_log("Done.",SID_LOG_CLOSE);
  }
  else
    SID_trap_error("Could not find file with root {%s}",ERROR_IO_OPEN,filename_root_in);
}


int main(int argc, char *argv[]){
  plist_info  plist;
  char        filename_PHKs_root[256];
  char        filename_snapshot_root[256];
  char        filename_snapshot[256];
  char        filename_output_properties_dir[256];
  char        filename_output_properties[256];
  char        filename_output_profiles_dir[256];
  char        filename_output_profiles[256];
  char        filename_output_properties_temp[256];
  char        filename_output_profiles_temp[256];
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
  size_t      n_particles_cumulative;
  GBPREAL    *x_array;
  GBPREAL    *y_array;
  GBPREAL    *z_array;
  size_t     *ids_particles;
  size_t     *ids_particles_index;
  size_t     *ids_groups;
  int        *n_particles_groups_process;
  int        *n_particles_groups;
  int        *n_particles_subgroups;
  int        *group_offset;
  size_t      n_particles_in_groups;
  size_t     *ids_snapshot;
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
  double      box_size;
  double      particle_mass;
  
  int         r_val;
  struct      stat file_stats;
  size_t      n_bytes;
  size_t      n_bytes_buffer;
  void       *buffer;
  
  FILE       *fp_PHKs;
  FILE       *fp_profiles;
  FILE       *fp_PHKs_temp;
  FILE       *fp_profiles_temp;
  cosmo_info *cosmo;
  halo_properties_info  properties;
  halo_profile_info     profile;
  int                   n_temp;
  int                   n_truncated;
  int                   largest_truncated;
  int                   largest_truncated_local;
  int                   n_bits=N_BITS_MIN;

  SID_init(&argc,&argv,NULL);
  SID_profile_start("make_group_PHKs",SID_PROFILE_NOTMPIENABLED);

  // Fetch user inputs
  strcpy(filename_snapshot_root,argv[1]);
  strcpy(filename_PHKs_root,    argv[2]);
  i_file_lo  =atoi(argv[3]);
  i_file_hi  =atoi(argv[4]);
  i_file_skip=atoi(argv[5]);

  SID_log("Generating group PH keys for files #%d->#%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file_lo,i_file_hi);

  char *filename_number;
  for(i_file=i_file_lo;i_file<=i_file_hi;i_file+=i_file_skip){
    SID_log("Processing file #%03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file);
    //SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);

    // Read group and particle info
    int    *PHK_group      =NULL;
    size_t *PHK_group_index=NULL;
    init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
    filename_number=(char *)SID_malloc(sizeof(char)*10);
    sprintf(filename_number,"%03d", i_file);
    ADaPS_store(&(plist.data),(void *)filename_number,"read_catalog",ADaPS_DEFAULT);
    read_groups(filename_PHKs_root,i_file,READ_GROUPS_ALL|READ_GROUPS_MBP_IDS_ONLY,&plist,filename_number);
    read_gadget_binary_local(filename_snapshot_root,i_file,&plist);

    // Fetch some header information
    n_particles_cumulative = ((size_t  *)ADaPS_fetch(plist.data,"n_particles_all_%s",filename_number))[0];
    n_groups_all           = ((int     *)ADaPS_fetch(plist.data,"n_groups_all_%s",   filename_number))[0];
    n_groups               = ((int     *)ADaPS_fetch(plist.data,"n_groups_%s",       filename_number))[0];
    box_size               = ((double  *)ADaPS_fetch(plist.data,"box_size",          filename_number))[0];

    // Determine the number of bits to use for the PHKs
    double dx;
    dx=10.*M_PER_MPC;
    for(n_bits=N_BITS_MIN;(box_size/pow(2.,(double)(n_bits+1)))>dx && n_bits<=20;) n_bits++;
 
    // If there's any groups to analyze ...
    if(n_groups>0){
       // Fetch some needed data
       n_particles_groups =  (int     *)ADaPS_fetch(plist.data,"n_particles_group_%s",filename_number);
       x_array            =  (GBPREAL *)ADaPS_fetch(plist.data,"x_dark");
       y_array            =  (GBPREAL *)ADaPS_fetch(plist.data,"y_dark");
       z_array            =  (GBPREAL *)ADaPS_fetch(plist.data,"z_dark");
       ids_particles      =  (size_t  *)ADaPS_fetch(plist.data,"id_dark");
       ids_groups         =  (size_t  *)ADaPS_fetch(plist.data,"particle_ids_%s",filename_number);
       PHK_group          =  (int     *)SID_malloc(sizeof(int)*n_groups);

       // Sort the GADGET arrays by their particle IDs
       merge_sort(ids_particles,n_groups,&ids_particles_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);

       // Compute PHKs
       SID_log("Computing PHKs (using %d bits per dimension)...",SID_LOG_OPEN,n_bits);
       for(i_group=0;i_group<n_groups;i_group++){

         // The particles are not stored in the same order as the groups.  We need
         //    to find the index for the MBP of each particle as a result.
         size_t particle_index;
         particle_index=ids_particles_index[find_index(ids_particles,ids_groups[i_group],(size_t)n_groups,ids_particles_index)];

         // Compute the key for this particle
         PHK_group[i_group]=compute_PHK_from_Cartesian(n_bits,3,(double)x_array[particle_index]/box_size,
                                                                (double)y_array[particle_index]/box_size,
                                                                (double)z_array[particle_index]/box_size);
       }
       SID_free(SID_FARG ids_particles_index);
       SID_log("Done.",SID_LOG_CLOSE);

       // Sort PHKs
       SID_log("Sorting PHKs...",SID_LOG_OPEN);
       merge_sort((void *)PHK_group,n_groups,&PHK_group_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);
       SID_log("Done.",SID_LOG_CLOSE);

    }

    // Count the number of particles
    for(i_group=0,n_particles_cumulative=0;i_group<n_groups;i_group++)
       n_particles_cumulative+=n_particles_groups[PHK_group_index[i_group]];

    // Write results
    int index_temp;
    SID_log("Writing results...",SID_LOG_OPEN);
    sprintf(filename_output_properties,"%s_%s.catalog_PHKs",filename_PHKs_root,filename_number);
    fp_PHKs=fopen(filename_output_properties,"w");
    fwrite(&n_groups,              sizeof(int),   1,fp_PHKs);
    fwrite(&n_bits,                sizeof(int),   1,fp_PHKs);
    fwrite(&n_particles_cumulative,sizeof(size_t),1,fp_PHKs);
    for(i_group=0,n_particles_cumulative=0;i_group<n_groups;i_group++){
       index_temp             =(int)PHK_group_index[i_group];
       n_particles_cumulative+=n_particles_groups[index_temp];
       fwrite(&(PHK_group[index_temp]),sizeof(int),   1,fp_PHKs);
       fwrite(&index_temp,             sizeof(int),   1,fp_PHKs);
       fwrite(&n_particles_cumulative, sizeof(size_t),1,fp_PHKs);
    }
    fclose(fp_PHKs);
    SID_log("Done.",SID_LOG_CLOSE);

    // Clean-up
    free_plist(&plist);
    if(n_groups>0){
       SID_free(SID_FARG PHK_group);
       SID_free(SID_FARG PHK_group_index);
    }

    SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
    SID_log("Done.",SID_LOG_CLOSE);
  }

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

