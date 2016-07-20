#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <gbpLib.h>
#include <gbpSPH.h>

#define MAX_MARK_BUFFER_SIZE SIZE_OF_MEGABYTE

void read_mark_file(plist_info *plist,
                    const char *mark_name,
                    const char *filename_in,
                    int         mode){
  int      i_species;
  size_t   i_particle;
  size_t   j_particle;
  size_t   k_particle;
  int      i_rank;
  size_t   i_mark;
  size_t   n_particles_local;
  size_t  *mark_list_buffer;
  int     *mark_list;
  size_t  *ids_local;
  size_t  *mark_list_local;
  size_t   n_mark_total;
  size_t   n_mark_total_check;
  size_t   n_mark_type_local[N_GADGET_TYPE];
  size_t   n_mark_local;
  size_t   n_particle_local;
  SID_fp   fp_mark_file;
  size_t   i_start_local[N_GADGET_TYPE];
  size_t   n_mark_bcast;
  size_t  *ids_local_index;
  size_t   n_buffer;
  int      flag_allocate;
  int      flag_read_mode;
  int      flag_mark_mode;
  int      flag_op_mode;
  markfile_header_info header={N_GADGET_TYPE};
  
  SID_log("Reading mark file...",SID_LOG_OPEN);
  
  // Interpret run mode
  if(check_mode_for_flag(mode,MARK_READ_ALL))
    flag_read_mode=MARK_READ_ALL;
  else
    flag_read_mode=MARK_DEFAULT;
  if(check_mode_for_flag(mode,MARK_LIST_ONLY))
    flag_mark_mode=MARK_LIST_ONLY;
  else
    flag_mark_mode=MARK_DEFAULT;
  if(check_mode_for_flag(mode,MARK_INIT) || check_mode_for_flag(mode,MARK_OR))
    flag_op_mode=MARK_DEFAULT;
  else
    flag_op_mode=MARK_AND;
  
  // Open mark list and read header
  SID_fopen_chunked(filename_in,
                    "r",
                    &fp_mark_file,
                    &header);
  if(header.n_type!=N_GADGET_TYPE)
    SID_trap_error("Inconsistant number of species in mark file (ie. %d!=%d)!",ERROR_LOGIC,header.n_type,N_GADGET_TYPE);

  // List numbers of particles in the log output
  size_t n_particles_all;
  int    n_non_zero;
  for(i_species=0,n_particles_all=0,n_non_zero=0;i_species<header.n_type;i_species++){
    if(header.n_mark_species[i_species]>0){
      n_particles_all+=header.n_mark_species[i_species];
      n_non_zero++;
    }
  }
  SID_log("%lld",SID_LOG_CONTINUE,n_particles_all);
  if(n_non_zero>0)
    SID_log(" (",SID_LOG_CONTINUE,n_particles_all);
  for(i_species=0;i_species<N_GADGET_TYPE;i_species++){
    if(header.n_mark_species[i_species]>0){
      if(i_species==n_non_zero-1){
        if(n_non_zero>1)
          SID_log("and %lld %s",SID_LOG_CONTINUE,header.n_mark_species[i_species],plist->species[i_species]);
        else
          SID_log("%lld %s",SID_LOG_CONTINUE,header.n_mark_species[i_species],plist->species[i_species]);
      }
      else{
        if(n_non_zero>1)
          SID_log("%lld %s, ",SID_LOG_CONTINUE,header.n_mark_species[i_species],plist->species[i_species]);
        else
          SID_log("%lld %s",SID_LOG_CONTINUE,header.n_mark_species[i_species],plist->species[i_species]);
      }
    }
  }
  if(n_non_zero>0)
    SID_log(") particles...",SID_LOG_CONTINUE);
  else
    SID_log(" particles...",SID_LOG_CONTINUE);


  // Set list sizes and prep offsets for reading
  for(i_species=0,n_mark_local=0,n_mark_total_check=0;i_species<header.n_type;i_species++){
    if(header.n_mark_species[i_species]>0){
      ADaPS_store(&(plist->data),(void *)(&(header.n_mark_species[i_species])),"n_%s_%s",ADaPS_SCALAR_SIZE_T,mark_name,plist->species[i_species]);
      switch(flag_read_mode){
      case MARK_READ_ALL:
        n_mark_type_local[i_species]=header.n_mark_species[i_species];
        i_start_local[i_species]    =0;
        break;
      default:
        n_mark_type_local[i_species]=header.n_mark_species[i_species]/SID.n_proc;
        i_start_local[i_species]    =(SID.My_rank)*n_mark_type_local[i_species];
        if(SID.I_am_last_rank)
          n_mark_type_local[i_species]=header.n_mark_species[i_species]-i_start_local[i_species];
        break;
      }
      ADaPS_store(&(plist->data),(void *)(&(n_mark_type_local[i_species])),"n_local_%s_%s",ADaPS_SCALAR_SIZE_T,mark_name,plist->species[i_species]);
      n_mark_local      +=n_mark_type_local[i_species];
      n_mark_total_check+=header.n_mark_species[i_species];
    }
  }

  // Sanity check
  SID_Allreduce(&n_mark_local,&n_mark_total,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
  if(n_mark_total!=n_mark_total_check)
    SID_trap_error("Particle numbers don't add-up right in read_mark_file!",ERROR_LOGIC);

  // Read file and create/store mark arrays
  switch(flag_mark_mode){
  case MARK_LIST_ONLY:
    for(i_species=0;i_species<header.n_type;i_species++){
      if(header.n_mark_species[i_species]>0){
        // Allocate array
        if(n_mark_type_local[i_species]>0)
          mark_list_local=(size_t *)SID_malloc(sizeof(size_t)*n_mark_type_local[i_species]);
        else
          mark_list_local=NULL;

        // Perform read
        SID_fread_chunked(mark_list_local,
                          n_mark_type_local[i_species],
                          i_start_local[i_species],
                          &fp_mark_file);

        // Sort marked particles
        if(n_mark_type_local[i_species]>0){
          merge_sort(mark_list_local,n_mark_type_local[i_species],NULL,SID_SIZE_T,SORT_INPLACE_ONLY,SORT_COMPUTE_INPLACE);
          ADaPS_store(&(plist->data),(void *)(mark_list_local),"%s_%s",ADaPS_DEFAULT,mark_name,plist->species[i_species]);
        }
      }
    }
    break;
  default:
    mark_list_buffer=(size_t *)SID_malloc(sizeof(size_t)*MAX_MARK_BUFFER_SIZE);
    for(i_species=0;i_species<header.n_type;i_species++){
      if(header.n_mark_species[i_species]>0){
        n_particles_local=((size_t *)ADaPS_fetch(plist->data,"n_%s",plist->species[i_species]))[0];
        // Initialize arrays
        ids_local=(size_t *)ADaPS_fetch(plist->data,"id_%s",plist->species[i_species]);
        if(ADaPS_exist(plist->data,"%s_%s",mark_name,plist->species[i_species])){
          mark_list=(int *)ADaPS_fetch(plist->data,"%s_%s",mark_name,plist->species[i_species]);
          flag_allocate=FALSE;
        }
        else{
          mark_list=(int *)SID_malloc(sizeof(int)*n_particles_local);
          for(i_particle=0;i_particle<n_particles_local;i_particle++)
            mark_list[i_particle]=FALSE;
          flag_allocate=TRUE;
        }      
        merge_sort(ids_local,n_particles_local,&ids_local_index,SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
        // Use a buffer to increase speed
        for(i_particle=0;i_particle<header.n_mark_species[i_species];){
          n_buffer=MIN(header.n_mark_species[i_species]-i_particle,MAX_MARK_BUFFER_SIZE);
          SID_fread_chunked_all(mark_list_local,
                                n_buffer,
                                &fp_mark_file);
          merge_sort(mark_list_local,n_buffer,NULL,SID_SIZE_T,SORT_INPLACE_ONLY,SORT_COMPUTE_INPLACE);
          for(j_particle=0,k_particle=find_index(ids_local,mark_list_buffer[0],n_particles_local,ids_local_index);
              j_particle<n_buffer;
              j_particle++,i_particle++){
            while(ids_local[ids_local_index[k_particle]]<mark_list_local[j_particle] && k_particle<n_buffer-1) k_particle++;
            if(ids_local[ids_local_index[k_particle]]==mark_list_local[j_particle]){
              switch(flag_op_mode){
              case MARK_INIT:
              case MARK_AND:
              case MARK_OR:
                mark_list[i_particle]=TRUE;
                break;
              }
            }
          }
        }
        SID_free((void **)&ids_local_index);
        ADaPS_store(&(plist->data),(void *)mark_list,"%s_%s",ADaPS_DEFAULT,mark_name,plist->species[i_species]);
      }
    }
    SID_free((void **)&mark_list_buffer);
    break;
  }
  SID_fclose_chunked(&fp_mark_file);

  SID_log("Done.",SID_LOG_CLOSE);
}
