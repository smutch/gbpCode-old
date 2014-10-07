#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpSPH.h>

#define WRITE_BUFFER_SIZE_LOCAL (1024*1024)

size_t set_n_of_type_file_left_local(size_t n_of_type_all,size_t n_of_type_processed,int i_file,int n_files);
size_t set_n_of_type_file_left_local(size_t n_of_type_all,size_t n_of_type_processed,int i_file,int n_files){
   if(n_files==1)
      return(n_of_type_all);
   else if(i_file==(n_files-1))
      return(n_of_type_all-n_of_type_processed);
   else
      return(n_of_type_all/(size_t)n_files);
}

size_t n_particles_for_file(size_t n_particles_all,
                            int    n_files,
                            int    i_rank,
                            int    n_ranks);
size_t n_particles_for_file(size_t n_particles_all,
                            int    n_files,
                            int    i_rank,
                            int    n_ranks){
   if(n_files==1)
      return(n_particles_all);
}

void write_gadget_binary_block(plist_info  *plist,
                               int          n_files,
                               char        *filename_out_root,
                               int          n_variables,
                               size_t       storage_size[3],
                               char         variable[3][8],
                               int          convert_mode[3],
                               int          mode);
void write_gadget_binary_block(plist_info  *plist,
                               int          n_files,
                               char        *filename_out_root,
                               int          n_variables,
                               size_t       storage_size[3],
                               char         variable[3][8],
                               int          convert_mode[3],
                               int          mode){

  // Write a log message
  if(n_variables==1)
    SID_log("Writing GADGET binary block for variable {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,variable[0]);
  else{
    SID_log("Writing GADGET binary block for vector {",SID_LOG_OPEN|SID_LOG_TIMER,variable[0]);
    for(int i_variable=0;i_variable<n_variables;i_variable++){
       if(i_variable>0)
          SID_log(",",SID_LOG_CONTINUE);
       SID_log("%s",SID_LOG_CONTINUE,variable[i_variable]);
    }
    SID_log("}...",SID_LOG_CONTINUE);
  }

  // Create write buffers
  int     write_size[3];
  char   *local_array[3];
  char   *buffer_alloc[3];
  char   *buffer_alloc_vector;
  size_t  max_variable_size=sizeof(size_t);
  for(int i_coord=0;i_coord<n_variables;i_coord++)
     buffer_alloc[i_coord]=(char *)SID_malloc(max_variable_size*WRITE_BUFFER_SIZE_LOCAL);
  if(n_variables>1)
     buffer_alloc_vector=(char *)SID_malloc(3*max_variable_size*WRITE_BUFFER_SIZE_LOCAL);

  // Write block
  for(int i_type=0;i_type<N_GADGET_TYPE;i_type++){
     if(ADaPS_exist(plist->data,"n_all_%s",plist->species[i_type])){
        size_t  n_of_type_all=((size_t *)ADaPS_fetch(plist->data,"n_all_%s",plist->species[i_type]))[0];
        size_t  n_of_type_processed=0;
        int     i_file=0;
        size_t  n_of_type_file_left=set_n_of_type_file_left_local(n_of_type_all,n_of_type_processed,i_file,n_files);
        int     block_size;
        FILE   *fp_out=NULL;
        for(int i_rank=0;i_rank<SID.n_proc;i_rank++){
           // Fetch and communicate this rank's particle count
           size_t n_of_type    =((size_t *)ADaPS_fetch(plist->data,"n_%s",    plist->species[i_type]))[0];
           SID_Bcast(&n_of_type,sizeof(size_t),i_rank,SID.COMM_WORLD);
           // Fetch the local data pointers
           for(int i_variable=0;i_variable<n_variables;i_variable++)
              local_array[i_variable]=(char *)ADaPS_fetch(plist->data,"%s_%s",variable[i_variable],plist->species[i_type]);

           // We only need the sending and writing (MASTER) rank
           if(SID.My_rank==i_rank || SID.I_am_Master){
              size_t n_buffer;
              for(size_t j_particle=0;j_particle<n_of_type;){
                 // Set the buffer size
                 size_t n_buffer=MIN(WRITE_BUFFER_SIZE_LOCAL,n_of_type-j_particle);
                 if(n_buffer>n_of_type_file_left)
                    n_buffer=n_of_type_file_left;
                 // Communicate the buffer; MASTER->MASTER
                 if(i_rank==MASTER_RANK){
                    for(int i_variable=0;i_variable<n_variables;i_variable++)
                       memcpy(buffer_alloc[i_variable],&(local_array[i_variable][j_particle*storage_size[i_variable]]),n_buffer*storage_size[i_variable]);
                 }
                 // Communicate the buffer; MASTER->OTHER
                 else{
                    for(int i_variable=0;i_variable<n_variables;i_variable++){
                       // Communicate the buffer
                       SID_Sendrecv(&(local_array[i_variable][j_particle*storage_size[i_variable]]),
                                    n_buffer*storage_size[i_variable],
                                    SID_BYTE,
                                    MASTER_RANK,
                                    191273,
                                    buffer_alloc[i_variable],
                                    n_buffer*storage_size[i_variable],
                                    SID_BYTE,
                                    i_rank,
                                    191273,
                                    SID.COMM_WORLD);
                    }
                 }
                 // Write the buffer
                 if(SID.I_am_Master){
                    // Perform any desired conversion and set the variable write size
                    size_t buffer_item_size=0;
                    for(size_t i_variable=0;i_variable<n_variables;i_variable++){
                       // Convert IDs to original byte size
                       if(convert_mode[i_variable]==1){
                          size_t *buffer_size_t=(size_t *)buffer_alloc[i_variable];
                          int    *buffer_int   =(int    *)buffer_alloc[i_variable];
                          for(size_t i_buffer=0;i_buffer<n_buffer;i_buffer++)
                             buffer_int[i_buffer]=(int)(buffer_size_t[i_buffer]);
                          write_size[i_variable]=sizeof(int);
                       }
                       else
                          write_size[i_variable]=storage_size[i_variable];
                       buffer_item_size+=write_size[i_variable];
                    }

                    // Open files and decide how many particles it gets
                    if(n_of_type_file_left==0 || fp_out==NULL){
                       if(fp_out!=NULL){
                          n_of_type_file_left=set_n_of_type_file_left_local(n_of_type_all,n_of_type_processed,i_file,n_files);
                          fwrite(&block_size,sizeof(int),1,fp_out);
                          fclose(fp_out);
                       }
                       char  filename_out[MAX_FILENAME_LENGTH];
                       if(n_files==1)
                          sprintf(filename_out,"%s",filename_out_root);
                       else
                          sprintf(filename_out,"%s.%d",filename_out_root,i_file);
                       i_file++;
                       fp_out=fopen(filename_out,"a");
                       // Write opening length
                       if(n_of_type_processed==0){
                          block_size=n_of_type_file_left*buffer_item_size;
                          fwrite(&block_size,sizeof(int),1,fp_out);
                       }
                    }

                    size_t n_write_file=n_buffer;
                    size_t n_buffer_bytes_written=0;
                    size_t n_bytes_write;
                    for(size_t k_particle=0;
                        k_particle<n_buffer;
                        k_particle+=n_write_file,
                        n_buffer_bytes_written+=n_bytes_write,
                        j_particle+=n_write_file,
                        n_of_type_file_left-=n_write_file,
                        n_of_type_processed+=n_write_file){
                       // Open files and decide how many particles it gets
                       if(n_of_type_file_left==0 || fp_out==NULL){
                          if(fp_out!=NULL){
                             n_of_type_file_left=set_n_of_type_file_left_local(n_of_type_all,n_of_type_processed,i_file,n_files);
                             fwrite(&block_size,sizeof(int),1,fp_out);
                             fclose(fp_out);
                          }
                          char  filename_out[MAX_FILENAME_LENGTH];
                          if(n_files==1)
                             sprintf(filename_out,"%s",filename_out_root);
                          else
                             sprintf(filename_out,"%s.%d",filename_out_root,i_file);
                          i_file++;
                          fp_out=fopen(filename_out,"a");
                          // Write opening length
                          if(n_of_type_processed==0){
                             block_size=n_of_type_file_left*buffer_item_size;
                             fwrite(&block_size,sizeof(int),1,fp_out);
                          }
                       }
                       n_write_file =MIN(n_buffer,n_of_type-j_particle);
                       n_bytes_write=n_write_file*buffer_item_size;
                       
                       // Write vectors
                       if(n_variables>1){
                          size_t i_buffer=0;
                          size_t j_buffer[3]={0,0,0};
                          for(size_t k_buffer=0;k_buffer<n_buffer;k_buffer++){
                             for(size_t i_variable=0;i_variable<n_variables;i_buffer+=write_size[i_variable],j_buffer[i_variable]+=write_size[i_variable],i_variable++)
                                memcpy(&(buffer_alloc_vector[i_buffer]),
                                       &(buffer_alloc[i_variable][j_buffer[i_variable]]),
                                       write_size[i_variable]);
                          }
                          fwrite(&(buffer_alloc_vector[n_buffer_bytes_written]),buffer_item_size,n_write_file,fp_out);
                       }
                       // Write scalars
                       else
                          fwrite(&(buffer_alloc[0][n_buffer_bytes_written]),buffer_item_size,n_write_file,fp_out);
                    } // Loop over buffer
                 } // if I'm master 
              } // loop over n_particles_type
           } // only i_rank and MASTER
           SID_Bcast(&i_file,             sizeof(int),   i_rank,SID.COMM_WORLD);
           SID_Bcast(&n_of_type_processed,sizeof(size_t),i_rank,SID.COMM_WORLD);
           SID_Bcast(&n_of_type_file_left,sizeof(size_t),i_rank,SID.COMM_WORLD);
        } // i_rank
        if(fp_out!=NULL){
           fwrite(&block_size,sizeof(int),1,fp_out);
           fclose(fp_out);
        }
     } // if type is defined
  } // i_type

  // Clean-up
  for(int i_coord=0;i_coord<n_variables;i_coord++)
     SID_free(SID_FARG buffer_alloc[i_coord]);
  if(n_variables>1)
     SID_free(SID_FARG buffer_alloc_vector);

  SID_log("Done.",SID_LOG_CLOSE);
}

void write_gadget_binary_new(plist_info  *plist,
                             char        *filename_out_root,
                             int          n_files,
                             int          mode){

  SID_log("Writing snapshot to {%s} in %d parts...",SID_LOG_OPEN,filename_out_root,n_files);

  // Fetch total particle counts
  size_t n_all[N_GADGET_TYPE];
  for(int i_type=0;i_type<N_GADGET_TYPE;i_type++){
     if(ADaPS_exist(plist->data,"n_all_%s",plist->species[i_type]))
        n_all[i_type]=((size_t *)ADaPS_fetch(plist->data,"n_all_%s",plist->species[i_type]))[0];
     else
        n_all[i_type]=0;
  }

  // Write headers
  gadget_header_info *header;
  header=(gadget_header_info *)SID_calloc(sizeof(gadget_header_info));
  for(int i_type=0;i_type<N_GADGET_TYPE;i_type++){
     int    lo_word=(unsigned int)(((size_t)(n_all[i_type])<<32)>>32);
     int    hi_word=(unsigned int) ((size_t)(n_all[i_type])>>32);
     if(ADaPS_exist(plist->data,"mass_array_%s",plist->species[i_type]))
        header->mass_array[i_type]=((double *)ADaPS_fetch(plist->data,"mass_array_%s",plist->species[i_type]))[0];
     else
        header->mass_array[i_type]=0.;
     header->n_all_lo_word[i_type]=lo_word;
     header->n_all_hi_word[i_type]=hi_word;
  }
  header->n_files        =n_files;
  header->time           =((double *)ADaPS_fetch(plist->data,"time"))[0];
  header->redshift       =((double *)ADaPS_fetch(plist->data,"redshift"))[0];
  header->box_size       =((double *)ADaPS_fetch(plist->data,"box_size"))[0];
  header->Omega_M        =((double *)ADaPS_fetch(plist->data,"Omega_M"))[0];
  header->Omega_Lambda   =((double *)ADaPS_fetch(plist->data,"Omega_Lambda"))[0];
  header->h_Hubble       =((double *)ADaPS_fetch(plist->data,"h_Hubble"))[0];
  header->flag_SFr       =((int    *)ADaPS_fetch(plist->data,"flag_Sfr"))[0];
  header->flag_feedback  =((int    *)ADaPS_fetch(plist->data,"flag_feedback"))[0];
  header->flag_cooling   =((int    *)ADaPS_fetch(plist->data,"flag_cooling"))[0];
  header->flag_ages      =((int    *)ADaPS_fetch(plist->data,"flag_ages"))[0];
  header->flag_metals    =((int    *)ADaPS_fetch(plist->data,"flag_metals"))[0];
  header->flag_entropyICs=((int    *)ADaPS_fetch(plist->data,"flag_entropyICs"))[0];
  for(int i_file=0;i_file<n_files;i_file++){
     // Set particle counts for the file
     for(int i_type=0;i_type<N_GADGET_TYPE;i_type++)
        header->n_file[i_type]=(int)set_n_of_type_file_left_local(n_all[i_type],0,i_file,n_files);

     // Set filename
     char filename_out[MAX_FILENAME_LENGTH];
     if(n_files>1)
        sprintf(filename_out,"%s.%d",filename_out_root,i_file);
     else
        sprintf(filename_out,"%s",filename_out_root);

     // Open file and write header
     int   record_length=sizeof(gadget_header_info);
     FILE *fp_out=fopen(filename_out,"w");
     fwrite(&record_length,sizeof(int),               1,fp_out);
     fwrite(header,        sizeof(gadget_header_info),1,fp_out);
     fwrite(&record_length,sizeof(int),               1,fp_out);
     fclose(fp_out);
  }
  SID_free(SID_FARG header);

  // Initialize some arrays
  size_t  storage_size[3];
  char    variable[3][8];
  int     convert_mode[3];
  int     n_variables;

  // Write positions
  n_variables=3;
  sprintf(variable[0],"x");
  sprintf(variable[1],"y");
  sprintf(variable[2],"z");
  storage_size[0]=sizeof(GBPREAL);
  storage_size[1]=sizeof(GBPREAL);
  storage_size[2]=sizeof(GBPREAL);
  convert_mode[0]=0;
  convert_mode[1]=0;
  convert_mode[2]=0;
  write_gadget_binary_block(plist,
                            n_files,
                            filename_out_root,
                            n_variables,
                            storage_size,
                            variable,
                            convert_mode,
                            mode);

  // Write velocities
  sprintf(variable[0],"vx");
  sprintf(variable[1],"vy");
  sprintf(variable[2],"vz");
  write_gadget_binary_block(plist,
                            n_files,
                            filename_out_root,
                            n_variables,
                            storage_size,
                            variable,
                            convert_mode,
                            mode);

  // Write IDs
  int flag_long_IDs=((int *)ADaPS_fetch(plist->data,"flag_long_IDs"))[0];
  n_variables=1;
  sprintf(variable[0],"id");
  storage_size[0]=sizeof(size_t);
  if(!flag_long_IDs)
     convert_mode[0]=1;
  else
     convert_mode[0]=0;
  write_gadget_binary_block(plist,
                            n_files,
                            filename_out_root,
                            n_variables,
                            storage_size,
                            variable,
                            convert_mode,
                            mode);
      
  SID_log("Done.",SID_LOG_CLOSE);      
}

