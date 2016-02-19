#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpSPH.h>

void process_gadget_file(const char *status_message,
                         char   *filename_read_root,
                         int     snapshot_number,
                         int     select_function(gadget_read_info *fp_gadget,
                                                 void             *params,
                                                 size_t            i_particle,
                                                 size_t            i_particle_type,
                                                 int               i_type,
                                                 GBPREAL          *pos,
                                                 GBPREAL          *vel,
                                                 size_t            ID_i),
                         void    action_function(gadget_read_info *fp_gadget,
                                                 void             *params,
                                                 size_t            i_particle,
                                                 size_t            i_particle_type,
                                                 int               i_type,
                                                 GBPREAL          *pos,
                                                 GBPREAL          *vel,
                                                 size_t            ID_i),
                         void   *params,
                         size_t *n_particles_type_local_pass,
                         size_t *n_particles_type_pass,
                         int    *flag_long_IDs,
                         int     mode){
   SID_log(status_message,SID_LOG_OPEN|SID_LOG_TIMER);
 
   // Initialize particle counters
   size_t *n_particles_type_local;
   size_t *n_particles_type;
   size_t  n_particles_local=0;
   int flag_n_particles_type_local_allocate=FALSE;
   if(n_particles_type_local_pass==NULL){
      n_particles_type_local=(size_t *)SID_malloc(sizeof(size_t)*N_GADGET_TYPE);
      flag_n_particles_type_local_allocate=TRUE;
   }
   else
      n_particles_type_local=n_particles_type_local_pass;
   int flag_n_particles_type_allocate=FALSE;
   if(n_particles_type_pass==NULL){
      n_particles_type=(size_t *)SID_malloc(sizeof(size_t)*N_GADGET_TYPE);
      flag_n_particles_type_allocate=TRUE;
   }
   else
      n_particles_type=n_particles_type_pass;
   for(int i_type=0;i_type<N_GADGET_TYPE;i_type++){
      n_particles_type_local[i_type]=0;
      n_particles_type[i_type]      =0;
   }
 
   // Read the header and determine the input file-format
   gadget_read_info   fp_gadget;
   int                flag_filefound=init_gadget_read(filename_read_root,snapshot_number,&fp_gadget);
   int                flag_multifile=fp_gadget.flag_multifile;
   int                flag_file_type=fp_gadget.flag_file_type;
   gadget_header_info header        =fp_gadget.header;
   if(!flag_filefound)
      SID_trap_error("File not found.",ERROR_LOGIC);

   // Create a count of all the particles in the snapshot
   size_t  n_particles_snap=fp_gadget.n_particles;
 
   // Decide which files each rank will process
   int i_file;
   int i_file_skip;
   int flag_all_to_all;
   int flag_I_am_reader;
   int all_to_all_reader_rank=MASTER_RANK;
   if(check_mode_for_flag(mode,PROCESS_GADGET_BINARY_ALL_TO_ALL)){
      i_file     =0;
      i_file_skip=1;
      flag_all_to_all=TRUE;
      if(SID.My_rank==all_to_all_reader_rank)
         flag_I_am_reader=TRUE;
      else
         flag_I_am_reader=FALSE;
   }
   else{
      i_file     =SID.My_rank;
      i_file_skip=SID.n_proc;
      flag_all_to_all =FALSE;
      flag_I_am_reader=TRUE;
   }

   // Allocate buffers
   GBPREAL   *pos_buffer=(GBPREAL *)SID_malloc(3*sizeof(GBPREAL)*PROCESS_GADGET_FILE_BUFFER_SIZE);
   GBPREAL   *vel_buffer=(GBPREAL *)SID_malloc(3*sizeof(GBPREAL)*PROCESS_GADGET_FILE_BUFFER_SIZE);
   char      *ids_buffer=(char    *)SID_malloc(  sizeof(size_t) *PROCESS_GADGET_FILE_BUFFER_SIZE);
   uint64_t  *ids_buffer_long=(uint64_t *)ids_buffer;
   uint32_t  *ids_buffer_int =(uint32_t *)ids_buffer;
 
   // Loop over the files
   pcounter_info pcounter;
   int           flag_active                =TRUE;
   int           n_ranks_active             =SID.n_proc;
   size_t        n_particles_processed_local=0;
   size_t        n_particles_processed      =0;
   SID_init_pcounter(&pcounter,n_particles_snap,10);
   for(int j_file=0;i_file<fp_gadget.header.n_files;i_file+=i_file_skip,j_file++){
 
      // Open file pointers
      size_t record_length;
      int    record_length_open;
      int    record_length_close;
      char   filename[MAX_FILENAME_LENGTH];
      set_gadget_filename(&fp_gadget,i_file,filename);
      FILE *fp_pos=fopen(filename,"r");
      FILE *fp_vel=fopen(filename,"r");
      FILE *fp_ids=fopen(filename,"r");
 
      // Initialize positions pointer and read header
      record_length=sizeof(gadget_header_info);
      fread_verify(&record_length_open,4,1,fp_pos);
      fread_verify(&header,sizeof(gadget_header_info),1,fp_pos);
      fread_verify(&record_length_close,4,1,fp_pos);
      //if((size_t)record_length_open!=record_length_close)
      //  SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
      //if((size_t)record_length_open!=sizeof(gadget_header_info))
      //  SID_log_warning("Problem with GADGET header size",ERROR_LOGIC);
      fseeko(fp_pos,(off_t)(4),SEEK_CUR);
      fseeko(fp_vel,(off_t)(2*4+sizeof(gadget_header_info)),SEEK_SET);
      fseeko(fp_ids,(off_t)(2*4+sizeof(gadget_header_info)),SEEK_SET);

      // Count the number of particles in the file
      size_t n_particles_file=0;
      for(int i_type=0;i_type<N_GADGET_TYPE;i_type++)
         n_particles_file+=(size_t)(header.n_file[i_type]);

      // Initialize velocities pointer
      record_length=3*n_particles_file*sizeof(GBPREAL);
      fread_verify(&record_length_open,4,1,fp_vel);
      fseeko(fp_vel,(off_t)(record_length),SEEK_CUR);
      fread_verify(&record_length_close,4,1,fp_vel);
      //if(record_length_open!=record_length_close)
      //  SID_log_warning("Problem with GADGET record size (close of positions)",ERROR_LOGIC);
      fseeko(fp_vel,(off_t)(4),SEEK_CUR);
      fseeko(fp_ids,(off_t)(2*4+record_length),SEEK_CUR);

      // Initialize IDs pointer
      record_length=3*n_particles_file*sizeof(GBPREAL);
      fread_verify(&record_length_open,4,1,fp_ids);
      fseeko(fp_ids,(off_t)(record_length),SEEK_CUR);
      fread_verify(&record_length_close,4,1,fp_ids);
      //if(record_length_open!=record_length_close)
      //  SID_log_warning("Problem with GADGET record size (close of velocities)",ERROR_LOGIC);
      fread_verify(&record_length_open,4,1,fp_ids);

      // Determine what kind of IDs we have
      size_t ID_byte_size;
      int    flag_force_long_IDs=FALSE;
      if(j_file==0){
         if((size_t)record_length_open==(n_particles_file*sizeof(uint64_t)) || flag_force_long_IDs){
            (*flag_long_IDs)=TRUE;
            ID_byte_size    =sizeof(uint64_t);
            SID_log("Assuming long IDs.",SID_LOG_COMMENT);
         }
         else if((size_t)record_length_open==(n_particles_file*sizeof(int))){
            (*flag_long_IDs)=FALSE;
            ID_byte_size    =sizeof(int);
            SID_log("Assuming integer IDs.",SID_LOG_COMMENT);
         }
         else{
            (*flag_long_IDs)=TRUE;
            ID_byte_size    =sizeof(uint64_t);
            SID_log_warning("Header does not match ID block size.  Could not determine the IDs type.  Assuming long IDs.",ERROR_LOGIC);
         }
      }
 
      // Process the file in buffered chunks
      int    i_type;
      size_t i_particle;
      size_t j_particle;
      size_t n_buffer=0;
      for(i_type=0,i_particle=0;i_type<N_GADGET_TYPE;i_type++){
         for(size_t j_particle=0;j_particle<header.n_file[i_type];i_particle+=n_buffer,j_particle+=n_buffer){
            n_buffer=MIN(PROCESS_GADGET_FILE_BUFFER_SIZE,header.n_file[i_type]-j_particle);
            if(flag_I_am_reader){
               fread_verify(pos_buffer,sizeof(GBPREAL),3*n_buffer,fp_pos);
               fread_verify(vel_buffer,sizeof(GBPREAL),3*n_buffer,fp_vel);
               fread_verify(ids_buffer,ID_byte_size,     n_buffer,fp_ids);
            }
            if(flag_all_to_all){
               SID_Bcast(pos_buffer,sizeof(GBPREAL)*3*n_buffer,all_to_all_reader_rank,SID.COMM_WORLD);
               SID_Bcast(vel_buffer,sizeof(GBPREAL)*3*n_buffer,all_to_all_reader_rank,SID.COMM_WORLD);
               SID_Bcast(ids_buffer,ID_byte_size     *n_buffer,all_to_all_reader_rank,SID.COMM_WORLD);
            }
            for(int i_buffer=0;i_buffer<n_buffer;i_buffer++){
               size_t ID_i;
               switch((*flag_long_IDs)){
                  case TRUE:
                     ID_i=(size_t)(ids_buffer_long[i_buffer]);
                     break;
                  case FALSE:
                     ID_i=(size_t)(ids_buffer_int[i_buffer]);
                     break;
               }
               if(select_function(&fp_gadget,
                                  params,
                                  n_particles_local,
                                  n_particles_type_local[i_type],
                                  i_type,
                                  &(pos_buffer[3*i_buffer]),
                                  &(vel_buffer[3*i_buffer]),
                                  ID_i)){
                  action_function(&fp_gadget,
                                  params,
                                  n_particles_local,
                                  n_particles_type_local[i_type],
                                  i_type,
                                  &(pos_buffer[3*i_buffer]),
                                  &(vel_buffer[3*i_buffer]),
                                  ID_i);
                  n_particles_type_local[i_type]++;
                  n_particles_local++;
                  fp_gadget.first_action_call=FALSE;
               }
               fp_gadget.first_select_call=FALSE;
            }
            n_particles_processed_local+=n_buffer;
            SID_Allreduce(&flag_active,&n_ranks_active,1,SID_INT,SID_SUM,SID.COMM_WORLD);
            SID_Allreduce(&n_particles_processed_local,&n_particles_processed,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
            SID_check_pcounter(&pcounter,n_particles_processed);
         }
      }
      fclose(fp_pos);
      fclose(fp_vel);
      fclose(fp_ids);
   }
 
   // Ranks that are finished still need to communicate
   //   particle counts and update the progress counter
   //   until all ranks are done
   flag_active=FALSE;
   do{
      SID_Allreduce(&flag_active,&n_ranks_active,1,SID_INT,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(&n_particles_processed_local,&n_particles_processed,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
      SID_check_pcounter(&pcounter,n_particles_processed);
   } while(n_ranks_active>0);

   // Create global counts
   SID_Allreduce(n_particles_type_local,n_particles_type,N_GADGET_TYPE,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);

   SID_log("%lld particles processed.",SID_LOG_COMMENT,n_particles_processed);
 
   // Clean-up
   if(flag_n_particles_type_local_allocate)
      SID_free(SID_FARG n_particles_type_local);
   if(flag_n_particles_type_allocate)
      SID_free(SID_FARG n_particles_type);
   SID_free(SID_FARG pos_buffer);
   SID_free(SID_FARG vel_buffer);
   SID_free(SID_FARG ids_buffer);
 
   SID_log("Done.",SID_LOG_CLOSE);
}

