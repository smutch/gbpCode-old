#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpSPH.h>

#define READ_BUFFER_SIZE_LOCAL (1024*1024)
#define WRITE_GADGET_BINARY_DEFAULT 0
typedef struct select_gadget_volume_params_info select_gadget_volume_params_info;
struct select_gadget_volume_params_info{
   plist_info *plist;
   GBPREAL     cen[3];
   GBPREAL     size;
   GBPREAL     size2;
   GBPREAL     box_size;
};

void allocate_gadget_particles(plist_info *plist, 
                               size_t     *n_particles_type_local,
                               size_t     *n_particles_type,
                               int         flag_long_IDs);
void allocate_gadget_particles(plist_info *plist, 
                               size_t     *n_particles_type_local,
                               size_t     *n_particles_type,
                               int         flag_long_IDs){
   SID_log("Allocating for particles...",SID_LOG_OPEN);
   for(int i_type=0;i_type<N_GADGET_TYPE;i_type++){
      if(n_particles_type[i_type]>0){
         GBPREAL *x_array =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_particles_type_local[i_type]);
         GBPREAL *y_array =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_particles_type_local[i_type]);
         GBPREAL *z_array =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_particles_type_local[i_type]);
         GBPREAL *vx_array=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_particles_type_local[i_type]);
         GBPREAL *vy_array=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_particles_type_local[i_type]);
         GBPREAL *vz_array=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_particles_type_local[i_type]);
         size_t  *id_array=(size_t  *)SID_malloc(sizeof(size_t) *n_particles_type_local[i_type]);
         ADaPS_store(&(plist->data),(void *)x_array, "x_%s", ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)y_array, "y_%s", ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)z_array, "z_%s", ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)vx_array,"vx_%s",ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)vy_array,"vy_%s",ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)vz_array,"vz_%s",ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)id_array,"id_%s",ADaPS_DEFAULT,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)(&(n_particles_type_local[i_type])),"n_%s",    ADaPS_SCALAR_SIZE_T,plist->species[i_type]);
         ADaPS_store(&(plist->data),(void *)(&(n_particles_type[i_type])),      "n_all_%s",ADaPS_SCALAR_SIZE_T,plist->species[i_type]);
      }
   }
   ADaPS_store(&(plist->data),(void *)(&flag_long_IDs),"flag_long_IDs",ADaPS_SCALAR_INT);
   SID_log("Done.",SID_LOG_CLOSE);
}

int select_gadget_cube(gadget_read_info *fp_gadget,
                       void             *params,
                       size_t            i_particle,
                       size_t            i_particle_type,
                       int               i_type,
                       GBPREAL          *pos,
                       GBPREAL          *vel,
                       size_t            ID_i);
int select_gadget_cube(gadget_read_info *fp_gadget,
                       void             *params,
                       size_t            i_particle,
                       size_t            i_particle_type,
                       int               i_type,
                       GBPREAL          *pos,
                       GBPREAL          *vel,
                       size_t            ID_i){
   static GBPREAL  half_cube_size;
   static GBPREAL  half_box_size;
   static GBPREAL  box_size;
   static GBPREAL *cen;
   if(fp_gadget->first_select_call){
      cen           =    ((select_gadget_volume_params_info *)params)->cen;
      half_cube_size=0.5*((select_gadget_volume_params_info *)params)->size;
      half_box_size =0.5*((select_gadget_volume_params_info *)params)->box_size;
      box_size      =    ((select_gadget_volume_params_info *)params)->box_size;
   }
   int flag_select=TRUE;
   for(int i_coord=0;i_coord<3 && flag_select;i_coord++){
      GBPREAL coord=pos[i_coord]-cen[i_coord];
      force_periodic(&coord,-half_box_size,box_size);
      if(coord<(-half_cube_size) || coord>half_cube_size)
         flag_select=FALSE;
   }
   return(flag_select);
}

int select_gadget_sphere(gadget_read_info *fp_gadget,
                         void             *params,
                         size_t            i_particle,
                         size_t            i_particle_type,
                         int               i_type,
                         GBPREAL          *pos,
                         GBPREAL          *vel,
                         size_t            ID_i);
int select_gadget_sphere(gadget_read_info *fp_gadget,
                         void             *params,
                         size_t            i_particle,
                         size_t            i_particle_type,
                         int               i_type,
                         GBPREAL          *pos,
                         GBPREAL          *vel,
                         size_t            ID_i){
   static GBPREAL  sphere_radius;
   static GBPREAL  sphere_radius2;
   static GBPREAL  half_box_size;
   static GBPREAL  box_size;
   static GBPREAL *cen;
   if(fp_gadget->first_select_call){
      cen           =    ((select_gadget_volume_params_info *)params)->cen;
      sphere_radius =    ((select_gadget_volume_params_info *)params)->size;
      sphere_radius2=    ((select_gadget_volume_params_info *)params)->size2;
      half_box_size =0.5*((select_gadget_volume_params_info *)params)->box_size;
      box_size      =    ((select_gadget_volume_params_info *)params)->box_size;
   }
   int flag_select=TRUE;
   GBPREAL radius2=0;
   for(int i_coord=0;i_coord<3 && flag_select;i_coord++){
      GBPREAL coord=pos[i_coord]-cen[i_coord];
      force_periodic(&coord,-half_box_size,box_size);
      if(coord<(-sphere_radius) || coord>sphere_radius)
         flag_select=FALSE;
      else
         radius2+=(coord*coord);
   }
   if(radius2>sphere_radius2)
      flag_select=FALSE;

   return(flag_select);
}

void count_gadget_particles(gadget_read_info *fp_gadget,
                            void             *params,
                            size_t            i_particle,
                            size_t            i_particle_type,
                            int               i_type,
                            GBPREAL          *pos,
                            GBPREAL          *vel,
                            size_t            ID_i);
void count_gadget_particles(gadget_read_info *fp_gadget,
                            void             *params,
                            size_t            i_particle,
                            size_t            i_particle_type,
                            int               i_type,
                            GBPREAL          *pos,
                            GBPREAL          *vel,
                            size_t            ID_i){
}

void store_gadget_particles(gadget_read_info *fp_gadget,
                            void             *params,
                            size_t            i_particle,
                            size_t            i_particle_type,
                            int               i_type,
                            GBPREAL          *pos,
                            GBPREAL          *vel,
                            size_t            ID_i);
void store_gadget_particles(gadget_read_info *fp_gadget,
                            void             *params,
                            size_t            i_particle,
                            size_t            i_particle_type,
                            int               i_type,
                            GBPREAL          *pos,
                            GBPREAL          *vel,
                            size_t            ID_i){
   static GBPREAL *x_particle[N_GADGET_TYPE];
   static GBPREAL *y_particle[N_GADGET_TYPE];
   static GBPREAL *z_particle[N_GADGET_TYPE];
   static GBPREAL *vx_particle[N_GADGET_TYPE];
   static GBPREAL *vy_particle[N_GADGET_TYPE];
   static GBPREAL *vz_particle[N_GADGET_TYPE];
   static size_t  *id_particle[N_GADGET_TYPE];
   if(fp_gadget->first_action_call){
      plist_info *plist=((select_gadget_volume_params_info *)params)->plist;

      // Store header stuff
      ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.box_size)),        "box_size",        ADaPS_SCALAR_DOUBLE);
      ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.time)),            "time",            ADaPS_SCALAR_DOUBLE);
      ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.redshift)),        "redshift",        ADaPS_SCALAR_DOUBLE);
      ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.flag_SFr)),        "flag_Sfr",        ADaPS_SCALAR_INT);
      ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.flag_feedback)),   "flag_feedback",   ADaPS_SCALAR_INT);
      ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.flag_cooling)),    "flag_cooling",    ADaPS_SCALAR_INT);
      ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.flag_metals)),     "flag_metals",     ADaPS_SCALAR_INT);
      ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.flag_ages)),       "flag_ages",       ADaPS_SCALAR_INT);
      ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.flag_entropyICs)), "flag_entropyICs", ADaPS_SCALAR_INT);
      ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.Omega_M)),         "Omega_M",         ADaPS_SCALAR_DOUBLE);
      ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.Omega_Lambda)),    "Omega_Lambda",    ADaPS_SCALAR_DOUBLE);
      ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.h_Hubble)),        "h_Hubble",        ADaPS_SCALAR_DOUBLE);

      // Fetch particle arrays
      for(int j_type=0;j_type<N_GADGET_TYPE;j_type++){
         if(ADaPS_exist(plist->data,"n_%s",plist->species[j_type])){
            x_particle [j_type]=(GBPREAL *)ADaPS_fetch(plist->data,"x_%s", plist->species[j_type]);
            y_particle [j_type]=(GBPREAL *)ADaPS_fetch(plist->data,"y_%s", plist->species[j_type]);
            z_particle [j_type]=(GBPREAL *)ADaPS_fetch(plist->data,"z_%s", plist->species[j_type]);
            vx_particle[j_type]=(GBPREAL *)ADaPS_fetch(plist->data,"vx_%s",plist->species[j_type]);
            vy_particle[j_type]=(GBPREAL *)ADaPS_fetch(plist->data,"vy_%s",plist->species[j_type]);
            vz_particle[j_type]=(GBPREAL *)ADaPS_fetch(plist->data,"vz_%s",plist->species[j_type]);
            id_particle[j_type]=(size_t  *)ADaPS_fetch(plist->data,"id_%s",plist->species[j_type]);
            if(fp_gadget->header.mass_array[i_type]>0.)
               ADaPS_store(&(plist->data),(void *)(&(fp_gadget->header.mass_array[i_type])),"mass_array_%s",ADaPS_SCALAR_DOUBLE,plist->species[j_type]);
         }
      }
   }
   x_particle [i_type][i_particle_type]=pos[0];
   y_particle [i_type][i_particle_type]=pos[1];
   z_particle [i_type][i_particle_type]=pos[2];
   vx_particle[i_type][i_particle_type]=vel[0];
   vy_particle[i_type][i_particle_type]=vel[1];
   vz_particle[i_type][i_particle_type]=vel[2];
   id_particle[i_type][i_particle_type]=ID_i;
}

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
                         int     mode);
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
   GBPREAL   *pos_buffer=(GBPREAL *)SID_malloc(3*sizeof(GBPREAL)*READ_BUFFER_SIZE_LOCAL);
   GBPREAL   *vel_buffer=(GBPREAL *)SID_malloc(3*sizeof(GBPREAL)*READ_BUFFER_SIZE_LOCAL);
   char      *ids_buffer=(char    *)SID_malloc(  sizeof(size_t) *READ_BUFFER_SIZE_LOCAL);
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
      fread(&record_length_open,4,1,fp_pos);
      fread(&header,sizeof(gadget_header_info),1,fp_pos);
      fread(&record_length_close,4,1,fp_pos);
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
      fread(&record_length_open,4,1,fp_vel);
      fseeko(fp_vel,(off_t)(record_length),SEEK_CUR);
      fread(&record_length_close,4,1,fp_vel);
      //if(record_length_open!=record_length_close)
      //  SID_log_warning("Problem with GADGET record size (close of positions)",ERROR_LOGIC);
      fseeko(fp_vel,(off_t)(4),SEEK_CUR);
      fseeko(fp_ids,(off_t)(2*4+record_length),SEEK_CUR);

      // Initialize IDs pointer
      record_length=3*n_particles_file*sizeof(GBPREAL);
      fread(&record_length_open,4,1,fp_ids);
      fseeko(fp_ids,(off_t)(record_length),SEEK_CUR);
      fread(&record_length_close,4,1,fp_ids);
      //if(record_length_open!=record_length_close)
      //  SID_log_warning("Problem with GADGET record size (close of velocities)",ERROR_LOGIC);
      fread(&record_length_open,4,1,fp_ids);

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
            n_buffer=MIN(READ_BUFFER_SIZE_LOCAL,header.n_file[i_type]-j_particle);
            if(flag_I_am_reader){
               fread(pos_buffer,sizeof(GBPREAL),3*n_buffer,fp_pos);
               fread(vel_buffer,sizeof(GBPREAL),3*n_buffer,fp_vel);
               fread(ids_buffer,ID_byte_size,     n_buffer,fp_ids);
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

int main(int argc, char *argv[]){
  SID_init(&argc,&argv,NULL,NULL);

  // Parse command line
  select_gadget_volume_params_info select_gadget_volume_params;
  int  snapshot;
  int  n_files_out;
  int  select_mode;
  char filename_in_root[MAX_FILENAME_LENGTH];
  char filename_out_root[MAX_FILENAME_LENGTH];
  GBPREAL cen_select[3];
  GBPREAL select_size;
  strcpy(filename_in_root, argv[1]);
  strcpy(filename_out_root,argv[7]);
  snapshot                          =         atoi(argv[2]);
  select_gadget_volume_params.cen[0]=(GBPREAL)atof(argv[3]);
  select_gadget_volume_params.cen[1]=(GBPREAL)atof(argv[4]);
  select_gadget_volume_params.cen[2]=(GBPREAL)atof(argv[5]);
  select_gadget_volume_params.size  =(GBPREAL)atof(argv[6]);
  n_files_out                       =         atoi(argv[8]);
  select_mode                       =         atoi(argv[9]);

  // Check that the selection mode is valid and set function pointer
  int (*select_function)(gadget_read_info *fp_gadget,
                         void             *params,
                         size_t            i_particle,
                         size_t            i_particle_type,
                         int               i_type,
                         GBPREAL          *pos,
                         GBPREAL          *vel,
                         size_t            ID_i);
  if(select_mode==1)
     select_function=select_gadget_cube;
  else if(select_mode==2){
     select_function=select_gadget_sphere;
     select_gadget_volume_params.size2=pow(select_gadget_volume_params.size,2.);
  }
  else
     SID_trap_error("Invalid selection mode (%d) given.",ERROR_SYNTAX,select_mode);

  SID_log("Excising volume from Gadget binary file {%s;snapshot=%d}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in_root,snapshot);

  // Initialize the plist data structure 
  plist_info plist;
  select_gadget_volume_params.plist=&plist;
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);

  // Read the header and determine the input file-format
  gadget_read_info   fp_gadget;
  int                flag_filefound=init_gadget_read(filename_in_root,snapshot,&fp_gadget);
  int                flag_multifile=fp_gadget.flag_multifile;
  int                flag_file_type=fp_gadget.flag_file_type;
  gadget_header_info header        =fp_gadget.header;
  if(!flag_filefound)
     SID_trap_error("File not found.",ERROR_LOGIC);
  select_gadget_volume_params.box_size=fp_gadget.header.box_size;

  // Count the particles
  size_t n_particles_type_local[N_GADGET_TYPE];
  size_t n_particles_type[N_GADGET_TYPE];
  int    flag_long_IDs;
  process_gadget_file("Counting particles in selection...",
                      filename_in_root,
                      snapshot,
                      select_function,
                      count_gadget_particles,
                      &select_gadget_volume_params,
                      n_particles_type_local,
                      n_particles_type,
                      &flag_long_IDs,
                      PROCESS_GADGET_BINARY_DEFAULT);

  // Allocate RAM for the particles
  allocate_gadget_particles(&plist,n_particles_type_local,n_particles_type,flag_long_IDs);

  // Read the particles
  process_gadget_file("Performing read/select/write...",
                      filename_in_root,
                      snapshot,
                      select_function,
                      store_gadget_particles,
                      &select_gadget_volume_params,
                      NULL,
                      NULL,
                      &flag_long_IDs,
                      PROCESS_GADGET_BINARY_DEFAULT);

  // Write the snapshot
  char filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_out,"%s_%03d",filename_out_root,snapshot);
  write_gadget_binary_new(&plist,filename_out,n_files_out,WRITE_GADGET_BINARY_DEFAULT);

  // Clean-up 
  free_plist(&plist);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

