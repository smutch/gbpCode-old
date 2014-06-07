#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpSPH.h>

typedef select_gadget_volume_params_info struct select_gadget_volume_params_info;
struct select_gadget_volume_params_info{
   plist_info *plist;
   GBPREAL     cen[3];
   GBPREAL     size;
};

void allocate_gadget_particles(plist_info *plist, 
                               size_t     *n_particles_type_local,
                               size_t     *n_particles_type);
void allocate_gadget_particles(plist_info *plist, 
                               size_t     *n_particles_type_local,
                               size_t     *n_particles_type){
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
      }
      ADaPS_store(&(plist->data),(void *)(&(n_particles_type_local[i_type])),"n_%s",    ADaPS_SCALAR_SIZE_T,plist->species[i_type]);
      ADaPS_store(&(plist->data),(void *)(&(n_particles_type[i_type])),      "n_all_%s",ADaPS_SCALAR_SIZE_T,plist->species[i_type]);
   }
}

int select_gadget_cube(fp_gadget_info *fp_gadget,
                       void           *params,
                       size_t          i_particle,
                       size_t          i_particle_type,
                       int             i_type,
                       GBPREAL        *pos,
                       GBPREAL        *vel,
                       size_t          ID_i);
int select_gadget_cube(fp_gadget_info *fp_gadget,
                       void           *params,
                       size_t          i_particle,
                       size_t          i_particle_type,
                       int             i_type,
                       GBPREAL        *pos,
                       GBPREAL        *vel,
                       size_t          ID_i){
}

int store_gadget_particle(fp_gadget_info *fp_gadget,
                          void           *params,
                          size_t          i_particle,
                          size_t          i_particle_type,
                          int             i_type,
                          GBPREAL        *pos,
                          GBPREAL        *vel,
                          size_t          ID_i){
   static GBPREAL *x_particle[N_GADGET_TYPE];
   static GBPREAL *y_particle[N_GADGET_TYPE];
   static GBPREAL *z_particle[N_GADGET_TYPE];
   static GBPREAL *vx_particle[N_GADGET_TYPE];
   static GBPREAL *vy_particle[N_GADGET_TYPE];
   static GBPREAL *vz_particle[N_GADGET_TYPE];
   static size_t  *id_particle[N_GADGET_TYPE];
   if(fp_gadget->flag_first_select){
   }
   x_particle[i_type][i_particle_type] =pos[0];
   y_particle[i_type][i_particle_type] =pos[1];
   z_particle[i_type][i_particle_type] =pos[2];
   vx_particle[i_type][i_particle_type]=vel[0];
   vy_particle[i_type][i_particle_type]=vel[1];
   vz_particle[i_type][i_particle_type]=vel[2];
   id_particle[i_type][i_particle_type]=ID_i;
}



void process_gadget_file(char   *filename_read_root,
                         int     snapshot_number,
                         int     select_function(),
                         void    action_function(),
                         void   *params,
                         size_t *n_particles_type_local_pass,
                         size_t *n_particles_type_pass,
                         int     mode){
  SID_log("Processing Gadget binary file...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Allocate counter arrays if we need to
  int flag_n_particles_type_local_allocate=FALSE;
  if(n_particles_type_local_pass==NULL){
     n_particles_type_local=(size_t *)SID_malloc(sizeof(size_t)*N_GADGET_TYPE);
     flag_n_particles_type_local_allocate=TRUE;
  }
  else
     n_particles_type_local=n_particles_type_local_pass;
  int flag_n_particles_type_local_allocate=FALSE;
  if(n_particles_type_local_pass==NULL){
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
  int                flag_filefound=init_gadget_read(filename_read_root,snapshot,&fp_gadget);
  int                flag_multifile=fp_gadget.flag_multifile;
  int                flag_file_type=fp_gadget.flag_file_type;
  gadget_header_info header        =fp_gadget.header;
  if(!flag_filefound)
     SID_trap_error("File not found.",ERROR_LOGIC);

  // Decide which files each rank will process
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

  // Loop over the files
  pcounter_info pcounter;
  int           flag_active=TRUE;
  SID_init_pcounter(&pcounter,n_particles_snap,10);
  for(;i_file<fp_gadget.header.n_files;i_file+=i_file_skip){

     // Open file pointers
     set_gadget_filename(fp_gadget,i_file,filename);
     FILE *fp_pos=fopen(filename,"r");
     FILE *fp_vel=fopen(filename,"r");
     FILE *fp_IDs=fopen(filename,"r");

     // Header
     fread(&record_length_open,4,1,fp_pos);
     fread(&header,sizeof(gadget_header_info),1,fp_pos);
     fread(&record_length_close,4,1,fp_pos);
     if(record_length_open!=record_length_close)
       SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
     if(record_length_open!=sizeof(gadget_header_info))
       SID_log_warning("Problem with GADGET header size",ERROR_LOGIC);
     fread(&record_length_open,4,1,fp_pos);
     fseeko(fp_vel,(off_t)(2*4+record_length_open),SEEK_CUR);
     fseeko(fp_IDs,(off_t)(2*4+record_length_open),SEEK_CUR);

     // Positions

     // Velocities

     // IDs
     int    flag_long_ids;
     size_t ID_byte_size;

     // Process the file in buffered chunks
     int    i_type;
     size_t i_particle;
     size_t j_particle;
     for(i_type=0,i_particle=0;i_type<N_GADGET_TYPE;i_type++){
        for(j_particle=0;j_particle<header.n_file[i_type];i_particle+=n_buffer){
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
              if(select_function()){
                 action_function();
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

  // Clean-up
  if(flag_n_particles_type_local_allocate)
     SID_free(SID_FARG n_particles_type_local);
  if(flag_n_particles_type_allocate)
     SID_free(SID_FARG n_particles_type);

  SID_log("Done.",SID_LOG_CLOSE);
}

int main(int argc, char *argv[]){
  SID_init(&argc,&argv,NULL);

  // Parse command line
  select_gadget_volume_params_info select_gadget_volume_params;
  int  snapshot;
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

  SID_log("Excising cube from Gadget bindary file {%s;snapshot=%d}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in_root,snapshot);

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
  process_gadget_file(filename_in_root,
                      snapshot,
                      select_gadget_cube,
                      count_gadget_particles,
                      &select_gadget_volume_params,
                      n_particles_type_local,
                      n_particles_type,
                      PROCESS_GADGET_BINARY_ALL_TO_ALL);

  // Allocate RAM for the particles
  allocate_gadget_particles(&plist,n_particles_type_local,n_particles_type);

  // Read the particles
  process_gadget_file(filename_in_root,
                      snapshot,
                      select_gadget_cube,
                      store_gadget_particles,
                      &select_gadget_volume_params,
                      NULL,
                      NULL,
                      PROCESS_GADGET_BINARY_ALL_TO_ALL);

  // Write the snapshot
  char filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_out,"%s_%03d",filename_out_root,snapshot);
  write_gadget_binary(filename_out,&plist);

  // Clean-up 
  free_plist(&plist);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

