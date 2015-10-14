#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpSPH.h>

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
  strcpy(filename_in_root,                         argv[1]);
  snapshot                          =         atoi(argv[2]);
  select_gadget_volume_params.cen[0]=(GBPREAL)atof(argv[3]);
  select_gadget_volume_params.cen[1]=(GBPREAL)atof(argv[4]);
  select_gadget_volume_params.cen[2]=(GBPREAL)atof(argv[5]);
  select_gadget_volume_params.size  =(GBPREAL)atof(argv[6]);
  strcpy(filename_out_root,                        argv[7]);
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
                      process_gadget_file_fctn_null,
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

