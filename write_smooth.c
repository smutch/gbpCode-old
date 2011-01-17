#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>

void write_smooth(plist_info *plist,
                  char       *species_name,
                  char       *filename,
                  int         n_chunks,
                  int         mode){
  void   *local_array;
  size_t  n_particles_local;
  size_t  n_particles_bcast;
  size_t  n_particles_total;
  size_t  i_particle_start_local;
  int     i_quantity;
  int     n_quantities;
  size_t  n_write;
  int     i_rank;
  double  array_mean;
  SID_fp                 fp_write_smoothfile;
  smoothfile_header_info write_smoothfile_header;

  SID_log("Writing {%s} particles to smooth file {%s}...",SID_LOG_OPEN,species_name,filename);

  // Count how many (and which) arrays are to be written
  write_smoothfile_header.n_quantities     =SMOOTH_N_QUANTITIES;
  write_smoothfile_header.n_quantities_used=0;
  for(i_quantity=0;i_quantity<write_smoothfile_header.n_quantities;i_quantity++){
    write_smoothfile_header.used[i_quantity]=FALSE;
    switch(i_quantity){
    case 0:
      if(ADaPS_exist(plist->data,"r_smooth_%s",species_name)){
        write_smoothfile_header.used[i_quantity]=TRUE;
      }
      break;
    case 1:
      if(ADaPS_exist(plist->data,"rho_%s",species_name)){
        write_smoothfile_header.used[i_quantity]=TRUE;
      }
      break;
    case 2:
      if(ADaPS_exist(plist->data,"sigma_v_%s",species_name)){
        write_smoothfile_header.used[i_quantity]=TRUE;
      }
      break;
    default:
      SID_trap_error("Unspecified case in write_smooth",ERROR_LOGIC);
      break;
    }
    if(write_smoothfile_header.used[i_quantity])
      write_smoothfile_header.n_quantities_used++;
  }

  // Fetch the number of particles
  n_particles_local=((size_t *)ADaPS_fetch(plist->data,"n_%s",    species_name))[0]; 
  n_particles_total=((size_t *)ADaPS_fetch(plist->data,"n_all_%s",species_name))[0]; 

  // Compute local file offsets
  i_particle_start_local=0;
#ifdef USE_MPI
  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
    if(SID.My_rank==i_rank)
      n_particles_bcast=n_particles_local;
    MPI_Bcast(&n_particles_bcast,1,MPI_SIZE_T,i_rank,MPI_COMM_WORLD);
    if(i_rank<SID.My_rank)
      i_particle_start_local+=n_particles_bcast;
  }
#endif

  // Open smooth file
  strcpy(write_smoothfile_header.species,species_name);
  n_write=write_smoothfile_header.n_quantities_used*n_particles_total;
  SID_fopen_chunked(filename,
                    "w",
                    &fp_write_smoothfile,
                    &write_smoothfile_header,
                    sizeof(smoothfile_header_info),
                    n_write,
                    sizeof(REAL),
                    n_chunks);

  // Write arrays (for backwards compatibility, add new quantities to the end of the list)
  for(i_quantity=0;i_quantity<write_smoothfile_header.n_quantities;i_quantity++){
    if(write_smoothfile_header.used[i_quantity]){
      switch(i_quantity){
      case 0:
        SID_log("Writing smoothing lengths...",SID_LOG_OPEN);
        local_array=ADaPS_fetch(plist->data,"r_smooth_%s",species_name);
        SID_fwrite_chunked(local_array,
                           n_particles_local,
                           i_particle_start_local,
                           &fp_write_smoothfile);
        calc_stat(local_array,NULL,n_particles_local,ADaPS_REAL,CALC_STAT_MEAN|CALC_STAT_RETURN_DOUBLE,&array_mean);
        SID_log("Done. (mean=%le kpc)",SID_LOG_CLOSE,array_mean/M_PER_KPC);
        break;
      case 1:
        SID_log("Writting densities...",SID_LOG_OPEN);
        local_array=ADaPS_fetch(plist->data,"rho_%s",species_name);
        SID_fwrite_chunked(local_array,
                           n_particles_local,
                           i_particle_start_local,
                           &fp_write_smoothfile);
        calc_stat(local_array,NULL,n_particles_local,ADaPS_REAL,CALC_STAT_MEAN|CALC_STAT_RETURN_DOUBLE,&array_mean);
        SID_log("Done. (mean=%le M_sol/Mpc^3)",SID_LOG_CLOSE,array_mean/(M_SOL/pow(M_PER_MPC,3.)));
        break;
      case 2:
        SID_log("Writing sigma_v's...",SID_LOG_OPEN);
        local_array=ADaPS_fetch(plist->data,"sigma_v_%s",species_name);
        SID_fwrite_chunked(local_array,
                           n_particles_local,
                           i_particle_start_local,
                           &fp_write_smoothfile);
        calc_stat(local_array,NULL,n_particles_local,ADaPS_REAL,CALC_STAT_MEAN|CALC_STAT_RETURN_DOUBLE,&array_mean);
        SID_log("Done. (mean=%le km/s)",SID_LOG_CLOSE,array_mean/1e3);
        break;
      }
    }
  }
  SID_fclose_chunked(&fp_write_smoothfile);
  SID_log("Done.",SID_LOG_CLOSE);
}
