#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpClustering.h>
#include <gbpHalos.h>
#include <gbpSSimPL.h>

void swap_endian_gadget_header_local(gadget_header_info *header);
void swap_endian_gadget_header_local(gadget_header_info *header){
   swap_endian((char *)(&(header->n_file)),         N_GADGET_TYPE,sizeof(int));
   swap_endian((char *)(&(header->mass_array)),     N_GADGET_TYPE,sizeof(double));
   swap_endian((char *)(&(header->time)),           1,            sizeof(double));
   swap_endian((char *)(&(header->redshift)),       1,            sizeof(double));
   swap_endian((char *)(&(header->flag_SFr)),       1,            sizeof(int));
   swap_endian((char *)(&(header->flag_feedback)),  1,            sizeof(int));
   swap_endian((char *)(&(header->n_all_lo_word)),  N_GADGET_TYPE,sizeof(unsigned int));
   swap_endian((char *)(&(header->flag_cooling)),   1,            sizeof(int));
   swap_endian((char *)(&(header->n_files)),        1,            sizeof(int));
   swap_endian((char *)(&(header->box_size)),       1,            sizeof(double));
   swap_endian((char *)(&(header->Omega_M)),        1,            sizeof(double));
   swap_endian((char *)(&(header->Omega_Lambda)),   1,            sizeof(double));
   swap_endian((char *)(&(header->h_Hubble)),       1,            sizeof(double));
   swap_endian((char *)(&(header->flag_ages)),      1,            sizeof(int));
   swap_endian((char *)(&(header->flag_metals)),    1,            sizeof(int));
   swap_endian((char *)(&(header->n_all_hi_word)),  N_GADGET_TYPE,sizeof(unsigned int));
   swap_endian((char *)(&(header->flag_entropyICs)),1,            sizeof(int));
}

int swap_endian_snapshot(const char *filename_in_root,const char *filename_out_root,int region_number,int snap_number,int mode,int *IDs_byte_size_return){

  if(region_number<0)
     SID_log("Swapping endian of full snapshot...",SID_LOG_OPEN);
  else
     SID_log("Swapping endian of region #%03d snapshot...",SID_LOG_OPEN,region_number);

  // Sanity check
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE) && check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE))
     SID_trap_error("Invalid mode flag (%d) in swap_endian_catalogs_properties_local().",ERROR_LOGIC,mode);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  char filename_out[MAX_FILENAME_LENGTH];
  if(region_number<0){
     sprintf(filename_in, "%s/snapshots/snapshot_%03d/snapshot_%03d.%d",filename_in_root, snap_number,snap_number,0);
     sprintf(filename_out,"%s/snapshots/snapshot_%03d/snapshot_%03d.%d",filename_out_root,snap_number,snap_number,0);
  }
  else{
     sprintf(filename_in, "%s/snapshots/snapshot_region%03d_%03d/snapshot_region%03d_%03d.%d",filename_in_root, region_number,snap_number,region_number,snap_number,0);
     sprintf(filename_out,"%s/snapshots/snapshot_region%03d_%03d/snapshot_region%03d_%03d.%d",filename_out_root,region_number,snap_number,region_number,snap_number,0);
  }

  // Open input and output files
  FILE *fp_in =NULL;
  FILE *fp_out=NULL;
  int   flag_type=FALSE;
  if((fp_in=fopen(filename_in,"r"))==NULL){
     if(region_number<0)
        sprintf(filename_in, "%s/snapshots/snapshot_%03d",filename_in_root,snap_number);
     else
        sprintf(filename_in, "%s/snapshots/snapshot_region%03d_%03d",filename_in_root,region_number,snap_number);
     if((fp_in=fopen(filename_in,"r"))==NULL){
        SID_log("not present.",SID_LOG_CLOSE);
        return(FALSE);
     }
     flag_type=TRUE;
  }
  else{
     char filename_dir[MAX_FILENAME_LENGTH];
     if(region_number<0)
        sprintf(filename_dir,"%s/snapshots/snapshot_%03d",filename_out_root,snap_number);
     else
        sprintf(filename_dir,"%s/snapshots/snapshot_region%03d_%03d",filename_out_root,region_number,snap_number);
     mkdir(filename_dir,02755);
  }

  // Read the needed header information
  gadget_header_info header;
  int block_size_in;
  int block_size_out;
  int i_file_in;
  int n_files;
  fread(&block_size_in, sizeof(int),               1,fp_in);
  fread(&header,        sizeof(gadget_header_info),1,fp_in);
  fread(&block_size_out,sizeof(int),               1,fp_in);
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE))
     swap_endian_gadget_header_local(&header);
  n_files=header.n_files;
  fclose(fp_in);

  // Figure-out what sort of IDs we're working with (loop over files
  //    because the first few may have no particles in them to figure
  //    this out with)
  int IDs_byte_size   =0;
  int n_particles_file=0;
  int flag_ID_size_not_set=TRUE;
  for(int i_file=0;i_file<n_files && flag_ID_size_not_set;i_file++){
     if(flag_type){
        if(region_number<0){
           sprintf(filename_in, "%s/snapshots/snapshot_%03d",filename_in_root, snap_number);
           sprintf(filename_out,"%s/snapshots/snapshot_%03d",filename_out_root,snap_number);
        }
        else{
           sprintf(filename_in, "%s/snapshots/snapshot_region%03d_%03d",filename_in_root, region_number,snap_number);
           sprintf(filename_out,"%s/snapshots/snapshot_region%03d_%03d",filename_out_root,region_number,snap_number);
        }
     }
     else{
        if(region_number<0){
           sprintf(filename_in, "%s/snapshots/snapshot_%03d/snapshot_%03d.%d",filename_in_root, snap_number,snap_number,i_file);
           sprintf(filename_out,"%s/snapshots/snapshot_%03d/snapshot_%03d.%d",filename_out_root,snap_number,snap_number,i_file);
        }
        else{
           sprintf(filename_in, "%s/snapshots/snapshot_region%03d_%03d/snapshot_region%03d_%03d.%d",filename_in_root,
                   region_number,snap_number,region_number,snap_number,i_file);
           sprintf(filename_out,"%s/snapshots/snapshot_region%03d_%03d/snapshot_region%03d_%03d.%d",filename_out_root,
                   region_number,snap_number,region_number,snap_number,i_file);
        }
     }
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);
     fread(&block_size_in, sizeof(int),               1,fp_in);
     fread(&header,        sizeof(gadget_header_info),1,fp_in);
     fread(&block_size_out,sizeof(int),               1,fp_in);
     if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE))
        swap_endian_gadget_header_local(&header);
     for(int i_type=0;i_type<N_GADGET_TYPE;i_type++)
        n_particles_file+=header.n_file[i_type];
     if(n_particles_file>0){
        fseeko(fp_in,(off_t)(6*sizeof(int)+sizeof(gadget_header_info)+n_particles_file*6*sizeof(GBPREAL)),SEEK_SET);
        fread(&block_size_in, sizeof(int),1,fp_in);
        if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE))
           swap_endian((char *)(&block_size_in),1,sizeof(int));
        if(n_particles_file==(block_size_in/sizeof(int)))
           IDs_byte_size=sizeof(int);
        else if(n_particles_file==(block_size_in/sizeof(long long)))
           IDs_byte_size=sizeof(long long);
        else
           SID_trap_error("Can not determine particle ID size from IDs block size (%d).",ERROR_LOGIC,block_size_in);
        SID_log("%d-byte IDs...",SID_LOG_CONTINUE,IDs_byte_size);
        flag_ID_size_not_set=FALSE;
     }
     fclose(fp_in);
  }

  // Sanity check
  if(flag_type && n_files!=1) SID_trap_error("Single file format dataset {%s} has n_files=%d (ie !=1) in its header.",ERROR_LOGIC,filename_in,n_files);

  int i_file;
  for(i_file=0;i_file<n_files;i_file++){
     // Open the files
     if(flag_type){
        if(region_number<0){
           sprintf(filename_in, "%s/snapshots/snapshot_%03d",filename_in_root, snap_number);
           sprintf(filename_out,"%s/snapshots/snapshot_%03d",filename_out_root,snap_number);
        }
        else{
           sprintf(filename_in, "%s/snapshots/snapshot_region%03d_%03d",filename_in_root, region_number,snap_number);
           sprintf(filename_out,"%s/snapshots/snapshot_region%03d_%03d",filename_out_root,region_number,snap_number);
        }
     }
     else{
        if(region_number<0){
           sprintf(filename_in, "%s/snapshots/snapshot_%03d/snapshot_%03d.%d",filename_in_root, snap_number,snap_number,i_file);
           sprintf(filename_out,"%s/snapshots/snapshot_%03d/snapshot_%03d.%d",filename_out_root,snap_number,snap_number,i_file);
        }
        else{
           sprintf(filename_in, "%s/snapshots/snapshot_region%03d_%03d/snapshot_region%03d_%03d.%d",filename_in_root, 
                   region_number,snap_number,region_number,snap_number,i_file);
           sprintf(filename_out,"%s/snapshots/snapshot_region%03d_%03d/snapshot_region%03d_%03d.%d",filename_out_root,
                   region_number,snap_number,region_number,snap_number,i_file);
        }
     }
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);
     if((fp_out=fopen(filename_out,"w"))==NULL)
        SID_trap_error("Could not open {%s} for writing.",ERROR_IO_OPEN,filename_out);

     // Allocate buffer
     char *buffer=(char *)SID_malloc(sizeof(double)*3);

     // Process the header
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),buffer);
     fread(&header,sizeof(gadget_header_info),1,fp_in);
     if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE))
        swap_endian_gadget_header_local(&header);
     fwrite(&header,sizeof(gadget_header_info),1,fp_out);
     n_particles_file=0;
     for(int i_type=0;i_type<N_GADGET_TYPE;i_type++)
        n_particles_file+=header.n_file[i_type];
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),buffer);

     // Process positions
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),buffer);
     for(int i_halo=0;i_halo<n_particles_file;i_halo++)
        rewrite_swap_endian(fp_in,fp_out,3,sizeof(GBPREAL),buffer);
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),buffer);

     // Process velocities
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),buffer);
     for(int i_halo=0;i_halo<n_particles_file;i_halo++)
        rewrite_swap_endian(fp_in,fp_out,3,sizeof(GBPREAL),buffer);
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),buffer);

     // Process IDs
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),buffer);
     for(int i_halo=0;i_halo<n_particles_file;i_halo++)
        rewrite_swap_endian(fp_in,fp_out,1,IDs_byte_size,buffer);
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),buffer);

     // Free the buffer
     SID_free(SID_FARG buffer);

     // Close files
     fclose(fp_in);
     fclose(fp_out);
  }
 
  (*IDs_byte_size_return)=IDs_byte_size; 
  SID_log("Done.",SID_LOG_CLOSE);
  return(TRUE);
}

