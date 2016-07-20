#include <string.h>
#include <gbpLib.h>
#include <gbpHalos.h>

#define _FILE_OFFSET_BITS 64

int fread_catalog_file(fp_catalog_info *fp_in,halo_properties_SHORT_info *properties_short_out,halo_properties_SAGE_info *properties_out,halo_properties_info *properties_all_out,halo_profile_info *profiles_out,int halo_index){
  int n_skip;
  int r_val=0;

  if(halo_index<0||halo_index>=fp_in->n_halos_total)
     SID_trap_error("halo_index (%d) is out of range (0->%d) in fread_catalog_file().",ERROR_LOGIC,halo_index,fp_in->n_halos_total-1);

  // Skip to the right place (if need-be)
  if(halo_index!=fp_in->i_halo || halo_index>fp_in->i_halo_stop || halo_index<fp_in->i_halo_start){
     // We always have to scan forward, so if we're going backwards, we have to start from scratch
     if(halo_index<fp_in->i_halo)         fopen_nth_catalog_file(fp_in,0);
     while(halo_index>fp_in->i_halo_stop) fopen_nth_catalog_file(fp_in,fp_in->i_file+1);
     n_skip=halo_index-fp_in->i_halo;
     if(n_skip>0){
        if(fp_in->flag_read_properties)
           fseeko(fp_in->fp_properties,(off_t)(sizeof(halo_properties_info)*n_skip),SEEK_CUR);
        if(fp_in->flag_read_profiles){
           int i_profile;
           int n_bins;
           for(i_profile=0;i_profile<n_skip;i_profile++){
              fread_verify(&n_bins,sizeof(int),1,fp_in->fp_profiles);
              fseeko(fp_in->fp_profiles,(off_t)(n_bins*sizeof(halo_profile_bin_info)),SEEK_CUR);
           }
        }
     }
     else if(n_skip<0)
        SID_trap_error("Negative skips (%d) not supported in fread_catalog_file().",ERROR_LOGIC,n_skip);
     fp_in->i_halo+=n_skip;
  }

  // We must insist that something be read else the i_halo pointer will not work
  if(!(fp_in->flag_read_properties) && !(fp_in->flag_read_profiles))
     SID_trap_error("Nothing is being read in fread_catalog_file().",ERROR_LOGIC);

  // Read properties
  if(fp_in->flag_read_properties){
     // Perform Read
     halo_properties_info properties_in;

     fread_verify(&properties_in,sizeof(halo_properties_info),1,fp_in->fp_properties);
     if(properties_all_out!=NULL)
        memcpy(properties_all_out,&properties_in,sizeof(halo_properties_info));

     // Store the quantities we need to keep if alternate pointers are given
     if(properties_short_out!=NULL){
        properties_short_out->n_particles     =(int)   properties_in.n_particles;
        properties_short_out->M_vir           =(double)properties_in.M_vir;
        properties_short_out->R_vir           =(float) properties_in.R_vir;
        properties_short_out->pos[0]          =(float) properties_in.position_MBP[0];
        properties_short_out->pos[1]          =(float) properties_in.position_MBP[1];
        properties_short_out->pos[2]          =(float) properties_in.position_MBP[2];
        properties_short_out->vel[0]          =(float) properties_in.velocity_COM[0];
        properties_short_out->vel[1]          =(float) properties_in.velocity_COM[1];
        properties_short_out->vel[2]          =(float) properties_in.velocity_COM[2];
        properties_short_out->sigma_v         =(float) properties_in.sigma_v;
        properties_short_out->V_max           =(float) properties_in.V_max;
     }
     if(properties_out!=NULL){
        properties_out->descendant      =0;
        properties_out->progenitor_first=0;
        properties_out->progenitor_next =0;
        properties_out->group_halo_first=0;
        properties_out->group_halo_next =0;             
        properties_out->n_particles     =properties_in.n_particles;
        properties_out->M_Mean200       =(float)(properties_in.M_vir/1e10);
        properties_out->M_vir           =(float)(properties_in.M_vir/1e10);
        properties_out->M_TopHat        =(float)(properties_in.M_vir/1e10);
        properties_out->pos[0]          =properties_in.position_MBP[0];
        properties_out->pos[1]          =properties_in.position_MBP[1];
        properties_out->pos[2]          =properties_in.position_MBP[2];
        properties_out->vel[0]          =properties_in.velocity_COM[0];
        properties_out->vel[1]          =properties_in.velocity_COM[1];
        properties_out->vel[2]          =properties_in.velocity_COM[2];
        properties_out->sigma_v         =properties_in.sigma_v;
        properties_out->v_max           =properties_in.V_max;
        properties_out->spin[0]         =properties_in.spin[0];
        properties_out->spin[1]         =properties_in.spin[1];
        properties_out->spin[2]         =properties_in.spin[2];
        properties_out->most_bound_id   =properties_in.id_MBP;
        properties_out->snap_num        =fp_in->snap_num;
        properties_out->file_num        =0;
        properties_out->halo_index      =halo_index;
        properties_out->half_mass       =0.;
     }
  }
  else if(fp_in->flag_read_properties)
     SID_trap_error("File pointer not initialized while reading halo properties.",ERROR_LOGIC);

  // Read profiles
  if(fp_in->flag_read_profiles){
     fread_verify(&(profiles_out->n_bins),sizeof(int),                  1,                   fp_in->fp_profiles);
     fread_verify(&(profiles_out->bins),  sizeof(halo_profile_bin_info),profiles_out->n_bins,fp_in->fp_profiles);
  }
  else if(fp_in->flag_read_profiles)
     SID_trap_error("File pointer not initialized while reading halo profiles.",ERROR_LOGIC);

  // Set counter
  fp_in->i_halo++;

  return(r_val);
}

int fread_catalog_raw(fp_catalog_info *fp_in,halo_properties_info *properties_out,halo_profile_info *profiles_out,int halo_index){
  int n_skip;
  int r_val=0;

  if(halo_index<0||halo_index>=fp_in->n_halos_total)
     SID_trap_error("halo_index (%d) is out of range (0->%d) in fread_catalog_file().",ERROR_LOGIC,halo_index,fp_in->n_halos_total-1);

  // Skip to the right place (if need-be)
  if(halo_index!=fp_in->i_halo || halo_index>fp_in->i_halo_stop || halo_index<fp_in->i_halo_start){
     // We always have to scan forward, so if we're going backwards, we have to start from scratch
     if(halo_index<fp_in->i_halo_start)   fopen_nth_catalog_file(fp_in,0);
     while(halo_index>fp_in->i_halo_stop) fopen_nth_catalog_file(fp_in,fp_in->i_file+1);
     n_skip=halo_index-fp_in->i_halo;
     if(n_skip>0){
        if(fp_in->flag_read_properties)
           fseeko(fp_in->fp_properties,(off_t)(sizeof(halo_properties_info)*n_skip),SEEK_CUR);
        if(fp_in->flag_read_profiles){
           int i_profile;
           int n_bins;
           for(i_profile=0;i_profile<n_skip;i_profile++){
              fread_verify(&n_bins,sizeof(int),1,fp_in->fp_profiles);
              fseeko(fp_in->fp_profiles,(off_t)(n_bins*sizeof(halo_profile_bin_info)),SEEK_CUR);
           }
        }
     }
     fp_in->i_halo+=n_skip;
  }

  // We must insist that something be read else the i_halo pointer will not work
  if(!(fp_in->flag_read_properties) && !(fp_in->flag_read_profiles))
     SID_trap_error("Nothing is being read in fread_catalog_file().",ERROR_LOGIC);

  // Read properties
  if(fp_in->flag_read_properties){
     // Perform Read
     fread_verify(properties_out,sizeof(halo_properties_info),1,fp_in->fp_properties);
  }
  else if(fp_in->flag_read_properties)
     SID_trap_error("File pointer not initialized while reading halo properties.",ERROR_LOGIC);

  // Read profiles
  if(fp_in->flag_read_profiles){
     fread_verify(&(profiles_out->n_bins),sizeof(int),                  1,                   fp_in->fp_profiles);
     fread_verify(&(profiles_out->bins),  sizeof(halo_profile_bin_info),profiles_out->n_bins,fp_in->fp_profiles);
  }
  else if(fp_in->flag_read_profiles)
     SID_trap_error("File pointer not initialized while reading halo profiles.",ERROR_LOGIC);

  // Set counter
  fp_in->i_halo++;

  return(r_val);
}

