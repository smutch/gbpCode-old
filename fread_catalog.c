#include <string.h>
#include <gbpLib.h>
#include <gbpHalos.h>

int fread_catalog_file(fp_catalog_info *fp_in,halo_info *properties_out,halo_profile_info *profiles_out,int halo_index){
  int n_skip;
  // Skip to the right place
  if(halo_index!=(fp_in->i_halo+1) || halo_index>fp_in->i_halo_stop){
     // We always have to scan forward, so if we're going backwards, we have to start from scratch
     if(halo_index<fp_in->i_halo_start)   fopen_nth_catalog_file(fp_in,0);
     while(halo_index>fp_in->i_halo_stop) fopen_nth_catalog_file(fp_in,fp_in->i_file+1);
     n_skip=halo_index-fp_in->i_halo_start;
     if(n_skip>0){
        if(fp_in->flag_read_properties)
           fseeko64(fp_in->fp_properties,(off64_t)(sizeof(halo_properties_info)*n_skip),SEEK_CUR);
        if(fp_in->flag_read_profiles){
           int    i_profile;
           int    n_bins;
           size_t n_bytes_skip;
           if(SID.I_am_Master){
              for(i_profile=0,n_bytes_skip=0;i_profile<n_skip;i_profile++){
                 fread(&n_bins,sizeof(int),1,fp_in->fp_profiles);
                 fseeko64(fp_in->fp_profiles,(off64_t)(n_bins*sizeof(halo_profile_bin_info)),SEEK_CUR);
                 n_bytes_skip+=sizeof(int)+(size_t)n_bins*sizeof(halo_profile_bin_info);
              }
           }
           SID_Bcast(&n_bytes_skip,sizeof(size_t),MASTER_RANK,SID.COMM_WORLD);
           if(!SID.I_am_Master)
              fseeko64(fp_in->fp_profiles,(off64_t)(n_bytes_skip),SEEK_CUR);
        }
     }
  }

  // Read properties
  if(fp_in->flag_read_properties){
     // Perform Read
     halo_properties_info properties_in;
     fread(&properties_in,sizeof(properties_in),1,fp_in->fp_properties);

     // Store the quantities we need to keep
     properties_out->descendant      =0;
     properties_out->progenitor_first=0;
     properties_out->progenitor_next =0;
     properties_out->group_halo_first=0;
     properties_out->group_halo_next =0;             
     properties_out->n_particles     =properties_in.n_particles;
     properties_out->M_vir           =(float)(properties_in.M_vir/1e10);
     properties_out->R_vir           =properties_in.R_vir;
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
     properties_out->halo_index      =halo_index;
     properties_out->halo_id         =0;
     properties_out->group_id        =0;
  }

  // Read profiles
  if(fp_in->flag_read_profiles){
     fread(&(profiles_out->n_bins),sizeof(int),                  1,                   fp_in->fp_profiles);
     fread(&(profiles_out->bins),  sizeof(halo_profile_bin_info),profiles_out->n_bins,fp_in->fp_profiles);
  }

  // Set counter
  fp_in->i_halo=halo_index;
}

