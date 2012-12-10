#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHalos.h>

void read_group_properties(FILE *fp,halo_info *properties_out,int halo_index,int i_read){
  halo_properties_info properties_in;

  fread(&properties_in,sizeof(properties_in),1,fp);

  properties_out->descendant      =0;
  properties_out->progenitor_first=0;
  properties_out->progenitor_next =0;
  properties_out->group_halo_first=0;
  properties_out->group_halo_next =0;             
  properties_out->match_type      =0;             
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
  properties_out->snap_num        =i_read;
  properties_out->halo_index      =halo_index;
  properties_out->halo_id         =0;
  properties_out->group_id        =0;
                
}
