#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpSPH.h>

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

void process_gadget_file_fctn_null(gadget_read_info *fp_gadget,
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

int select_gadget_ids(gadget_read_info *fp_gadget,
                      void             *params,
                      size_t            i_particle,
                      size_t            i_particle_type,
                      int               i_type,
                      GBPREAL          *pos,
                      GBPREAL          *vel,
                      size_t            ID_i){
   static int      flag_long_ids;
   static int      n_ids;
   static size_t  *id_list;
   if(fp_gadget->first_select_call){
      n_ids  =((select_gadget_ids_params_info *)params)->n_ids;
      id_list=((select_gadget_ids_params_info *)params)->id_list;
   }
   return(is_a_member(&ID_i,id_list,n_ids,SID_SIZE_T));
}

int select_gadget_all(gadget_read_info *fp_gadget,
                      void             *params,
                      size_t            i_particle,
                      size_t            i_particle_type,
                      int               i_type,
                      GBPREAL          *pos,
                      GBPREAL          *vel,
                      size_t            ID_i){
   return(TRUE);
}

