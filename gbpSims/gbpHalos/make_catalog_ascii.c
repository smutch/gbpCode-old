#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHalos.h>

int main(int argc, char *argv[]){
  char    filename_in_root[256];
  char    filename_in_base[256];
  char    filename_properties[256];
  char    filename_profiles[256];
  char    filename_out_root[256];
  char    filename_out[256];
  char    prefix_text[5];
  FILE   *fp_out       =NULL;
  int     snap_number;
  int     snap_number_start;
  int     snap_number_stop;
  int     snap_number_step;
  float   lambda,v_c;
  float   offset_COM;
  float   r_min,r_max;
  double  box_size;

  SID_init(&argc,&argv,NULL);

  strcpy(filename_in_root,    argv[1]);
  box_size           =atof(argv[2]);
  snap_number_start  =atoi(argv[3]);
  snap_number_stop   =atoi(argv[4]);
  snap_number_step   =atoi(argv[5]);
  strcpy(filename_out_root,argv[6]);

  int flag_use_profiles=FALSE;

  sprintf(filename_in_base,"%s",filename_in_root);
  strip_path(filename_in_base);

  if(SID.I_am_Master){
    int i_type;
    SID_log("Processing catalogs for snaps %d->%d...",SID_LOG_OPEN|SID_LOG_TIMER,snap_number_start,snap_number_stop);
    SID_log("Properties structure size=%lld",SID_LOG_COMMENT,sizeof(halo_properties_info));
    for(snap_number=snap_number_start;snap_number<=snap_number_stop;snap_number++){
       for(i_type=0;i_type<2;i_type++){
         fp_catalog_info fp_catalog;
         int             catalog_mode;
         switch(i_type){
            case 0:
               sprintf(prefix_text,"");
               catalog_mode=READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES;
               break;
            case 1:
               sprintf(prefix_text,"sub");
               catalog_mode=READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES;
               break;
         }
         if(flag_use_profiles)
            catalog_mode|=READ_CATALOG_PROPERTIES;
         fopen_catalog(filename_in_root,
                       snap_number,
                       catalog_mode,
                       &fp_catalog);

         // Process halos
         sprintf(filename_out,"%s_%03d_%sgroups.ascii",filename_out_root,snap_number,prefix_text);
         SID_log("Writing snap #%03d %sgroup properties->to file {%s}...",SID_LOG_OPEN,snap_number,prefix_text,filename_out);
         SID_log("(%d %sgroups)...",SID_LOG_CONTINUE,fp_catalog.n_halos_total,prefix_text);
         fp_out=fopen(filename_out,"w");
         fprintf(fp_out,"# ASCII %sgroup catalog of snap #%d of %s\n",prefix_text,snap_number,filename_in_root);
         fprintf(fp_out,"# File columns: (1)     %sgroup number\n",prefix_text);
         fprintf(fp_out,"#               (2)     # of particles\n");
         fprintf(fp_out,"#               (3)     id_MBP\n");
         fprintf(fp_out,"#               (4)     x_COM      [Mpc/h] \n");
         fprintf(fp_out,"#               (5)     y_COM      [Mpc/h] \n");
         fprintf(fp_out,"#               (6)     z_COM      [Mpc/h] \n");
         fprintf(fp_out,"#               (7)     v_x_COM    [km/s]\n");
         fprintf(fp_out,"#               (8)     v_y_COM    [km/s]\n");
         fprintf(fp_out,"#               (9)     v_z_COM    [km/s]\n");
         fprintf(fp_out,"#               (10)    x_MBP      [Mpc/h] \n");
         fprintf(fp_out,"#               (11)    y_MBP      [Mpc/h] \n");
         fprintf(fp_out,"#               (12)    z_MBP      [Mpc/h] \n");
         fprintf(fp_out,"#               (13)    v_x_MBP    [km/s]\n");
         fprintf(fp_out,"#               (14)    v_y_MBP    [km/s]\n");
         fprintf(fp_out,"#               (15)    v_z_MBP    [km/s]\n");
         fprintf(fp_out,"#               (16)    M_vir      [M_sol/h]\n");
         fprintf(fp_out,"#               (17)    R_vir      [Mpc/h]\n");
         fprintf(fp_out,"#               (18)    R_halo     [Mpc/h]\n");
         fprintf(fp_out,"#               (19)    R_max      [Mpc/h]\n");
         fprintf(fp_out,"#               (20)    V_max      [km/s]\n");
         fprintf(fp_out,"#               (21)    V_vir      [km/s]\n");
         fprintf(fp_out,"#               (22)    sigma_v    [km/s]\n");
         fprintf(fp_out,"#               (23-25) spin       [Mpc/h km/s]\n");
         fprintf(fp_out,"#               (26)    lambda (spin parameter)\n");
         fprintf(fp_out,"#               (27)    q_triaxial\n");
         fprintf(fp_out,"#               (28)    s_triaxial\n");
         fprintf(fp_out,"#               (29)    offset_COM [kpc/h]\n");
         if(flag_use_profiles)
            fprintf(fp_out,"#               (30)    overdensity(R_halo)\n");

         halo_properties_info *properties=NULL;
         halo_profile_info    *profile   =NULL;
         properties = (halo_properties_info *)SID_malloc(sizeof(halo_properties_info));
         if(flag_use_profiles)
            profile = (halo_profile_info *)SID_malloc(sizeof(halo_profile_info));
         for(int i_group=0;i_group<fp_catalog.n_halos_total;i_group++){
           fread_catalog_file(&fp_catalog,NULL,properties,profile,i_group);

           float overdensity;
           if(flag_use_profiles){
              int n_bins=profile->n_bins;
              overdensity=profile->bins[n_bins-1].overdensity;
           }
           else
              overdensity=-1.;

           // Create a few properties
           double dx,dy,dz;
           v_c       =sqrt(G_NEWTON*properties->M_vir*M_SOL/(properties->R_vir*M_PER_MPC))*1e-3;
           lambda    =sqrt(properties->spin[0]*properties->spin[0]+properties->spin[1]*properties->spin[1]+properties->spin[2]*properties->spin[2])/(sqrt(2.)*properties->R_vir*v_c);
           dx        =d_periodic(properties->position_COM[0]-properties->position_MBP[0],box_size);
           dy        =d_periodic(properties->position_COM[1]-properties->position_MBP[1],box_size);
           dz        =d_periodic(properties->position_COM[2]-properties->position_MBP[2],box_size);
           offset_COM=sqrt(dx*dx+dy*dy+dz*dz);

           // Perform write
           fprintf(fp_out,"%9d %9d %9lld  %12.5e %12.5e %12.5e %11.5f %11.5f %11.5f  %12.5e %12.5e %12.5e %11.5f %11.5f %11.5f  %10.5le %11.5f %11.5f %11.5f  %11.5f %11.5f %11.5f  %12.5e %12.5e %12.5e %11.5f  %11.5f %11.5f  %11.5f",
                   i_group,properties->n_particles,properties->id_MBP,
                   properties->position_COM[0],properties->position_COM[1],properties->position_COM[2],
                   properties->velocity_COM[0],properties->velocity_COM[1],properties->velocity_COM[2],
                   properties->position_MBP[0],properties->position_MBP[1],properties->position_MBP[2],
                   properties->velocity_MBP[0],properties->velocity_MBP[1],properties->velocity_MBP[2],
                   properties->M_vir,
                   properties->R_vir,
                   properties->R_halo,
                   properties->R_max,
                   properties->V_max,
                   v_c,
                   properties->sigma_v,
                   properties->spin[0],properties->spin[1],properties->spin[2],
                   lambda,
                   properties->q_triaxial,
                   properties->s_triaxial,
                   offset_COM*1e3); // converts to kpc/h
            if(flag_use_profiles)
               fprintf(fp_out,"  %11.5f",overdensity);
            fprintf(fp_out,"\n");
         }

         // Close catalog
         SID_free(SID_FARG properties);
         if(flag_use_profiles)
            SID_free(SID_FARG profile);
         fclose(fp_out);
         fclose_catalog(&fp_catalog);
         SID_log("Done.",SID_LOG_CLOSE);
       }
     }
     SID_log("Done.",SID_LOG_CLOSE);
  }  

  SID_exit(ERROR_NONE);
}
