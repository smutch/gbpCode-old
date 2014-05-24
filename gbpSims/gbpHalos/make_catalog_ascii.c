#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHalos.h>

int main(int argc, char *argv[]){
  char    filename_root[256];
  char    filename_properties[256];
  char    filename_profiles[256];
  char    filename_out[256];
  char    prefix_text[5];
  FILE   *fp_profiles  =NULL;
  FILE   *fp_properties=NULL;
  FILE   *fp_out       =NULL;
  int     i_file;
  int     i_file_profiles;
  int     n_files;
  int     n_files_profiles;
  int     n_groups_all;
  int     i_group;
  int     i_group_selected;
  int     snap_number;
  int     snap_number_start;
  int     snap_number_stop;
  int     snap_number_step;
  int     n_groups_properties;
  int     n_groups_profiles;
  int     n_groups_profiles_all;
  halo_properties_info properties;
  halo_profile_info    profile;
  float  lambda,v_c;
  float  offset_COM;
  float  r_min,r_max;
  double box_size;

  SID_init(&argc,&argv,NULL);

  strcpy(filename_root,  argv[1]);
  box_size         =atof(argv[2]);
  snap_number_start=atoi(argv[3]);
  snap_number_stop =atoi(argv[4]);
  snap_number_step =atoi(argv[5]);

  int flag_use_profiles=FALSE;

  if(SID.I_am_Master){
    int i_type;
    SID_log("Processing catalogs for snaps %d->%d...",SID_LOG_OPEN|SID_LOG_TIMER,snap_number_start,snap_number_stop);
    SID_log("Properties structure size=%lld",SID_LOG_COMMENT,sizeof(halo_properties_info));
    for(snap_number=snap_number_start;snap_number<=snap_number_stop;snap_number++){
       for(i_type=0;i_type<2;i_type++){
         switch(i_type){
            case 0:
               sprintf(prefix_text,"");
               break;
            case 1:
               sprintf(prefix_text,"sub");
               break;
         }
         int flag_multifile=FALSE;
         int flag_notfound =TRUE;
         sprintf(filename_properties,"%s_%03d.catalog_%sgroups_properties",filename_root,snap_number,prefix_text);
         sprintf(filename_out,       "%s.ascii",filename_properties);
         if(flag_use_profiles)
            sprintf(filename_profiles,  "%s_%03d.catalog_%sgroups_profiles",  filename_root,snap_number,prefix_text);
         fp_properties=fopen(filename_properties,"r");
         SID_log("Writing snap #%03d %sgroup properties to file {%s}...",SID_LOG_OPEN,snap_number,prefix_text,filename_out);

         fread(&i_file,             sizeof(int),1,fp_properties);
         fread(&n_files,            sizeof(int),1,fp_properties);
         fread(&n_groups_properties,sizeof(int),1,fp_properties);
         fread(&n_groups_all,       sizeof(int),1,fp_properties);
         SID_log("(%d %sgroups)...",SID_LOG_CONTINUE,n_groups_properties,prefix_text);

         if(flag_use_profiles){
            fp_profiles=fopen(filename_profiles,"r");
            fread(&i_file_profiles,       sizeof(int),1,fp_profiles);
            fread(&n_files_profiles,      sizeof(int),1,fp_profiles);
            fread(&n_groups_profiles,     sizeof(int),1,fp_profiles);
            fread(&n_groups_profiles_all, sizeof(int),1,fp_profiles);
         }

         // Process halos
         fp_out=fopen(filename_out,"w");
         fprintf(fp_out,"# ASCII %sgroup catalog of snap #%d of %s\n",prefix_text,snap_number,filename_root);
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
         for(i_group=0;i_group<n_groups_properties;i_group++){
           fread(&properties,sizeof(halo_properties_info),1,fp_properties);

           float overdensity;
           if(flag_use_profiles){
              fread(&(profile.n_bins),sizeof(int),                 1,              fp_profiles);
              fread(&(profile.bins),  sizeof(halo_profile_bin_info),profile.n_bins,fp_profiles);
              overdensity=profile.bins[profile.n_bins-1].overdensity;
           }
           else
              overdensity=-1.;

           // Create a few properties
           double dx,dy,dz;
           v_c       =sqrt(G_NEWTON*properties.M_vir*M_SOL/(properties.R_vir*M_PER_MPC))*1e-3;
           lambda    =sqrt(properties.spin[0]*properties.spin[0]+properties.spin[1]*properties.spin[1]+properties.spin[2]*properties.spin[2])/(sqrt(2.)*properties.R_vir*v_c);
           dx        =d_periodic(properties.position_COM[0]-properties.position_MBP[0],box_size);
           dy        =d_periodic(properties.position_COM[1]-properties.position_MBP[1],box_size);
           dz        =d_periodic(properties.position_COM[2]-properties.position_MBP[2],box_size);
           offset_COM=sqrt(dx*dx+dy*dy+dz*dz);

           // Perform write
           fprintf(fp_out,"%9d %9d %9lld  %12.5e %12.5e %12.5e %11.5f %11.5f %11.5f  %12.5e %12.5e %12.5e %11.5f %11.5f %11.5f  %10.5le %11.5f %11.5f %11.5f  %11.5f %11.5f %11.5f  %12.5e %12.5e %12.5e %11.5f  %11.5f %11.5f  %11.5f",
                   i_group,properties.n_particles,properties.id_MBP,
                   properties.position_COM[0],properties.position_COM[1],properties.position_COM[2],
                   properties.velocity_COM[0],properties.velocity_COM[1],properties.velocity_COM[2],
                   properties.position_MBP[0],properties.position_MBP[1],properties.position_MBP[2],
                   properties.velocity_MBP[0],properties.velocity_MBP[1],properties.velocity_MBP[2],
                   properties.M_vir,
                   properties.R_vir,
                   properties.R_halo,
                   properties.R_max,
                   properties.V_max,
                   v_c,
                   properties.sigma_v,
                   properties.spin[0],properties.spin[1],properties.spin[2],
                   lambda,
                   properties.q_triaxial,
                   properties.s_triaxial,
                   offset_COM*1e3); // converts to kpc/h
            if(flag_use_profiles)
               fprintf(fp_out,"  %11.5f",overdensity);
            fprintf(fp_out,"\n");
         }

         // Close files
         if(flag_use_profiles)
            fclose(fp_profiles);
         fclose(fp_properties);
         fclose(fp_out);
         SID_log("Done.",SID_LOG_CLOSE);
       }
     }
     SID_log("Done.",SID_LOG_CLOSE);
  }  

  SID_exit(ERROR_NONE);
}
