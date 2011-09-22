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

  SID_init(&argc,&argv,NULL);

  strcpy(filename_root,  argv[1]);
  snap_number_start=atoi(argv[2]);
  snap_number_stop =atoi(argv[3]);
  snap_number_step =atoi(argv[4]);

  sprintf(prefix_text,"sub");
  sprintf(prefix_text,"");

  if(SID.I_am_Master){

    SID_log("Processing %sgroup catalogs for snaps %d->%d...",SID_LOG_OPEN|SID_LOG_TIMER,prefix_text,snap_number_start,snap_number_stop);
    for(snap_number=snap_number_start;snap_number<=snap_number_stop;snap_number++){
      sprintf(filename_properties,"%s_%03d.catalog_%sgroups_properties",filename_root,snap_number,prefix_text);
      sprintf(filename_profiles,  "%s_%03d.catalog_%sgroups_profiles",  filename_root,snap_number,prefix_text);
      sprintf(filename_out,"%s.ascii",filename_properties);
      SID_log("Writing snap #%03d %sgroup properties to file {%s}...",SID_LOG_OPEN,snap_number,prefix_text,filename_out);

      fp_properties=fopen(filename_properties,"r");
      fread(&i_file,             sizeof(int),1,fp_properties);
      fread(&n_files,            sizeof(int),1,fp_properties);
      fread(&n_groups_properties,sizeof(int),1,fp_properties);
      fread(&n_groups_all,       sizeof(int),1,fp_properties);
      SID_log("(%d %sgroups)...",SID_LOG_CONTINUE,n_groups_properties,prefix_text);

      fp_profiles=fopen(filename_profiles,"r");
      fread(&i_file_profiles,       sizeof(int),1,fp_profiles);
      fread(&n_files_profiles,      sizeof(int),1,fp_profiles);
      fread(&n_groups_profiles,     sizeof(int),1,fp_profiles);
      fread(&n_groups_profiles_all, sizeof(int),1,fp_profiles);

      // Process halos
      fp_out=fopen(filename_out,"w");
      fprintf(fp_out,"# Ascii %sgroup catalog of snap #%d of %s\n",prefix_text,snap_number,filename_root);
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
      fprintf(fp_out,"#               (19)    R_max      [kpc/h]\n");
      fprintf(fp_out,"#               (20)    V_max      [km/s]\n");
      fprintf(fp_out,"#               (21)    V_vir      [km/s]\n");
      fprintf(fp_out,"#               (22)    sigma_v    [km/s]\n");
      fprintf(fp_out,"#               (23-25) spin       [Mpc/h km/s]\n");
      fprintf(fp_out,"#               (26)    lambda\n");
      fprintf(fp_out,"#               (27)    q_triaxial\n");
      fprintf(fp_out,"#               (28)    s_triaxial\n");
      fprintf(fp_out,"#               (29)    offset_COM [Mpc/h]\n");
      fprintf(fp_out,"#               (30)    overdensity(R_halo)\n");
      for(i_group=0;i_group<n_groups_properties;i_group++){
        fread(&properties,sizeof(halo_properties_info),1,fp_properties);

        fread(&(profile.n_bins),sizeof(int),                 1,              fp_profiles);
        fread(&(profile.bins),  sizeof(halo_profile_bin_info),profile.n_bins,fp_profiles);

        // Create a few properties
        v_c       =sqrt(G_NEWTON*properties.M_vir*M_SOL/(properties.R_vir*M_PER_MPC))*1e-3;
        lambda    =sqrt(properties.spin[0]*properties.spin[0]+properties.spin[1]*properties.spin[1]+properties.spin[2]*properties.spin[2])/(sqrt(2.)*properties.R_vir*v_c);
        offset_COM=sqrt(pow(properties.position_COM[0]-properties.position_MBP[0],2.)+pow(properties.position_COM[1]-properties.position_MBP[1],2.)+pow(properties.position_COM[2]-properties.position_MBP[2],2.));

        // Perform write
        fprintf(fp_out,"%9d %9d %9lld  %12.5e %12.5e %12.5e %11.5f %11.5f %11.5f  %12.5e %12.5e %12.5e %11.5f %11.5f %11.5f  %10.5le %11.5f %11.5f %11.5f  %11.5f %11.5f %11.5f  %12.5e %12.5e %12.5e %11.5f  %11.5f %11.5f  %11.5f  %11.5f\n",
                i_group,properties.n_particles,properties.id_MBP,
                properties.position_COM[0],properties.position_COM[1],properties.position_COM[2],
                properties.velocity_COM[0],properties.velocity_COM[1],properties.velocity_COM[2],
                properties.position_MBP[0],properties.position_MBP[1],properties.position_MBP[2],
                properties.velocity_MBP[0],properties.velocity_MBP[1],properties.velocity_MBP[2],
                properties.M_vir,
                properties.R_vir,
                properties.R_halo,
                properties.R_max*1e3,
                properties.V_max,
                v_c,
                properties.sigma_v,
                properties.spin[0],properties.spin[1],properties.spin[2],
                lambda,
                properties.q_triaxial,
                properties.s_triaxial,
                offset_COM*1e3,
                profile.bins[profile.n_bins-1].overdensity);
      }

      // Close files
      fclose(fp_profiles);
      fclose(fp_properties);
      fclose(fp_out);
      SID_log("Done.",SID_LOG_CLOSE);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }  

  SID_exit(ERROR_NONE);
}
