#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>

int main(int argc, char *argv[]){
  char    filename_root[256];
  char    filename_properties[256];
  char    filename_profiles[256];
  char    filename_out[256];
  char    prefix_text[5];
  FILE   *fp_properties=NULL;
  FILE   *fp_profiles  =NULL;
  FILE   *fp_out       =NULL;
  int     i_file;
  int     n_files;
  int     n_groups_all;
  int     i_group;
  int     i_group_selected;
  int     i_profile;
  int     flag_process_group;
  int     snap_number;
  int     n_groups_properties;
  int     n_groups_profiles;  
  halo_properties_info properties;
  halo_profile_info    profile;
  float  lambda,v_c;
  float  offset_COM;
  float  r_min,r_max;

  SID_init(&argc,&argv,NULL);

  strcpy(filename_root, argv[1]);
  snap_number     =atoi(argv[2]);
  i_group_selected=atoi(argv[3]);

  if(i_group_selected<0){
    flag_process_group=TRUE;
    i_group_selected*=-1;
    sprintf(prefix_text,"");
  }
  else{
    flag_process_group=FALSE;
    sprintf(prefix_text,"sub");
  }

  if(SID.I_am_Master){
    sprintf(filename_properties,"%s_%03d.catalog_%sgroups_properties",filename_root,snap_number,prefix_text);
    sprintf(filename_profiles,  "%s_%03d.catalog_%sgroups_profiles",  filename_root,snap_number,prefix_text);
    SID_log("Selected %sgroup properties file:{%s}",SID_LOG_COMMENT,prefix_text,filename_properties);
    SID_log("Selected %sgroup profiles   file:{%s}",SID_LOG_COMMENT,prefix_text,filename_profiles);
    SID_log("Selected %sgroup number         : %d", SID_LOG_COMMENT,prefix_text,i_group_selected);

    fp_properties=fopen(filename_properties,"r");
    fp_profiles  =fopen(filename_profiles,  "r");
    fread(&i_file,             sizeof(int),1,fp_properties);
    fread(&n_files,            sizeof(int),1,fp_properties);
    fread(&n_groups_properties,sizeof(int),1,fp_properties);
    fread(&n_groups_all,       sizeof(int),1,fp_properties);
    fread(&i_file,             sizeof(int),1,fp_profiles);
    fread(&n_files,            sizeof(int),1,fp_profiles);
    fread(&n_groups_profiles,  sizeof(int),1,fp_profiles);
    fread(&n_groups_all,       sizeof(int),1,fp_profiles);

    if(n_groups_properties!=n_groups_profiles)
      SID_trap_error("There's a mismatch in the number (ie %d!=%d) of groups in the 2 files.",ERROR_LOGIC,n_groups_properties,n_groups_profiles);
    else if(i_group_selected>=n_groups_properties)
      SID_trap_error("You have requested a %sgroup exceeding the number in the file (max=%d)",ERROR_LOGIC,prefix_text,n_groups_properties);
    else
      SID_log("Number of %sgroups in file:     : %d",SID_LOG_COMMENT,prefix_text,n_groups_properties);

    // Skip unwanted halos
    SID_log("Skipping halos...",SID_LOG_OPEN);
    for(i_group=0;i_group<i_group_selected;i_group++){
      fread(&(profile.n_bins),sizeof(int),1,fp_profiles);
      fseeko(fp_properties,sizeof(halo_properties_info),SEEK_CUR);
      fseeko(fp_profiles,  sizeof(halo_profile_bin_info)*profile.n_bins,SEEK_CUR);
    }
    SID_log("Done.",SID_LOG_CLOSE);

    SID_log("Reading selected halo...",SID_LOG_OPEN);
    // Read properties
    fread(&properties,sizeof(halo_properties_info),1,fp_properties);

    // Read profiles
    fread(&(profile.n_bins),sizeof(int),                  1,             fp_profiles);
    fread(profile.bins,     sizeof(halo_profile_bin_info),profile.n_bins,fp_profiles);
    SID_log("Done.",SID_LOG_CLOSE);

    // Close files
    fclose(fp_properties);
    fclose(fp_profiles);

    // Write output
    sprintf(filename_out,"%s_%03d.catalog_%sgroup.%09d",filename_root,snap_number,prefix_text,i_group);
    fp_out=fopen(filename_out,"w");
    v_c       =sqrt(G_NEWTON*properties.M_vir*M_SOL/(properties.R_vir*M_PER_MPC))*1e-3;
    lambda    =sqrt(properties.spin[0]*properties.spin[0]+properties.spin[1]*properties.spin[1]+properties.spin[2]*properties.spin[2])/(sqrt(2.)*properties.R_vir*v_c);
    offset_COM=sqrt(pow(properties.position_COM[0]-properties.position_MBP[0],2.)+pow(properties.position_COM[1]-properties.position_MBP[1],2.)+pow(properties.position_COM[2]-properties.position_MBP[2],2.));
    SID_log("Analysis of %sgroup #%d in snap #%d of %s",SID_LOG_COMMENT,prefix_text,i_group,snap_number,filename_root);
    SID_log("   id_MBP      =%14lld",                      SID_LOG_COMMENT,properties.id_MBP);
    SID_log("   n_particles =%14d",                        SID_LOG_COMMENT,properties.n_particles);
    SID_log("   position_COM=%14.6e %14.6e %14.6e [Mpc/h]",SID_LOG_COMMENT,properties.position_COM[0],properties.position_COM[1],properties.position_COM[2]);
    SID_log("   position_MBP=%14.6e %14.6e %14.6e [Mpc/h]",SID_LOG_COMMENT,properties.position_MBP[0],properties.position_MBP[1],properties.position_MBP[2]);
    SID_log("   velocity_COM=%14.6e %14.6e %14.6e [km/s]", SID_LOG_COMMENT,properties.velocity_COM[0],properties.velocity_COM[1],properties.velocity_COM[2]);
    SID_log("   velocity_MBP=%14.6e %14.6e %14.6e [km/s]", SID_LOG_COMMENT,properties.velocity_MBP[0],properties.velocity_MBP[1],properties.velocity_MBP[2]);
    SID_log("   M_vir       =%14.6e [M_sol/h]",            SID_LOG_COMMENT,properties.M_vir);
    SID_log("   R_vir       =%14.6e [Mpc/h]",              SID_LOG_COMMENT,properties.R_vir);
    SID_log("   R_halo      =%14.6e [Mpc/h]",              SID_LOG_COMMENT,properties.R_halo);
    SID_log("   R_max       =%14.6e [kpc/h]",              SID_LOG_COMMENT,properties.R_max*1e3);
    SID_log("   V_max       =%14.6e [km/s]",               SID_LOG_COMMENT,properties.V_max);
    SID_log("   V_vir       =%14.6e [km/s]",               SID_LOG_COMMENT,v_c);
    SID_log("   sigma_v     =%14.6e [km/s]",               SID_LOG_COMMENT,properties.sigma_v);
    SID_log("   spin        =%14.6e %14.6e %14.6e [Mpc/h km/s]",SID_LOG_COMMENT,properties.spin[0],properties.spin[1],properties.spin[2]);
    SID_log("   lambda      =%14.6e",                      SID_LOG_COMMENT,lambda);
    SID_log("   q_trixaial  =%14.6e",                      SID_LOG_COMMENT,properties.q_triaxial);
    SID_log("   s_triaxial  =%14.6e",                      SID_LOG_COMMENT,properties.s_triaxial);
    SID_log("   offset_COM  =%14.6e [kpc/h]",              SID_LOG_COMMENT,offset_COM*1e3);
    fprintf(fp_out,"# Analysis of %sgroup #%d in snap #%d of %s",prefix_text,i_group,snap_number,filename_root);
    fprintf(fp_out,"#   id_MBP      =%14lld\n",                      properties.id_MBP);
    fprintf(fp_out,"#   n_particles =%14d\n",                        properties.n_particles);
    fprintf(fp_out,"#   position_COM=%14.6e %14.6e %14.6e [Mpc/h]\n",properties.position_COM[0],properties.position_COM[1],properties.position_COM[2]);
    fprintf(fp_out,"#   position_MBP=%14.6e %14.6e %14.6e [Mpc/h]\n",properties.position_MBP[0],properties.position_MBP[1],properties.position_MBP[2]);
    fprintf(fp_out,"#   velocity_COM=%14.6e %14.6e %14.6e [km/s]\n", properties.velocity_COM[0],properties.velocity_COM[1],properties.velocity_COM[2]);
    fprintf(fp_out,"#   velocity_MBP=%14.6e %14.6e %14.6e [km/s]\n", properties.velocity_MBP[0],properties.velocity_MBP[1],properties.velocity_MBP[2]);
    fprintf(fp_out,"#   M_vir       =%14.6e [M_sol/h]\n",            properties.M_vir);
    fprintf(fp_out,"#   R_vir       =%14.6e [Mpc/h]\n",              properties.R_vir);
    fprintf(fp_out,"#   R_halo      =%14.6e [Mpc/h]\n",              properties.R_halo);
    fprintf(fp_out,"#   R_max       =%14.6e [kpc/h]\n",              properties.R_max*1e3);
    fprintf(fp_out,"#   V_max       =%14.6e [km/s]\n",               properties.V_max);
    fprintf(fp_out,"#   V_vir       =%14.6e [km/s]\n",               v_c);
    fprintf(fp_out,"#   sigma_v     =%14.6e [km/s]\n",               properties.sigma_v);
    fprintf(fp_out,"#   spin        =%14.6e %14.6e %14.6e [Mpc/h km/s]\n",properties.spin[0],properties.spin[1],properties.spin[2]);
    fprintf(fp_out,"#   lambda      =%14.6e\n",                      lambda);
    fprintf(fp_out,"#   q_trixaial  =%14.6e\n",                      properties.q_triaxial);
    fprintf(fp_out,"#   s_triaxial  =%14.6e\n",                      properties.s_triaxial);
    fprintf(fp_out,"#   offset_COM  =%14.6e [kpc/h]\n",              offset_COM*1e3);
    fprintf(fp_out,"# Profile columns: (1)     r_min       [Mpc/h]\n");
    fprintf(fp_out,"#                  (2)     r_med       [Mpc/h]\n");
    fprintf(fp_out,"#                  (3)     r_max       [Mpc/h]\n");
    fprintf(fp_out,"#                  (4)     n_particles \n");
    fprintf(fp_out,"#                  (5)     overdensity \n");
    fprintf(fp_out,"#                  (6)     v_c(r_max)  [km/s] \n");
    fprintf(fp_out,"#                  (7)     x_COM       [Mpc/h]\n");
    fprintf(fp_out,"#                  (8)     y_COM       [Mpc/h]\n");
    fprintf(fp_out,"#                  (9)     z_COM       [Mpc/h]\n");
    fprintf(fp_out,"#                  (10)    offset_COM  [Mpc/h]\n");
    fprintf(fp_out,"#                  (11)    M(<r_max)   [M_sol]\n");
    fprintf(fp_out,"#                  (12)    rho         [M_sol/Mpc^3]\n");
    fprintf(fp_out,"#                  (13)    v_c(r_max)  [km/s]\n");
    fprintf(fp_out,"#                  (14)    sigma_rad   [km/s]\n");
    fprintf(fp_out,"#                  (15)    sigma_tan   [km/s]\n");
    fprintf(fp_out,"#                  (16)    sigma_tot   [km/s]\n");
    fprintf(fp_out,"#                  (17)    beta\n");
    fprintf(fp_out,"#                  (18-20) spin        [Mpc/h km/s]\n");
    fprintf(fp_out,"#                  (21)    lambda\n");
    fprintf(fp_out,"#                  (22)    q_triaxial\n");
    fprintf(fp_out,"#                  (23)    s_triaxial\n");
    fprintf(fp_out,"#                  (24-32) shape_eigen_vectors\n");
    SID_log("Writing %d profile bins to output file {%s}...",SID_LOG_OPEN,profile.n_bins,filename_out);
    for(i_profile=0,r_max=0.;i_profile<profile.n_bins;i_profile++){
      r_min     =r_max;
      r_max     =profile.bins[i_profile].r_max;
      v_c       =sqrt(G_NEWTON*profile.bins[i_profile].M_r*M_SOL/(r_max*M_PER_MPC))*1e-3;
      lambda    =sqrt(profile.bins[i_profile].spin[0]*profile.bins[i_profile].spin[0]+profile.bins[i_profile].spin[1]*profile.bins[i_profile].spin[1]+profile.bins[i_profile].spin[2]*profile.bins[i_profile].spin[2])/(sqrt(2.)*v_c*profile.bins[i_profile].r_max);
      offset_COM=sqrt(pow(profile.bins[i_profile].position_COM[0],2.)+pow(profile.bins[i_profile].position_COM[1],2.)+pow(profile.bins[i_profile].position_COM[2],2.));
      fprintf(fp_out,"%e %e %e %9d %e %e %e %e %le %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
              r_min,
              profile.bins[i_profile].r_med,
              r_max,
              profile.bins[i_profile].n_particles,
              profile.bins[i_profile].overdensity,
              v_c,
              profile.bins[i_profile].position_COM[0],
              profile.bins[i_profile].position_COM[1],
              profile.bins[i_profile].position_COM[2],
              offset_COM,
              profile.bins[i_profile].M_r,
              profile.bins[i_profile].rho,
              sqrt(G_NEWTON*profile.bins[i_profile].M_r*M_SOL/(r_max*M_PER_MPC))*1e-3,
              profile.bins[i_profile].sigma_rad,
              profile.bins[i_profile].sigma_tan,
              profile.bins[i_profile].sigma_tot,
              1.-0.5*profile.bins[i_profile].sigma_tan/profile.bins[i_profile].sigma_rad,
              profile.bins[i_profile].spin[0],profile.bins[i_profile].spin[1],profile.bins[i_profile].spin[2],
              lambda,
              profile.bins[i_profile].q_triaxial,
              profile.bins[i_profile].s_triaxial,
              profile.bins[i_profile].shape_eigen_vectors[0][0],
              profile.bins[i_profile].shape_eigen_vectors[0][1],
              profile.bins[i_profile].shape_eigen_vectors[0][2],
              profile.bins[i_profile].shape_eigen_vectors[1][0],
              profile.bins[i_profile].shape_eigen_vectors[1][1],
              profile.bins[i_profile].shape_eigen_vectors[1][2],
              profile.bins[i_profile].shape_eigen_vectors[2][0],
              profile.bins[i_profile].shape_eigen_vectors[2][1],
              profile.bins[i_profile].shape_eigen_vectors[2][2]);
    }
    SID_log("Done.",SID_LOG_CLOSE);

    fclose(fp_out);
  }  

  SID_exit(ERROR_NONE);
}
