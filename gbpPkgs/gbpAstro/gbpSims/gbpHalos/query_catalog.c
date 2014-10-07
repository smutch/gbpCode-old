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

    fp_catalog_info fp_group_properties;
    if(flag_process_group)
      fopen_catalog(filename_root,
                    snap_number,
                    READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES,
                    &fp_group_properties);
    else
      fopen_catalog(filename_root,
                    snap_number,
                    READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES,
                    &fp_group_properties);

    SID_log("Number of %sgroups in file:     : %d",SID_LOG_COMMENT,prefix_text,fp_group_properties.n_halos_total);

    // Skip unwanted halos
    //SID_log("Skipping halos...",SID_LOG_OPEN);
    //for(i_group=0;i_group<i_group_selected;i_group++){
    //  fread(&(profile.n_bins),sizeof(int),1,fp_profiles);
    //  fseeko(fp_properties,sizeof(halo_properties_info),SEEK_CUR);
    //  fseeko(fp_profiles,  sizeof(halo_profile_bin_info)*profile.n_bins,SEEK_CUR);
    //}
    //SID_log("Done.",SID_LOG_CLOSE);

    // Perform read 
    SID_log("Reading selected halo...",SID_LOG_OPEN);
    fread_catalog_file(&fp_group_properties,NULL,&properties,NULL,i_group_selected);
    SID_log("Done.",SID_LOG_CLOSE);

    // Close file
    fclose_catalog(&fp_group_properties);

    // Write output
    v_c       =sqrt(G_NEWTON*properties.M_vir*M_SOL/(properties.R_vir*M_PER_MPC))*1e-3;
    lambda    =sqrt(properties.spin[0]*properties.spin[0]+properties.spin[1]*properties.spin[1]+properties.spin[2]*properties.spin[2])/(sqrt(2.)*properties.R_vir*v_c);
    offset_COM=sqrt(pow(properties.position_COM[0]-properties.position_MBP[0],2.)+pow(properties.position_COM[1]-properties.position_MBP[1],2.)+pow(properties.position_COM[2]-properties.position_MBP[2],2.));
    SID_log("Analysis of %sgroup #%d in snap #%d of %s",SID_LOG_COMMENT,prefix_text,i_group_selected,snap_number,filename_root);
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
    SID_log("Done.",SID_LOG_CLOSE);
  }  

  SID_exit(ERROR_NONE);
}
