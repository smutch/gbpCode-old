#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHalos.h>

void print_results_local(FILE                 *fp_out,
                         int                   flag_use_profiles,
                         double                box_size,
                         halo_properties_info *properties,
                         halo_properties_info *properties_group,
                         halo_profile_info    *profile,
                         int                   i_group,
                         int                   j_subgroup,
                         int                   j_group);
void print_results_local(FILE                 *fp_out,
                         int                   flag_use_profiles,
                         double                box_size,
                         halo_properties_info *properties,
                         halo_properties_info *properties_group,
                         halo_profile_info    *profile,
                         int                   i_group,
                         int                   j_subgroup,
                         int                   j_group){
   float overdensity;
   if(flag_use_profiles){
      int n_bins=profile->n_bins;
      overdensity=profile->bins[n_bins-1].overdensity;
   }
   else
      overdensity=-1.;

   // Create a few properties
   double dx,dy,dz,v_c,lambda;
   float  offset_COM;
   float  r_min,r_max;
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
    if(properties_group!=properties){
       float f_sub=(float)(properties->M_vir/properties_group->M_vir);
       fprintf(fp_out,"  %9d %3d %11.5f",j_group,j_subgroup,f_sub);
    }
    else{
       fprintf(fp_out,"  %3d",j_subgroup);
    }
    fprintf(fp_out,"\n");
}

int main(int argc, char *argv[]){
  char    filename_properties[256];
  char    filename_profiles[256];
  char    filename_out_root[256];
  char    filename_out[256];
  char    filename_SSimPL[MAX_FILENAME_LENGTH];
  char    filename_halo_type[MAX_FILENAME_LENGTH];
  FILE   *fp_out[2];
  int     snap_number;
  int     snap_number_start;
  int     snap_number_stop;
  int     snap_number_step;
  double  box_size;

  SID_init(&argc,&argv,NULL,NULL);

  strcpy(filename_SSimPL,   argv[1]);
  strcpy(filename_halo_type,argv[2]);
  box_size            =atof(argv[3]);
  snap_number_start   =atoi(argv[4]);
  snap_number_stop    =atoi(argv[5]);
  snap_number_step    =atoi(argv[6]);
  strcpy(filename_out_root, argv[7]);

  int flag_use_profiles=FALSE;

  if(SID.I_am_Master){
    int i_type;
    SID_log("Processing catalogs for snaps %d->%d...",SID_LOG_OPEN|SID_LOG_TIMER,snap_number_start,snap_number_stop);
    SID_log("Properties structure size=%lld",SID_LOG_COMMENT,sizeof(halo_properties_info));
    for(snap_number=snap_number_start;snap_number<=snap_number_stop;snap_number++){

         // Open halos
         char filename_halos[256];
         sprintf(filename_halos,"%s/halos/%s_%03d.catalog_groups",filename_SSimPL,filename_halo_type,snap_number);
         FILE *fp_halos=NULL;
         if((fp_halos=fopen(filename_halos,"r"))==NULL)
            SID_trap_error("Could not open halo file {%s} for reading.",ERROR_IO_OPEN,filename_halos);
         int n_groups_halos,group_offset_byte_size;
         fread(&n_groups_halos,        sizeof(int),1,fp_halos);
         fread(&group_offset_byte_size,sizeof(int),1,fp_halos);
         // Skip group sizes and offsets
         fseeko(fp_halos,(off_t)(n_groups_halos*(sizeof(int)+group_offset_byte_size)),SEEK_CUR);

         // Open catalogs
         char filename_cat_root[256];
         sprintf(filename_cat_root,"%s/catalogs/%s",filename_SSimPL,filename_halo_type);
         fp_catalog_info fp_catalog_groups;
         fp_catalog_info fp_catalog_subgroups;
         fopen_catalog(filename_cat_root,
                       snap_number,
                       READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES|READ_CATALOG_PROPERTIES,
                       &fp_catalog_groups);
         fopen_catalog(filename_cat_root,
                       snap_number,
                       READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES|READ_CATALOG_PROPERTIES,
                       &fp_catalog_subgroups);

         // Sanity check
         if(n_groups_halos!=fp_catalog_groups.n_halos_total)
            SID_trap_error("Group counts in halo and catalog files don't match (ie. %d!=%d).",ERROR_LOGIC,n_groups_halos,fp_catalog_groups.n_halos_total);

         // Process halos
         SID_log("Writing snap #%03d...",SID_LOG_OPEN,snap_number);
         SID_log("(%d groups, %d subgroups)...",SID_LOG_CONTINUE,fp_catalog_groups.n_halos_total,fp_catalog_subgroups.n_halos_total);
         for(int i_type=0;i_type<2;i_type++){
            if(i_type==0)
               sprintf(filename_out,"%s_%03d_groups.ascii",filename_out_root,snap_number);
            else
               sprintf(filename_out,"%s_%03d_subgroups.ascii",filename_out_root,snap_number);
            fp_out[i_type]=fopen(filename_out,"w");
            int i_column=1;
            if(i_type==0){
               fprintf(fp_out[i_type],"# ASCII group catalog of snap #%d of %s\n",snap_number,filename_SSimPL);
               fprintf(fp_out[i_type],"# File columns: (%02d)     group number\n",      i_column++);
            }
            else{
               fprintf(fp_out[i_type],"# ASCII subgroup catalog of snap #%d of %s\n",snap_number,filename_SSimPL);
               fprintf(fp_out[i_type],"# File columns: (%02d)     subgroup number\n",      i_column++);
            }
            fprintf(fp_out[i_type],"#               (%02d)     # of particles\n",      i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     id_MBP\n",              i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     x_COM      [Mpc/h]\n",  i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     y_COM      [Mpc/h]\n",  i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     z_COM      [Mpc/h]\n",  i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     v_x_COM    [km/s]\n",   i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     v_y_COM    [km/s]\n",   i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     v_z_COM    [km/s]\n",   i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     x_MBP      [Mpc/h]\n",  i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     y_MBP      [Mpc/h]\n",  i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     z_MBP      [Mpc/h]\n",  i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     v_x_MBP    [km/s]\n",   i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     v_y_MBP    [km/s]\n",   i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     v_z_MBP    [km/s]\n",   i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     M_vir      [M_sol/h]\n",i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     R_vir      [Mpc/h]\n",  i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     R_halo     [Mpc/h]\n",  i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     R_max      [Mpc/h]\n",  i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     V_max      [km/s]\n",   i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     V_vir      [km/s]\n",   i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     sigma_v    [km/s]\n",   i_column++);
            fprintf(fp_out[i_type],"#               (%02d,%02d,%02d) spin       [Mpc/h km/s]\n",i_column,i_column+1,i_column+2);i_column+=3;
            fprintf(fp_out[i_type],"#               (%02d)     lambda (spin parameter)\n",i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     q_triaxial\n",             i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     s_triaxial\n",             i_column++);
            fprintf(fp_out[i_type],"#               (%02d)     offset_COM [kpc/h]\n",     i_column++);
            if(flag_use_profiles)
               fprintf(fp_out[i_type],"#               (%02d)     overdensity(R_halo)\n",i_column++);
            if(i_type==0){
               fprintf(fp_out[i_type],"#               (%02d)     n_substructures\n",     i_column++);
            }
            else{
               fprintf(fp_out[i_type],"#               (%02d)     group index\n",           i_column++);
               fprintf(fp_out[i_type],"#               (%02d)     substructure rank\n",     i_column++);
               fprintf(fp_out[i_type],"#               (%02d)     fraction of group mass\n",i_column++);
            }
         }

         halo_properties_info *properties_group   =(halo_properties_info *)SID_malloc(sizeof(halo_properties_info));
         halo_properties_info *properties_subgroup=(halo_properties_info *)SID_malloc(sizeof(halo_properties_info));
         halo_profile_info    *profile_group   =NULL;
         halo_profile_info    *profile_subgroup=NULL;
         if(flag_use_profiles){
            profile_group    = (halo_profile_info *)SID_malloc(sizeof(halo_profile_info));
            profile_subgroup = (halo_profile_info *)SID_malloc(sizeof(halo_profile_info));
         }

         // Perform read and write
         for(int i_group=0,i_subgroup=0;i_group<fp_catalog_groups.n_halos_total;i_group++){
           int n_subgroups_group;
           fread_catalog_file(&fp_catalog_groups,NULL,properties_group,profile_group,i_group);
           fread(&n_subgroups_group,sizeof(int),1,fp_halos);
           print_results_local(fp_out[0],flag_use_profiles,box_size,properties_group,properties_group,profile_group,i_group,n_subgroups_group,-1);
           for(int j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
              fread_catalog_file(&fp_catalog_subgroups,NULL,properties_subgroup,profile_subgroup,i_subgroup);
              print_results_local(fp_out[1],flag_use_profiles,box_size,properties_subgroup,properties_group,profile_subgroup,i_subgroup,j_subgroup,i_group);
           }
         }

         // Clean-up
         SID_free(SID_FARG properties_group);
         SID_free(SID_FARG properties_subgroup);
         if(flag_use_profiles){
            SID_free(SID_FARG profile_group);
            SID_free(SID_FARG profile_subgroup);
         }
         fclose(fp_out[0]);
         fclose(fp_out[1]);
         fclose(fp_halos);
         fclose_catalog(&fp_catalog_groups);
         fclose_catalog(&fp_catalog_subgroups);
         SID_log("Done.",SID_LOG_CLOSE);
     }
     SID_log("Done.",SID_LOG_CLOSE);
  }  

  SID_exit(ERROR_NONE);
}
