#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>

int main(int argc, char *argv[]){
  char    filename_properties[256];
  char    filename_profiles[256];
  char    filename_out_root[256];
  char    filename_out[256];
  char    filename_SSimPL[MAX_FILENAME_LENGTH];
  char    filename_halo_type[MAX_FILENAME_LENGTH];
  int     snap_number;
  int     snap_number_start;
  int     snap_number_stop;
  int     snap_number_step;

  SID_init(&argc,&argv,NULL,NULL);

  strcpy(filename_SSimPL,   argv[1]);
  strcpy(filename_halo_type,argv[2]);
  snap_number_start   =atoi(argv[3]);
  snap_number_stop    =atoi(argv[4]);
  snap_number_step    =atoi(argv[5]);
  strcpy(filename_out_root, argv[6]);

  int flag_use_profiles=FALSE;

  if(SID.I_am_Master){
    SID_log("Processing catalogs for snaps %d->%d...",SID_LOG_OPEN|SID_LOG_TIMER,snap_number_start,snap_number_stop);
    for(snap_number=snap_number_start;snap_number<=snap_number_stop;snap_number++){

         // Open halos
         char filename_halos[256];
         sprintf(filename_halos,"%s/halos/%s_%03d.catalog_groups",filename_SSimPL,filename_halo_type,snap_number);
         FILE *fp_halos=NULL;
         if((fp_halos=fopen(filename_halos,"r"))==NULL)
            SID_trap_error("Could not open halo file {%s} for reading.",ERROR_IO_OPEN,filename_halos);
         int n_groups_halos,group_offset_byte_size;
         fread_verify(&n_groups_halos,        sizeof(int),1,fp_halos);
         fread_verify(&group_offset_byte_size,sizeof(int),1,fp_halos);

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

         // Open SO files if they're available
         fp_multifile_info fp_SO;
         int flag_use_SO=fopen_multifile("%s/catalogs/%s_%03d.catalog_groups_SO",sizeof(float),&fp_SO,filename_SSimPL,filename_halo_type,snap_number);
         if(flag_use_SO)
            SID_log("SO files present.",SID_LOG_COMMENT);

         // Sanity check
         if(n_groups_halos!=fp_catalog_groups.n_halos_total)
            SID_trap_error("Group counts in halo and catalog files don't match (ie. %d!=%d).",ERROR_LOGIC,n_groups_halos,fp_catalog_groups.n_halos_total);

         // Process halos
         SID_log("Processing snapshot #%03d...",SID_LOG_OPEN,snap_number);
         SID_log("(%d groups, %d subgroups)...",SID_LOG_CONTINUE,fp_catalog_groups.n_halos_total,fp_catalog_subgroups.n_halos_total);

         // Initialzie halo trend data structure
         halo_trend_info  halo_trend_data;
         char             filename_run[MAX_FILENAME_LENGTH];
         sprintf(filename_run,"%s/run/run.txt",filename_SSimPL);
         parameter_list_info *parameter_list=NULL;
         init_parameter_list(&parameter_list);
         add_parameter_to_list(parameter_list,"box_size",SID_DOUBLE,   PARAMETER_MODE_DEFAULT);
         add_parameter_to_list(parameter_list,"N_dark",  SID_SIZE_T,   PARAMETER_MODE_DEFAULT);
         add_parameter_to_list(parameter_list,"m_dark",  SID_DOUBLE,   PARAMETER_MODE_DEFAULT);
         read_gbpParam_file(filename_run,parameter_list);
         fetch_parameter_data(parameter_list,"box_size",&(halo_trend_data.box_size)); 
         fetch_parameter_data(parameter_list,"m_dark",  &(halo_trend_data.m_p)); 
         free_parameter_list(&parameter_list);
         char             filename_snaps[MAX_FILENAME_LENGTH];
         sprintf(filename_snaps,"%s/run/a_list.txt",filename_SSimPL);
         FILE *fp_snaps=fopen(filename_snaps,"r");
         size_t line_length=0;
         char  *line=NULL;
         halo_trend_data.n_snaps=count_lines_data(fp_snaps);
         halo_trend_data.z_list =(double *)SID_malloc(sizeof(double)*halo_trend_data.n_snaps);
         for (int i_snap=0;i_snap<halo_trend_data.n_snaps;i_snap++){
            double a_i;
            grab_next_line_data(fp_snaps,&line,&line_length);
            grab_double(line,1,&a_i);
            halo_trend_data.z_list[i_snap]=z_of_a(a_i);
         }
         SID_free(SID_FARG line);
         fclose(fp_snaps);

         // Initialize halo data structure
         halo_info        halo_data;
         halo_data.flag_use_profiles  = flag_use_profiles;
         halo_data.flag_use_SO        = flag_use_SO;
         halo_data.snapshot           = snap_number;
         halo_data.properties_group   =(halo_properties_info *)SID_malloc(sizeof(halo_properties_info));
         halo_data.properties_subgroup=(halo_properties_info *)SID_malloc(sizeof(halo_properties_info));
         halo_data.profiles_group     =(halo_profile_info    *)SID_malloc(sizeof(halo_profile_info));
         halo_data.profiles_subgroup  =(halo_profile_info    *)SID_malloc(sizeof(halo_profile_info));

         // Initialize trends
         trend_info *trend_M_FoF=NULL;
         init_trend(&trend_M_FoF,"SSFctn",&halo_trend_data,init_halo_trend_property_logM_FoF,free_halo_trend_property_logM_FoF,calc_halo_trend_property_index_logM_FoF);
         init_halo_trend_coordinate(&halo_trend_data,trend_M_FoF,"SSFctn");

         // Read halos and construct histograms
         for(int i_group=0,i_subgroup=0;i_group<fp_catalog_groups.n_halos_total;i_group++){
            int n_subgroups_group;
            // Read group catalog
            fread_catalog_file(&fp_catalog_groups,NULL,halo_data.properties_group,halo_data.profiles_group,i_group);
            // Read number of subgroups
            fread_verify(&n_subgroups_group,sizeof(int),1,fp_halos);
            // Read SO masses (if available)
            if(flag_use_SO)
               fread_multifile(&fp_SO,halo_data.SO_data_group,i_group);
            // Loop over subgroups
            halo_data.n_sub         =n_subgroups_group;
            halo_data.np_sub        =0;
            halo_data.np_sub_largest=0;
            for(int j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
               // Read subgroup properties
               fread_catalog_file(&fp_catalog_subgroups,NULL,halo_data.properties_subgroup,halo_data.profiles_subgroup,i_subgroup);
               int np_i=halo_data.properties_subgroup->n_particles;
               halo_data.np_sub+=np_i;
               if(np_i>halo_data.np_sub_largest)
                  halo_data.np_sub_largest=np_i;
               // Add halo to subgroup trends
            }
            // Add halo to group trends
            add_item_to_trend(trend_M_FoF,GBP_ADD_ITEM_TO_TREND_DEFAULT,&halo_data);
         }

         // Write results
         char filename_out[MAX_FILENAME_LENGTH];
         sprintf(filename_out,"%s_%03d",filename_out_root,snap_number);
         write_trend_ascii(trend_M_FoF,filename_out);
         free_trend(&trend_M_FoF);

         // Clean-up
         SID_free(SID_FARG halo_trend_data.z_list);
         SID_free(SID_FARG halo_data.properties_group);
         SID_free(SID_FARG halo_data.properties_subgroup);
         SID_free(SID_FARG halo_data.profiles_group);
         SID_free(SID_FARG halo_data.profiles_subgroup);
         fclose(fp_halos);
         fclose_catalog(&fp_catalog_groups);
         fclose_catalog(&fp_catalog_subgroups);
         fclose_multifile(&fp_SO);
         SID_log("Done.",SID_LOG_CLOSE);
     }
     SID_log("Done.",SID_LOG_CLOSE);
  }  

  SID_exit(ERROR_NONE);
}
