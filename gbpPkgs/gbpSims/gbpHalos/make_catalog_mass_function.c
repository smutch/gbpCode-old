#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHalos.h>

int main(int argc, char *argv[]){
  char    filename_properties[256];
  char    filename_profiles[256];
  char    filename_out_root[256];
  char    filename_out[256];
  char    filename_SSimPL[MAX_FILENAME_LENGTH];
  char    filename_halo_type[MAX_FILENAME_LENGTH];
  FILE   *fp_out[2];

  SID_init(&argc,&argv,NULL,NULL);

  strcpy(filename_SSimPL,               argv[1]);
  strcpy(filename_halo_type,            argv[2]);
  int    snap_number_start=(int)   atoi(argv[3]);
  int    snap_number_stop =(int)   atoi(argv[4]);
  int    snap_number_step =(int)   atoi(argv[5]);
  double lM_min           =(double)atof(argv[6]);
  double dlM              =(double)atof(argv[7]);
  int    n_bins           =(int)   atoi(argv[8]);
  strcpy(filename_out_root,             argv[9]);

  SID_log("Generating mass functions for catalogs {%s;%s}, snaps %d->%d...",SID_LOG_OPEN|SID_LOG_TIMER,
          filename_SSimPL,filename_halo_type,snap_number_start,snap_number_stop);

  // Initialize cosmology
  cosmo_info *cosmo=NULL;
  char filename_cosmology[MAX_FILENAME_LENGTH];
  sprintf(filename_cosmology,"%s/run",filename_SSimPL);
  read_gbpCosmo_file(&cosmo,filename_cosmology);
  double h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];

  // Read parameter file
  double box_size;
  double m_dark;
  char   filename_run[MAX_FILENAME_LENGTH];
  sprintf(filename_run,"%s/run/run.txt",filename_SSimPL);
  parameter_list_info *parameter_list=NULL;
  init_parameter_list(&parameter_list);
  add_parameter_to_list(parameter_list,"box_size",SID_DOUBLE,   PARAMETER_MODE_DEFAULT);
  add_parameter_to_list(parameter_list,"N_dark",  SID_SIZE_T,   PARAMETER_MODE_DEFAULT);
  add_parameter_to_list(parameter_list,"m_dark",  SID_DOUBLE,   PARAMETER_MODE_DEFAULT);
  read_parameter_file(filename_run,parameter_list);
  fetch_parameter_data(parameter_list,"box_size",&box_size);
  fetch_parameter_data(parameter_list,"m_dark",  &m_dark);
  free_parameter_list(&parameter_list);
  SID_log("Box size     = %.2lf Mpc/h",SID_LOG_COMMENT,box_size);

  // Read list of redshifts
  char   filename_alist[MAX_FILENAME_LENGTH];
  sprintf(filename_alist,"%s/run/a_list.txt",filename_SSimPL);
  FILE   *fp_snaps=fopen(filename_alist,"r");
  size_t  line_length=0;
  char   *line=NULL;
  int     n_snaps=count_lines_data(fp_snaps);
  double *z_list =(double *)SID_malloc(sizeof(double)*n_snaps);
  for (int i_snap=0;i_snap<n_snaps;i_snap++){
     double a_i;
     grab_next_line_data(fp_snaps,&line,&line_length);
     grab_double(line,1,&a_i);
     z_list[i_snap]=z_of_a(a_i);
  }
  SID_free(SID_FARG line);
  fclose(fp_snaps);
  SID_log("No. of snaps = %d (z=%.2lf to %.2lf)",SID_LOG_COMMENT,n_snaps,z_list[0],z_list[n_snaps-1]);

  int flag_use_profiles=FALSE;

  if(SID.I_am_Master){
    int i_type;
    for(int snap_number=snap_number_start;snap_number<=snap_number_stop;snap_number++){
         double redshift=z_list[snap_number];
         SID_log("Processing snap #%03d (z=%.2lf)...",SID_LOG_OPEN,snap_number,redshift);

         // Open halos
         SID_log("Reading...",SID_LOG_OPEN|SID_LOG_TIMER,snap_number,redshift);
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
         //sprintf(filename_cat_root,"%s/catalogs_reshape/%s",filename_SSimPL,filename_halo_type);
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
         float *group_SO_data=NULL;
         fp_multifile_info fp_SO;
         int flag_use_SO=fopen_multifile("%s/catalogs/%s_%03d.catalog_groups_SO",sizeof(float),&fp_SO,filename_SSimPL,filename_halo_type,snap_number);
         if(flag_use_SO){
            group_SO_data=(float *)SID_malloc(sizeof(float));
            SID_log("SO files present.",SID_LOG_COMMENT);
         }

         // Sanity check
         if(n_groups_halos!=fp_catalog_groups.n_halos_total)
            SID_trap_error("Group counts in halo and catalog files don't match (ie. %d!=%d).",ERROR_LOGIC,n_groups_halos,fp_catalog_groups.n_halos_total);

         // Process halos
         SID_log("(%d groups, %d subgroups)...",SID_LOG_CONTINUE,fp_catalog_groups.n_halos_total,fp_catalog_subgroups.n_halos_total);

         halo_properties_info *properties_group   =(halo_properties_info *)SID_malloc(sizeof(halo_properties_info));
         halo_properties_info *properties_subgroup=(halo_properties_info *)SID_malloc(sizeof(halo_properties_info));
         halo_profile_info    *profile_group   =NULL;
         halo_profile_info    *profile_subgroup=NULL;
         if(flag_use_profiles){
            profile_group    = (halo_profile_info *)SID_malloc(sizeof(halo_profile_info));
            profile_subgroup = (halo_profile_info *)SID_malloc(sizeof(halo_profile_info));
         }

         // Perform read
         int     n_data=fp_catalog_groups.n_halos_total;
         double *data  =(double *)SID_malloc(sizeof(double)*n_data);
         double *data_s=(double *)SID_calloc(sizeof(double)*n_data);
         for(int i_group=0,i_subgroup=0;i_group<n_data;i_group++){
            int n_subgroups_group;
            // Read group catalog
            fread_catalog_file(&fp_catalog_groups,NULL,properties_group,profile_group,i_group);
            // Read number of subgroups
            fread(&n_subgroups_group,sizeof(int),1,fp_halos);
            // Read SO masses (if available)
            if(group_SO_data!=NULL)
               fread_multifile(&fp_SO,group_SO_data,i_group);
            // Process group information
            //double M_vir_FoF=properties_group->M_vir;
            double M_vir_FoF=m_dark*(double)properties_subgroup->n_particles*(1.-pow((double)properties_subgroup->n_particles,-0.6));
            data[i_group]=take_log10(M_vir_FoF);

            // Loop over subgroups
            for(int j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
               // Read subgroup properties
               fread_catalog_file(&fp_catalog_subgroups,NULL,properties_subgroup,profile_subgroup,i_subgroup);
               // Process subgroup information
               double M_vir_sub=properties_subgroup->M_vir;
               //double M_vir_sub=m_dark*(double)properties_subgroup->n_particles*(1.-pow((double)properties_subgroup->n_particles,-0.6));
               data_s[i_group]+=M_vir_sub;
            }
            data_s[i_group]=take_log10(data_s[i_group]);
         }
         SID_log("Done.",SID_LOG_CLOSE);

         // Perform sort
         SID_log("Sorting data...",SID_LOG_OPEN|SID_LOG_TIMER);
         merge_sort(data,n_data,NULL,SID_DOUBLE,SORT_INPLACE_ONLY,SORT_COMPUTE_INPLACE);
         SID_log("Done.",SID_LOG_CLOSE);

         // Build mass functions
         SID_log("Computing mass function...",SID_LOG_OPEN|SID_LOG_TIMER);
         double *bin       =(double *)SID_malloc(sizeof(double)*(n_bins+1));
         double *bin_median=(double *)SID_malloc(sizeof(double)*n_bins);
         int    *hist      =(int    *)SID_calloc(sizeof(int)   *n_bins);
         double  lM_bin_min=lM_min;
         double  lM_bin_max=lM_min;
         int     i_data_lo=-1;
         int     i_data_hi=-1;
         int     i_bin =0;
         int     i_data=0;

         if(n_data>0) while(data[i_data]<lM_min && i_data<n_data) i_data++;

         for(i_bin=0;i_bin<n_bins;i_bin++){
           lM_bin_min=lM_bin_max;
           lM_bin_max=lM_min+((double)(i_bin+1))*dlM;
           bin[i_bin]=lM_bin_min;
           i_data_lo=i_data;
           i_data_hi=i_data;
           if(i_data<n_data){
              while(data[i_data]<lM_bin_max && i_data<n_data){
                 hist[i_bin]++;
                 i_data_hi=i_data;
                 i_data++;
                 if(i_data>=n_data) break;
              }
           }
           int i_data_mid=(i_data_lo+i_data_hi)/2;
           if(hist[i_bin]>0){
             if(hist[i_bin]%2)
               bin_median[i_bin]=data[i_data_mid];
             else
               bin_median[i_bin]=0.5*(data[i_data_mid]+data[i_data_mid+1]);
           }
           else
              bin_median[i_bin]=0.5*(lM_bin_max+lM_bin_min);
         }
         bin[i_bin]=lM_bin_max;
         SID_log("Done.",SID_LOG_CLOSE);

         // Write mass function
         sprintf(filename_out,"%s_%03d.txt",filename_out_root,snap_number);
         FILE *fp_out;
         if((fp_out=fopen(filename_out,"w"))==NULL){
           fprintf(stderr,"Error opening output file {%s}.\n",filename_out);
           SID_free(SID_FARG data);
           return(1);
         }
         SID_log("Writing results to {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_out);
         double box_volume=box_size*box_size*box_size;
         fprintf(fp_out,"# Mass function for snapshot %d of catalog {%s;%s}\n",snap_number,filename_SSimPL,filename_halo_type);
         fprintf(fp_out,"# Column (01): M_lo     [source units]\n");
         fprintf(fp_out,"#        (02): M_median [source units]\n");
         fprintf(fp_out,"#        (03): M_hi     [source units]\n");
         fprintf(fp_out,"#        (04): No. in bin\n");
         fprintf(fp_out,"#        (05): MFn (per unit volume, per dlogM)\n");
         fprintf(fp_out,"#        (06): +/- MFn\n");
         fprintf(fp_out,"#        (07): Sheth & Tormen MFn\n");
         fprintf(fp_out,"#        (08): Watson MFn\n");
         fprintf(fp_out,"#        (09): No. w/ M>M_lo\n");
         fprintf(fp_out,"#        (10): Cumulative MFn (per unit volume)\n");
         fprintf(fp_out,"#        (11): +/- Cumulative MFn\n");
         fprintf(fp_out,"#        (12): Sheth & Tormen Cumulative MFn\n");
         fprintf(fp_out,"#        (13): Watson Cumulative MFn\n");
         double M_sol_inv_h=M_SOL/h_Hubble;
         double Mpc_inv_h  =M_PER_MPC/h_Hubble;
         for(int i=0;i<n_bins;i++){
           double dn_dlogM_theory_1=mass_function(take_alog10(bin_median[i])*M_sol_inv_h,
                                                  redshift,
                                                  &cosmo,
                                                  MF_ST)*pow(Mpc_inv_h,3.0);
           double n_theory_1=mass_function_cumulative(take_alog10(bin[i])*M_sol_inv_h,
                                                      redshift,
                                                      &cosmo,
                                                      MF_ST)*pow(Mpc_inv_h,3.0);
           double dn_dlogM_theory_2=mass_function(take_alog10(bin_median[i])*M_sol_inv_h,
                                                  redshift,
                                                  &cosmo,
                                                  MF_WATSON)*pow(Mpc_inv_h,3.0);
           double n_theory_2=mass_function_cumulative(take_alog10(bin[i])*M_sol_inv_h,
                                                      redshift,
                                                      &cosmo,
                                                      MF_WATSON)*pow(Mpc_inv_h,3.0);
           // Compute cumulative histogram
           int cumulative_hist=0;
           for(int j_bin=i;j_bin<n_bins;j_bin++)
              cumulative_hist+=hist[j_bin];
           fprintf(fp_out,"%11.4le %11.4le %11.4le %6d %11.4le %11.4le %10.4le %10.4le %6d %10.4le %10.4le %10.4le %10.4le\n",
                   bin[i],
                   bin_median[i],
                   bin[i+1],
                   hist[i],
                   (double)(hist[i])/(box_volume*dlM),
                   sqrt((double)(hist[i]))/(box_volume*dlM),
                   dn_dlogM_theory_1,dn_dlogM_theory_2,
                   cumulative_hist,
                   (double)(cumulative_hist)/box_volume,
                   sqrt((double)(cumulative_hist))/box_volume,
                   n_theory_1,n_theory_2);
         }
         fclose(fp_out);
         SID_log("Done.",SID_LOG_CLOSE);

         // Clean-up
         SID_free(SID_FARG data);
         SID_free(SID_FARG data_s);
         SID_free(SID_FARG bin);
         SID_free(SID_FARG bin_median);
         SID_free(SID_FARG hist);
         SID_free(SID_FARG properties_group);
         SID_free(SID_FARG properties_subgroup);
         if(flag_use_profiles){
            SID_free(SID_FARG profile_group);
            SID_free(SID_FARG profile_subgroup);
         }
         SID_free(SID_FARG group_SO_data);
         fclose(fp_halos);
         fclose_catalog(&fp_catalog_groups);
         fclose_catalog(&fp_catalog_subgroups);
         fclose_multifile(&fp_SO);
         SID_log("Done.",SID_LOG_CLOSE);
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  SID_free(SID_FARG z_list); 
  free_cosmo(&cosmo);

  SID_exit(ERROR_NONE);
}
