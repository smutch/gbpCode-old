#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHalos.h>

typedef struct data_info_local data_in_info_local;
struct data_info_local{
   halo_properties_info group;   
   halo_properties_info subgroup;
   int n_particles_group;
   int n_particles_central;
   int n_particles_substructure;
};

void init_trend_f_local(trend_property_info *property,void *params,int i_hist,int *mode,gbp_va_list *vargs);
void init_trend_f_local(trend_property_info *property,void *params,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_FIXED;
   double x_min =0.;
   double dx    =0.01;
   int    n_x   =100;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double),&x_min);
   gbp_add_va_arg(vargs,sizeof(double),&dx);
   gbp_add_va_arg(vargs,sizeof(int),   &n_x);
}
void free_trend_f_local(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs);
void free_trend_f_local(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){}

void init_trend_from_hist_local(trend_property_info *property,void *parent_hist,int i_hist,int *mode,gbp_va_list *vargs);
void init_trend_from_hist_local(trend_property_info *property,void *parent_hist,int i_hist,int *mode,gbp_va_list *vargs){
   (*mode)=GBP_HISTOGRAM_FIXED;
   double x_min =((hist_info *)parent_hist)->x_min;
   double dx    =((hist_info *)parent_hist)->dx;
   int    n_x   =((hist_info *)parent_hist)->n_bins;
   gbp_va_start(vargs);
   gbp_add_va_arg(vargs,sizeof(double),&x_min);
   gbp_add_va_arg(vargs,sizeof(double),&dx);
   gbp_add_va_arg(vargs,sizeof(int),   &n_x);
}

void free_trend_from_hist_local(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs);
void free_trend_from_hist_local(trend_property_info *property,void *trees_in,int i_hist,int *mode,gbp_va_list *vargs){}

int calc_trend_M_group(trend_property_info *property,hist_info *hist,void *data_in);
int calc_trend_M_group(trend_property_info *property,hist_info *hist,void *data_in){
   return(calc_histogram_index(hist,take_log10(((data_info_local *)data_in)->group.M_vir)));
}

int calc_trend_f_cen(trend_property_info *property,hist_info *hist,void *data_in);
int calc_trend_f_cen(trend_property_info *property,hist_info *hist,void *data_in){
   int n_particles_central=((data_info_local *)data_in)->n_particles_central;
   int n_particles_group  =((data_info_local *)data_in)->n_particles_group;
   return(calc_histogram_index(hist,(double)(n_particles_central)/(double)(n_particles_group)));
}

int calc_trend_f_sub(trend_property_info *property,hist_info *hist,void *data_in);
int calc_trend_f_sub(trend_property_info *property,hist_info *hist,void *data_in){
   int n_particles_substructure=((data_info_local *)data_in)->n_particles_substructure;
   int n_particles_group       =((data_info_local *)data_in)->n_particles_group;
   return(calc_histogram_index(hist,(double)(n_particles_substructure)/(double)(n_particles_group)));
}

int main(int argc, char *argv[]){
  char    filename_SSimPL_root[256];
  char    filename_in_root[256];
  char    filename_halos_root[256];
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

  strcpy(filename_SSimPL_root,argv[1]);
  strcpy(filename_halos_root, argv[2]);
  box_size              =atof(argv[3]);
  snap_number_start     =atoi(argv[4]);
  snap_number_stop      =atoi(argv[5]);
  snap_number_step      =atoi(argv[6]);
  strcpy(filename_out_root,   argv[7]);
  sprintf(filename_in_root,"%s/catalogs/%s",filename_SSimPL_root,filename_halos_root);

  int flag_use_profiles=FALSE;

  // Create histograms
  hist_info   hist_logM_vir_group;
  hist_info   hist_logM_vir_subgroup;
  trend_info *trend_M_group;
  init_histogram(&hist_logM_vir_group,   GBP_HISTOGRAM_FIXED,7.,0.05,180);
  init_histogram(&hist_logM_vir_subgroup,GBP_HISTOGRAM_FIXED,7.,0.05,180);
  init_trend(&trend_M_group,"M_group",&hist_logM_vir_group,init_trend_from_hist_local,free_trend_from_hist_local,calc_trend_M_group);
  init_trend_coordinate(trend_M_group,"f_cen",NULL,init_trend_f_local,free_trend_f_local,calc_trend_f_cen);
  init_trend_coordinate(trend_M_group,"f_sub",NULL,init_trend_f_local,free_trend_f_local,calc_trend_f_sub);

  SID_log("Processing catalogs for snaps %d->%d...",SID_LOG_OPEN|SID_LOG_TIMER,snap_number_start,snap_number_stop);

  // Loop over the snapshots
  for(snap_number=snap_number_start;snap_number<=snap_number_stop;snap_number++){
     SID_log("Processing snapshot %03d...",SID_LOG_OPEN|SID_LOG_TIMER,snap_number);
     fp_catalog_info fp_catalog_groups;
     fp_catalog_info fp_catalog_subgroups;
     fopen_catalog(filename_in_root,
                   snap_number,
                   READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES,
                   &fp_catalog_groups);
     fopen_catalog(filename_in_root,
                   snap_number,
                   READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES,
                   &fp_catalog_subgroups);

     // Open the halo file (needed to get n_substructure)
     int n_groups_all,n_bytes_groups;
     char filename_groups[MAX_FILENAME_LENGTH];
     sprintf(filename_groups,"%s/halos/%s_%03d.catalog_groups",filename_SSimPL_root,filename_halos_root,snap_number);
     FILE *fp_groups=fopen(filename_groups,"r");
     fread(&n_groups_all,  sizeof(int),1,fp_groups);
     fread(&n_bytes_groups,sizeof(int),1,fp_groups);

     // Sanity checks
     if(n_bytes_groups!=4 && n_bytes_groups!=8)
        SID_trap_error("Invalid group offset byte length (%d).",ERROR_LOGIC,n_bytes_groups);

     // Skip to the number of subgroups
     fseeko(fp_groups,sizeof(int)   *n_groups_all,SEEK_CUR);
     fseeko(fp_groups,n_bytes_groups*n_groups_all,SEEK_CUR);

     // Initialize histograms
     clear_histogram(&hist_logM_vir_group);
     clear_histogram(&hist_logM_vir_subgroup);
     clear_trend(&trend_M_group);

     // Loop over groups
     data_info_local data;
     int i_group   =0;
     int i_subgroup=0;
     for(;i_group<fp_catalog_groups.n_halos_total;i_group++){
        // Read group from catalog
        fread_catalog_file(&fp_catalog_groups,NULL,&(data.group),NULL,i_group);

        // Read the number of subgroups for this group
        int n_subgroups;
        fread(&n_subgroups,sizeof(int),1,fp_groups);

        // Loop over subgroups
        data.n_particles_group       =data.group.n_particles;
        data.n_particles_central     =0;
        data.n_particles_substructure=0;
        for(int j_subgroup=0;j_subgroup<n_subgroups;j_subgroup++,i_subgroup++){
           fread_catalog_file(&fp_catalog_subgroups,NULL,&(data.subgroup),NULL,i_subgroup);

           if(data.subgroup.n_particles>data.n_particles_central)
              data.n_particles_central=data.subgroup.n_particles;
           data.n_particles_substructure+=data.subgroup.n_particles;

           // Populate the subgroup histograms/trends
           add_to_histogram(&hist_logM_vir_subgroup,take_log10(data.subgroup.M_vir));
        }

        // Populate the group histograms/trends
        data.n_particles_substructure-=data.n_particles_central;
        add_to_histogram(&hist_logM_vir_group,take_log10(data.group.M_vir));
        add_item_to_trend(trend_M_group,GBP_ADD_ITEM_TO_TREND_DEFAULT,&data);

     } // loop over groups

     // Finalize histograms and trends
     finalize_histogram(&hist_logM_vir_group);
     finalize_histogram(&hist_logM_vir_subgroup);
     finalize_trend    (trend_M_group);

     // Write results
     char filename_logM_vir_group[MAX_FILENAME_LENGTH];
     char filename_logM_vir_subgroup[MAX_FILENAME_LENGTH];
     char filename_f_cen[MAX_FILENAME_LENGTH];
     char filename_f_sub[MAX_FILENAME_LENGTH];
     sprintf(filename_logM_vir_group,   "%03d_logM_vir_group.txt",   snap_number);
     sprintf(filename_logM_vir_subgroup,"%03d_logM_vir_subgroup.txt",snap_number);
     sprintf(filename_f_cen,            "%s_%03d",filename_out_root,snap_number);
     sprintf(filename_f_sub,            "%s_%03d",filename_out_root,snap_number);
     //write_hist_ascii (hist_logM_vir_group   ,filename_logM_vir_group);
     //write_hist_ascii (hist_logM_vir_subgroup,filename_logM_vir_subgroup);
     write_trend_ascii(trend_M_group,           filename_f_cen);

     // Clean-up
     fclose_catalog(&fp_catalog_groups);
     fclose_catalog(&fp_catalog_subgroups);
     fclose(fp_groups);
     SID_log("Done.",SID_LOG_CLOSE);
  } // loop over snapshots

  // Clean-up
  free_histogram(&hist_logM_vir_group);
  free_histogram(&hist_logM_vir_subgroup);
  free_trend(&trend_M_group);

  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}
