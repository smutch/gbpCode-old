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
  float   lambda,v_c;
  float   r_min,r_max;
  double  box_size;
  double  subvolume_radius;
  double  M_cut;
  int     n_regions;

  SID_init(&argc,&argv,NULL,NULL);

  strcpy(filename_in_root,argv[1]);
  box_size          =atof(argv[2]);
  snap_number       =atoi(argv[3]);
  M_cut             =atof(argv[4]);
  subvolume_radius  =atof(argv[5]);
  n_regions         =atoi(argv[6]);

  int i_type=0;

  sprintf(filename_in_base,"%s",filename_in_root);
  strip_path(filename_in_base);

  if(SID.I_am_Master){
     int catalog_mode;
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
     SID_log("Processing %sgroup catalogs for snap %d...",SID_LOG_OPEN|SID_LOG_TIMER,prefix_text,snap_number);

     SID_log("Reading halos...",SID_LOG_OPEN|SID_LOG_TIMER);
     // Open catalog
     fp_catalog_info fp_catalog;
     fopen_catalog(filename_in_root,
                   snap_number,
                   catalog_mode,
                   &fp_catalog);
     int n_halos_total=fp_catalog.n_halos_total;

     // Allocate for arrays
     double *M_vir =(double *)SID_malloc(sizeof(double)*n_halos_total);
     double *x_halo=(double *)SID_malloc(sizeof(double)*n_halos_total);
     double *y_halo=(double *)SID_malloc(sizeof(double)*n_halos_total);
     double *z_halo=(double *)SID_malloc(sizeof(double)*n_halos_total);

     // Read halos
     halo_properties_info *properties=NULL;
     properties = (halo_properties_info *)SID_malloc(sizeof(halo_properties_info));
     int n_keep=0;
     for(int i_halo=0;i_halo<n_halos_total;i_halo++){
        // Read halo
        fread_catalog_file(&fp_catalog,NULL,properties,NULL,i_halo);

        // Store needed quantities
        if(TRUE){
           M_vir [n_keep]=(double)(properties->M_vir);
           x_halo[n_keep]=(double)(properties->position_COM[0]);
           y_halo[n_keep]=(double)(properties->position_COM[1]);
           z_halo[n_keep]=(double)(properties->position_COM[2]);
           n_keep++;
        }
     }

     // Close catalog
     SID_free(SID_FARG properties);
     fclose_catalog(&fp_catalog);
     SID_log("Done.",SID_LOG_CLOSE);

     // Perform analysis
     RNG_info RNG;
     double   subvolume_radius2=subvolume_radius*subvolume_radius;
     int      seed             =1926743;
     init_RNG(&seed,&RNG,RNG_DEFAULT);
     int n_occupied=0;
     for(int i_region=0;i_region<n_regions;i_region++){
        double x_cen=((double)random_number(&RNG))*box_size;
        double y_cen=((double)random_number(&RNG))*box_size;
        double z_cen=((double)random_number(&RNG))*box_size;
        double M_max=0.;
        double M_tot=0.;
        double x_max=0.;
        double y_max=0.;
        double z_max=0.;
        int    i_max=-1;
        int    N_subvolume=0;
        for(int i_halo=0;i_halo<n_keep;i_halo++){
           double dx     =d_periodic(x_halo[i_halo]-x_cen,box_size);
           double dy     =d_periodic(y_halo[i_halo]-y_cen,box_size);
           double dz     =d_periodic(z_halo[i_halo]-z_cen,box_size);
           double radius2=dx*dx+dy*dy+dz*dz;
           double radius =sqrt(dx*dx+dy*dy+dz*dz);
           if(radius2<subvolume_radius2){
              if(M_vir[i_halo]>=M_max && radius<0.2*subvolume_radius){
                 x_max=x_halo[i_halo];
                 y_max=y_halo[i_halo];
                 z_max=z_halo[i_halo];
                 M_max=M_vir[i_halo];
                 i_max=i_halo;
              }
              if(M_vir[i_halo]>=M_cut){
                 M_tot+=M_vir[i_halo];
                 N_subvolume++;
              }
           }
        }
        if(N_subvolume>0)
           n_occupied++;
        else
           printf("%6.2lf %6.2lf %6.2lf %le %le %6.2lf %6.2lf %6.2lf %d %d\n",x_cen,y_cen,z_cen,M_max,M_tot,x_max,y_max,z_max,i_max,N_subvolume);
     }

     // Clean-up
     free_RNG(&RNG);
     SID_free(SID_FARG M_vir );
     SID_free(SID_FARG x_halo);
     SID_free(SID_FARG y_halo);
     SID_free(SID_FARG z_halo);

     SID_log("Done.",SID_LOG_CLOSE);
  }  

  SID_exit(ERROR_NONE);
}
