#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>

void read_groupings(char *filename_root,int grouping_number,plist_info *plist,int mode,slab_info *slab,...){
   // Interpret variable arguments
   GBPREAL     box_size;
   double      redshift;
   cosmo_info *cosmo;
   double      h_Hubble;
   va_list     vargs;
   va_start(vargs,slab);

   // Interpret mode
   int flag_add_zspace_x;
   int flag_add_zspace_y;
   int flag_add_zspace_z;
   int n_z_dims;
   flag_add_zspace_x=check_mode_for_flag(mode,READ_GROUPING_ADD_VX);
   flag_add_zspace_y=check_mode_for_flag(mode,READ_GROUPING_ADD_VY);
   flag_add_zspace_z=check_mode_for_flag(mode,READ_GROUPING_ADD_VZ);
   n_z_dims=flag_add_zspace_x+flag_add_zspace_y+flag_add_zspace_z;
   if(n_z_dims>1)
      SID_trap_error("More then one redshift-space dimension (%d) has been specified.",ERROR_LOGIC,n_z_dims);
   else if(n_z_dims==1){
      box_size=(GBPREAL)va_arg(vargs,double);
      redshift=va_arg(vargs,double);
      cosmo   =va_arg(vargs,cosmo_info *);
      h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
   }
   else
      box_size=0.; // Should be unused in this case. Set to zero to catch bugs.

   // Set filename and open file
   char  filename_in[MAX_FILENAME_LENGTH];
   char  filename_name[MAX_FILENAME_LENGTH];
   FILE *fp_in;
   sprintf(filename_in,"%s_grouping_%03d.dat",filename_root,grouping_number);
   strcpy(filename_name,filename_in);
   strip_path(filename_name);
   SID_log("Reading grouping catalog {%s}...",SID_LOG_OPEN,filename_name);
   fp_in=fopen(filename_in,"r");

   // Set the columns of the file we're reading
   int x_column     =10;
   int y_column     =11;
   int z_column     =12;
   int vx_sub_column=14;
   int vy_sub_column=15;
   int vz_sub_column=16;
   int vx_FoF_column=17;
   int vy_FoF_column=18;
   int vz_FoF_column=19;
   int vx_sys_column=20;
   int vy_sys_column=21;
   int vz_sys_column=22;
   int vx_column    =vx_sub_column; // Used for z-space distortions
   int vy_column    =vy_sub_column; // Used for z-space distortions
   int vz_column    =vz_sub_column; // Used for z-space distortions

   // Count the number of halos that will be read (slab decomposed)
   size_t  n_halos;
   size_t  n_halos_allocate;
   int     i_halo; 
   char   *line=NULL; 
   size_t  line_length=0;
   GBPREAL x_in;
   GBPREAL y_in;
   GBPREAL z_in;
   GBPREAL vx_in;
   GBPREAL vy_in;
   GBPREAL vz_in;
   n_halos=(size_t)count_lines_data(fp_in);
   for(i_halo=0,n_halos_allocate=0;i_halo<n_halos;i_halo++){
      grab_next_line_data(fp_in,&line,&line_length);
      grab_real(line,x_column,&x_in);
      if(flag_add_zspace_x){
         grab_real(line,vx_column,&vx_in);
         x_in+=(GBPREAL)(1e3*h_Hubble*((double)vx_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
         force_periodic(&x_in,0.,box_size);
      }
      if(x_in>=slab->x_min_local && x_in<slab->x_max_local)
         n_halos_allocate++;
   }
   rewind(fp_in);

   // Allocate storage arrays
   GBPREAL *x_halos;
   GBPREAL *y_halos;
   GBPREAL *z_halos;
   GBPREAL *vx_halos_sub;
   GBPREAL *vy_halos_sub;
   GBPREAL *vz_halos_sub;
   GBPREAL *vx_halos_FoF;
   GBPREAL *vy_halos_FoF;
   GBPREAL *vz_halos_FoF;
   GBPREAL *vx_halos_sys;
   GBPREAL *vy_halos_sys;
   GBPREAL *vz_halos_sys;
   x_halos     =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
   y_halos     =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
   z_halos     =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
   vx_halos_sub=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
   vy_halos_sub=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
   vz_halos_sub=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
   vx_halos_FoF=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
   vy_halos_FoF=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
   vz_halos_FoF=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
   vx_halos_sys=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
   vy_halos_sys=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
   vz_halos_sys=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);

   // Set initial values for position bounds of the objects we will read
   GBPREAL x_min;
   GBPREAL x_max;
   GBPREAL y_min;
   GBPREAL y_max;
   GBPREAL z_min;
   GBPREAL z_max;
   GBPREAL vx_min;
   GBPREAL vx_max;
   GBPREAL vy_min;
   GBPREAL vy_max;
   GBPREAL vz_min;
   GBPREAL vz_max;
   x_min = 1e10;
   x_max =-1e10;
   y_min = 1e10;
   y_max =-1e10;
   z_min = 1e10;
   z_max =-1e10;
   vx_min= 1e10;
   vx_max=-1e10;
   vy_min= 1e10;
   vy_max=-1e10;
   vz_min= 1e10;
   vz_max=-1e10;

   // Read the halos one at a time, making decisions about where they will go
   size_t n_halos_local;
   for(i_halo=0,n_halos_local=0;i_halo<n_halos;i_halo++){
     grab_next_line_data(fp_in,&line,&line_length);
     grab_real(line,x_column, &x_in);
     if(flag_add_zspace_x){
        grab_real(line,vx_column,&vx_in);
        x_in+=(GBPREAL)(1e3*h_Hubble*((double)vx_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
        force_periodic(&x_in,0.,box_size);
     }
     if(x_in>=slab->x_min_local && x_in<slab->x_max_local){
        grab_real(line,y_column,&y_in);
        grab_real(line,z_column,&z_in);
        if(flag_add_zspace_y){
           grab_real(line,vy_column,&vy_in);
           y_in+=(GBPREAL)(1e3*h_Hubble*((double)vy_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
           force_periodic(&y_in,0.,box_size);
        }
        if(flag_add_zspace_z){
           grab_real(line,vz_column,&vz_in);
           z_in+=(GBPREAL)(1e3*h_Hubble*((double)vz_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
           force_periodic(&z_in,0.,box_size);
        }
        x_halos[n_halos_local] =(GBPREAL)x_in;
        y_halos[n_halos_local] =(GBPREAL)y_in;
        z_halos[n_halos_local] =(GBPREAL)z_in;
        grab_real(line,vx_sub_column,&vx_in);
        grab_real(line,vy_sys_column,&vy_in);
        grab_real(line,vz_sys_column,&vz_in);
        vx_halos_sub[n_halos_local]=(GBPREAL)vx_in;
        vy_halos_sub[n_halos_local]=(GBPREAL)vy_in;
        vz_halos_sub[n_halos_local]=(GBPREAL)vz_in;
        grab_real(line,vx_FoF_column,&vx_in);
        grab_real(line,vy_FoF_column,&vy_in);
        grab_real(line,vz_FoF_column,&vz_in);
        vx_halos_FoF[n_halos_local]=(GBPREAL)vx_in;
        vy_halos_FoF[n_halos_local]=(GBPREAL)vy_in;
        vz_halos_FoF[n_halos_local]=(GBPREAL)vz_in;
        grab_real(line,vx_sys_column,&vx_in);
        grab_real(line,vy_sys_column,&vy_in);
        grab_real(line,vz_sys_column,&vz_in);
        vx_halos_sys[n_halos_local]=(GBPREAL)vx_in;
        vy_halos_sys[n_halos_local]=(GBPREAL)vy_in;
        vz_halos_sys[n_halos_local]=(GBPREAL)vz_in;           
        // Keep track of ranges for a sanity check
        x_min=MIN(x_min,x_halos[n_halos_local]);
        x_max=MAX(x_max,x_halos[n_halos_local]);
        y_min=MIN(y_min,y_halos[n_halos_local]);
        y_max=MAX(y_max,y_halos[n_halos_local]);
        z_min=MIN(z_min,z_halos[n_halos_local]);
        z_max=MAX(z_max,z_halos[n_halos_local]);
        vx_min=MIN(vx_min,vx_halos_sub[n_halos_local]);
        vx_max=MAX(vx_max,vx_halos_sub[n_halos_local]);
        vy_min=MIN(vy_min,vy_halos_sub[n_halos_local]);
        vy_max=MAX(vy_max,vy_halos_sub[n_halos_local]);
        vz_min=MIN(vz_min,vz_halos_sub[n_halos_local]);
        vz_max=MAX(vz_max,vz_halos_sub[n_halos_local]);
        n_halos_local++;
     }
   }

   // Some sanity checks on what we have read
   size_t n_halos_read;
   SID_Allreduce(&n_halos_local,&n_halos_read,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&x_min, 1,SID_REAL,SID_MIN,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&x_max, 1,SID_REAL,SID_MAX,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&y_min, 1,SID_REAL,SID_MIN,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&y_max, 1,SID_REAL,SID_MAX,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&z_min, 1,SID_REAL,SID_MIN,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&z_max, 1,SID_REAL,SID_MAX,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&vx_min,1,SID_REAL,SID_MIN,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&vx_max,1,SID_REAL,SID_MAX,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&vy_min,1,SID_REAL,SID_MIN,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&vy_max,1,SID_REAL,SID_MAX,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&vz_min,1,SID_REAL,SID_MIN,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&vz_max,1,SID_REAL,SID_MAX,SID.COMM_WORLD);
   if(n_halos_read!=n_halos)
      SID_trap_error("The correct number of halos were read (ie. %lld!=%lld).  There must be a coordinate/box size problem.",
                     ERROR_LOGIC,n_halos_read,n_halos);
   SID_log("x_range =%le->%le",SID_LOG_COMMENT,x_min, x_max);
   SID_log("y_range =%le->%le",SID_LOG_COMMENT,y_min, y_max);
   SID_log("z_range =%le->%le",SID_LOG_COMMENT,z_min, z_max);
   SID_log("vx_range=%le->%le",SID_LOG_COMMENT,vx_min,vx_max);
   SID_log("vy_range=%le->%le",SID_LOG_COMMENT,vy_min,vy_max);
   SID_log("vz_range=%le->%le",SID_LOG_COMMENT,vz_min,vz_max);
   fclose(fp_in);

   // Store Arrays
   ADaPS_store(&(plist->data),(void *)&n_halos,      "n_all_halos", ADaPS_SCALAR_SIZE_T);
   ADaPS_store(&(plist->data),(void *)&n_halos_local,"n_halos",     ADaPS_SCALAR_SIZE_T);
   ADaPS_store(&(plist->data),(void *)x_halos,       "x_halos",     ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)y_halos,       "y_halos",     ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)z_halos,       "z_halos",     ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)vx_halos_sub,  "vx_halos",    ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)vy_halos_sub,  "vy_halos",    ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)vz_halos_sub,  "vz_halos",    ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)vx_halos_sys,  "vx_sys_halos",ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)vy_halos_sys,  "vy_sys_halos",ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)vz_halos_sys,  "vz_sys_halos",ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)vx_halos_FoF,  "vx_FoF_halos",ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)vy_halos_FoF,  "vy_FoF_halos",ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)vz_halos_FoF,  "vz_FoF_halos",ADaPS_DEFAULT);

   // Clean-up
   SID_free(SID_FARG line);
   va_end(vargs);
   SID_log("Done.",SID_LOG_CLOSE);
}

