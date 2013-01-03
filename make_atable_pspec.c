#define  _MAIN
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>

int main(int argc, char *argv[]){
   int     i,j,i_x,i_y,i_z;
   int     n_2D_total;
   char    species_name[256];
   double  h_Hubble;
   double  Omega_Lambda;
   double  Omega_M;
   double  Omega_b;
   double  f_gas;
   double  Omega_k;
   double  sigma_8;
   double  n_spec;
   double  redshift;
   int     i_rank;
   char    filename_in[MAX_FILENAME_LENGTH];
   char    filename_in_model[MAX_FILENAME_LENGTH];
   char    filename_out_root[MAX_FILENAME_LENGTH];
   char    grouping_name[6];
   char    filename_TF[256];
   char    n_string[64];
   int             n[3];
   double          x_in,y_in,z_in,vx_in,vy_in,vz_in;
   double          box_size;
   double          L[3];
   size_t          n_all;
   FILE           *fp_in;
   double  M_min,M_max,M_med;
   double  V_min,V_max,V_med;
   char   *line=NULL;
   size_t  line_length=0;
   int     i_bin,j_bin;
   int     i_file,n_files;
   int     n_1D,n_2D;
   int     grid_size;
   int     x_column;
   int     y_column;
   int     z_column;
   int     vx_column;
   int     vy_column;
   int     vz_column;
   double  k_min_1D;
   double  dk_1D;
   double  k_max_1D;
   double  k_min_2D;
   double  dk_2D;
   double  k_max_2D;
 
   // Initialization -- MPI etc.
   SID_init(&argc,&argv,NULL);
   if(argc!=12)
     SID_trap_error("Incorrect syntax.",ERROR_SYNTAX);

   // Parse arguments
   strcpy(filename_in, argv[1]);
   strcpy(filename_out_root,argv[2]);
   redshift         =(double)atof(argv[3]);
   grid_size        =(int)   atoi(argv[4]);
   box_size         =(double)atof(argv[5]);
   x_column         =(int)   atoi(argv[6]);
   y_column         =(int)   atoi(argv[7]);
   z_column         =(int)   atoi(argv[8]);
   vx_column        =(int)   atoi(argv[9]);
   vy_column        =(int)   atoi(argv[10]);
   vz_column        =(int)   atoi(argv[11]);

   SID_log("Producing power spectra for ascii file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in);
 
   // Set the k ranges
   double k_Nyq;
   k_Nyq   =(TWO_PI*(double)grid_size/box_size)/2.;
   k_min_1D=0.02;
   dk_1D   =0.02;
   k_max_1D=dk_1D*(float)((int)(k_Nyq/dk_1D));
   k_min_2D=0.0;
   dk_2D   =0.02;
   k_max_2D=dk_2D*(float)((int)(k_Nyq/dk_2D));

   // Initialize the objects structure
   plist_info plist;
   init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
 
   // Initialize the power spectrum
   pspec_info pspec;
   init_pspec(&pspec,
 	     MAP2GRID_DIST_DWT20,
             redshift,box_size,grid_size,
             k_min_1D,k_max_1D,dk_1D,
             k_min_2D,k_max_2D,dk_2D);
 
  
   // Loop over ithe real-space and 3 redshift-space frames
   int i_run;
   for(i_run=0;i_run<4;i_run++){
  
      // Read catalog
      switch(i_run){
      case 0:
         SID_log("Processing real-space ...",SID_LOG_OPEN|SID_LOG_TIMER);
         read_atable(filename_in,&plist,x_column,y_column,z_column,vx_column,vy_column,vz_column,"objs",READ_GROUPING_DEFAULT,&(pspec.FFT.slab),box_size);
         break;
      case 1:
         SID_log("Processing v_x redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
         read_atable(filename_in,&plist,x_column,y_column,z_column,vx_column,vy_column,vz_column,"objs",READ_GROUPING_DEFAULT|READ_GROUPING_ADD_VX,&(pspec.FFT.slab),box_size,redshift,pspec.cosmo);
         break;
      case 2:
         SID_log("Processing v_y redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
         read_atable(filename_in,&plist,x_column,y_column,z_column,vx_column,vy_column,vz_column,"objs",READ_GROUPING_DEFAULT|READ_GROUPING_ADD_VY,&(pspec.FFT.slab),box_size,redshift,pspec.cosmo);
         break;
      case 3:
         SID_log("Processing v_z redsift space...",SID_LOG_OPEN|SID_LOG_TIMER);
         read_atable(filename_in,&plist,x_column,y_column,z_column,vx_column,vy_column,vz_column,"objs",READ_GROUPING_DEFAULT|READ_GROUPING_ADD_VZ,&(pspec.FFT.slab),box_size,redshift,pspec.cosmo);
         break;
      }
  
      // Compute power spectrum
      compute_pspec(&plist,"objs",&pspec,i_run);
  
      SID_log("Done.",SID_LOG_CLOSE);
   } // Loop over 4 P(k)'s
  
   // Now that all 4 runs are done, let's write the results
   char filename_out_root_grouping[MAX_FILENAME_LENGTH];
   write_pspec(&pspec,filename_out_root,&plist,"objs");
 
   // Clean-up
   free_pspec(&pspec);
 
   SID_log("Done.",SID_LOG_CLOSE);
   SID_exit(ERROR_NONE);
}

