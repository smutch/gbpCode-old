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
   double  redshift;
   int     i_rank;
   int     i_grouping;
   int     i_grouping_start;
   int     i_grouping_stop;
   char    filename_in[MAX_FILENAME_LENGTH];
   char    filename_cosmology[MAX_FILENAME_LENGTH];
   char    filename_in_model[MAX_FILENAME_LENGTH];
   char    filename_out_root[MAX_FILENAME_LENGTH];
   char    grouping_name[6];
   char    filename_TF[256];
   char    n_string[64];
   int     n[3];
   double  box_size;
   double  L[3];
   size_t  n_all;
   FILE   *fp_in;
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
 
   // Initialization -- MPI etc.
   SID_init(&argc,&argv,NULL,NULL);
   if(argc!=12)
     SID_trap_error("Incorrect syntax.",ERROR_SYNTAX);

   // Parse arguments
   int n_jack;
   strcpy(filename_in,            argv[1]);
   strcpy(filename_out_root,      argv[2]);
   redshift         =(double)atof(argv[3]);
   strcpy(filename_cosmology,     argv[4]);
   n_jack           =(int)   atoi(argv[5]);
   box_size         =(double)atof(argv[6]);
   x_column         =(int)   atoi(argv[7]);
   y_column         =(int)   atoi(argv[8]);
   z_column         =(int)   atoi(argv[9]);
   vx_column        =(int)   atoi(argv[10]);
   vy_column        =(int)   atoi(argv[11]);
   vz_column        =(int)   atoi(argv[12]);

   // Sanity check
   int n_groupings;
   n_groupings=i_grouping_stop-i_grouping_start+1;
   if(n_groupings<1)
       SID_trap_error("No groupings have been selected (you chose start=%d, stop=%d).",ERROR_LOGIC,i_grouping_start,i_grouping_stop);
 
   SID_log("Producing correllation functions for ascii file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in);
 
   // Initialize the objects structure
   plist_info plist;
   init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);

   // Set the bin ranges and bin sizes we will use
   double r_min_l1D;
   double r_max_1D;
   double dr_1D;
   double r_min_2D;
   double r_max_2D;
   double dr_2D;
   r_min_l1D=  0.1;
   r_max_1D =200.0;
   dr_1D    = 10.0;
   r_min_2D =  0.0;
   r_max_2D = 50.0;
   dr_2D    =  2.0;

   // Set the factor by which the randoms will exceed the objects
   int F_random=5;

   // Set the PHK boundary width (must be a power of 2)
   int n_bits_PHK=5; 

   // Count the number of objects involved
   int   n_data;
   char  filename_count[MAX_FILENAME_LENGTH];
   FILE *fp_count;
   sprintf(filename_count,"%s",filename_in);
   fp_count=fopen(filename_count,"r");
   n_data=count_lines_data(fp_count);
   fclose(fp_count);
 
   // Initialize the power spectrum
   cfunc_info cfunc;
   init_cfunc(&cfunc,filename_cosmology,n_data,F_random*n_data,n_bits_PHK,
              redshift,box_size,n_jack,
              r_min_l1D,r_max_1D,dr_1D,
              r_min_2D, r_max_2D,dr_2D);

   // Loop over ithe real-space and 3 redshift-space frames
   int flag_compute_randoms=TRUE;
   int i_run;
   for(i_run=0;i_run<4;i_run++){
  
      // Read catalog
      switch(i_run){
      case 0:
         SID_log("Processing real-space ...",SID_LOG_OPEN|SID_LOG_TIMER);
         read_atable(filename_in,&plist,x_column,y_column,z_column,vx_column,vy_column,vz_column,"objs",READ_GROUPING_PHK,cfunc.n_bits_PHK,cfunc.PHK_width,box_size);
         break;
      case 1:
         SID_log("Processing v_x redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
         read_atable(filename_in,&plist,x_column,y_column,z_column,vx_column,vy_column,vz_column,"objs",READ_GROUPING_PHK|READ_GROUPING_ADD_VX,cfunc.n_bits_PHK,cfunc.PHK_width,box_size,redshift,cfunc.cosmo);
         break;
      case 2:
         SID_log("Processing v_y redshift space...",SID_LOG_OPEN|SID_LOG_TIMER);
         read_atable(filename_in,&plist,x_column,y_column,z_column,vx_column,vy_column,vz_column,"objs",READ_GROUPING_PHK|READ_GROUPING_ADD_VY,cfunc.n_bits_PHK,cfunc.PHK_width,box_size,redshift,cfunc.cosmo);
         break;
      case 3:
         SID_log("Processing v_z redsift space...",SID_LOG_OPEN|SID_LOG_TIMER);
         read_atable(filename_in,&plist,x_column,y_column,z_column,vx_column,vy_column,vz_column,"objs",READ_GROUPING_PHK|READ_GROUPING_ADD_VZ,cfunc.n_bits_PHK,cfunc.PHK_width,box_size,redshift,cfunc.cosmo);
         break;
      }

      // Generate randoms
      if(flag_compute_randoms){
         char filename_randoms[MAX_FILENAME_LENGTH];
         sprintf(filename_randoms,"%s_randoms.dat",filename_out_root);
         generate_randoms(&cfunc,&plist,"objs","randoms",filename_randoms);
         flag_compute_randoms=FALSE;
      }

      // Compute power spectrum
      compute_cfunc(&plist,"objs","randoms",&cfunc,i_run);

      // Now that all 4 runs are done, let's write the results
      write_cfunc(&cfunc,filename_out_root,&plist,"objs","randoms",i_run);

      SID_log("Done.",SID_LOG_CLOSE);
   } // Loop over 4 P(k)'s
 
   // Clean-up
   free_cfunc(&cfunc);
 
   SID_log("Done.",SID_LOG_CLOSE);
   SID_exit(ERROR_NONE);
}

