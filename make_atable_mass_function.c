#define  _MAIN
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpClustering.h>

int main(int argc, char *argv[]){
   double  redshift;
   char    filename_in[MAX_FILENAME_LENGTH];
   char    filename_out[MAX_FILENAME_LENGTH];
   double  box_size;
   double  M_min,M_max,M_med;
   char   *line=NULL;
   size_t  line_length=0;
   int     M_column;
   int     flag_log;
 
   // Initialization -- MPI etc.
   SID_init(&argc,&argv,NULL);
   if(argc!=7)
     SID_trap_error("Incorrect syntax.",ERROR_SYNTAX);

   // Parse arguments
   strcpy(filename_in, argv[1]);
   strcpy(filename_out,argv[2]);
   redshift=(double)atof(argv[3]);
   box_size=(double)atof(argv[4]);
   M_column=(int)   atoi(argv[5]);
   flag_log=(int)   atoi(argv[6]);

   SID_log("Producing a mass function for ascii file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in);
 
   // Initialize cosmology
   cosmo_info *cosmo;
   init_cosmo_std(&cosmo);
   double h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];

   // Open file
   FILE *fp_in;
   if((fp_in=fopen(filename_in,"r"))==NULL)
      SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);

   // Count number of lines in the input file 
   int n_data   =count_lines_data(fp_in);
   int n_bins   =100;
   
   // Allocate memory for the data. Read it and sort it in ascending order 
   double *data =(double *)malloc(sizeof(double)*n_data);
   for(int i=0;i<n_data;i++){
     grab_next_line_data(fp_in,&line,&line_length);
     grab_double(line,M_column,&(data[i]));
     if(!flag_log)
        data[i]=take_alog10(data[i]);
   }
   fclose(fp_in);

   // Perform sort
   merge_sort(data,n_data,NULL,SID_INT,SORT_INPLACE_ONLY,SORT_COMPUTE_INPLACE);

   // Compile histogram
   double *bin       =(double *)malloc(sizeof(double)*(n_bins+1));
   double *bin_median=(double *)malloc(sizeof(double)*n_bins);
   int    *hist      =(int    *)malloc(sizeof(int)*n_bins);
   int     i_data_lo=0;
   int     i_data_hi=0;
   int     n_target =n_data/n_bins+(!(n_data%n_bins!=0));
   for(int i_bin=0;i_bin<n_bins;i_bin++){
     if(i_data_hi<n_data-1){
       i_data_lo =i_data_hi;
       i_data_hi =MIN(i_data_lo+n_target,n_data);
       n_target =(n_data-i_data_lo)/(n_bins-i_bin)+(!((n_data-i_data_lo)%(n_bins-i_bin)!=0));
       while(data[i_data_hi]==data[i_data_hi-1] && i_data_hi<n_data-1) i_data_hi++;
       if(data[i_data_hi]==data[i_data_hi-1])                          i_data_hi++;
       int i_data_mid=(i_data_lo+i_data_hi)/2;
       hist[i_bin]=i_data_hi-i_data_lo;
       if(i_bin==0)
         bin[i_bin]=data[i_data_lo];
       bin[i_bin]  =data[i_data_lo];
       bin[i_bin+1]=data[i_data_hi];
       if((i_data_hi-i_data_lo+1)%2)
         bin_median[i_bin]=(data[i_data_mid]+data[i_data_mid+1])/2.;
       else
         bin_median[i_bin]=data[i_data_mid];
     }
     else
       n_bins=i_bin;
   }

   // Write mass function
   FILE *fp_out;
   if((fp_out=fopen(filename_out,"w"))==NULL){
     fprintf(stderr,"Error opening output file {%s}.\n",filename_out);
     free(data);
     return(1);
   }
   double box_volume=box_size*box_size*box_size;
   for(int i=0;i<n_bins;i++){
     double dlogM=take_log10(bin[i+1])-take_log10(bin[i]);
     double dn_dlogM_theory=mass_function(bin_median[i]*M_SOL/h_Hubble,
                                          redshift,
                                          &cosmo,
                                          PSPEC_LINEAR_TF,
                                          PSPEC_ALL_MATTER,
                                          MF_ST)*pow(M_PER_MPC/h_Hubble,3.0);
     fprintf(fp_out,"%11.4le %11.4le %11.4le %6d %11.4le %11.4le %10.4le\n",
             bin[i],
             bin_median[i],
             bin[i+1],
             hist[i],
             (double)(hist[i])/(box_volume*dlogM),
             sqrt((double)(hist[i]))/(box_volume*dlogM),
             dn_dlogM_theory);
   }
   fclose(fp_out);

   // Free allocated memory
   free(data);
   free(bin);
   free(bin_median);
   free(hist);
  
   SID_log("Done.",SID_LOG_CLOSE);
   SID_exit(ERROR_NONE);
}

