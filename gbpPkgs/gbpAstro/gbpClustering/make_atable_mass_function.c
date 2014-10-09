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
   char    filename_cosmology[MAX_FILENAME_LENGTH];
   double  box_size;
   double  lM_min,dlM;
   char   *line=NULL;
   size_t  line_length=0;
   int     M_column;
   int     n_bins;
   int     flag_log;
 
   // Initialization -- MPI etc.
   SID_init(&argc,&argv,NULL);
   if(argc!=11)
     SID_trap_error("Incorrect syntax.",ERROR_SYNTAX);

   // Parse arguments
   strcpy(filename_in, argv[1]);
   strcpy(filename_out,argv[2]);
   redshift=(double)atof(argv[3]);
   strcpy(filename_cosmology,argv[4]);
   box_size=(double)atof(argv[5]);
   M_column=(int)   atoi(argv[6]);
   flag_log=(int)   atoi(argv[7]);
   lM_min  =(double)atof(argv[8]);
   dlM     =(double)atof(argv[9]);
   n_bins  =(int)   atoi(argv[10]);

   SID_log("Producing a mass function for ascii file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in);
 
   // Initialize cosmology
   cosmo_info *cosmo;
   read_gbpCosmo_file(&cosmo,filename_cosmology);
   double h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];

   // Open file
   FILE *fp_in;
   if((fp_in=fopen(filename_in,"r"))==NULL)
      SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);

   // Allocate memory for the data. Read it and sort it in ascending order 
   SID_log("Reading data...",SID_LOG_OPEN|SID_LOG_TIMER);
   int     n_data=count_lines_data(fp_in);
   SID_log("(%d items)...",SID_LOG_CONTINUE,n_data);
   double *data  =(double *)malloc(sizeof(double)*n_data);
   for(int i=0;i<n_data;i++){
     grab_next_line_data(fp_in,&line,&line_length);
     grab_double(line,M_column,&(data[i]));
     if(!flag_log)
        data[i]=take_log10(data[i]);
   }
   fclose(fp_in);
   SID_free(SID_FARG line);
   SID_log("Done.",SID_LOG_CLOSE);

   // Perform sort
   SID_log("Sorting data...",SID_LOG_OPEN|SID_LOG_TIMER);
   merge_sort(data,n_data,NULL,SID_DOUBLE,SORT_INPLACE_ONLY,SORT_COMPUTE_INPLACE);
   SID_log("Done.",SID_LOG_CLOSE);

   // Compile histogram
   SID_log("Computing mass function...",SID_LOG_OPEN|SID_LOG_TIMER);
   int i_data = 0;
   while(data[i_data]<lM_min && i_data<(n_data-1)) i_data++;
   double *bin       =(double *)SID_malloc(sizeof(double)*(n_bins+1));
   double *bin_median=(double *)SID_malloc(sizeof(double)*n_bins);
   int    *hist      =(int    *)SID_calloc(sizeof(int)   *n_bins);
   double  lM_bin_min=lM_min;
   double  lM_bin_max=lM_min+dlM;
   int     i_data_lo=-1;
   int     i_data_hi=-1;
   int     i_bin;
   for(i_bin=0;i_bin<n_bins;i_bin++){
     bin[i_bin]=lM_bin_min;
     i_data_lo=i_data;
     while(data[i_data]<lM_bin_min && i_data<(n_data-1)){
        hist[i_bin]++;
        i_data++;
     }
     i_data_hi=i_data;
     int i_data_mid=(i_data_lo+i_data_hi)/2;
     if(hist[i_bin]>0){
       if(hist[i_bin]%2)
         bin_median[i_bin]=data[i_data_mid];
       else
         bin_median[i_bin]=0.5*(data[i_data_mid]+data[i_data_mid+1]);
     }
     else
        bin_median[i_bin]=0.5*(lM_bin_max+lM_bin_min);
     lM_bin_min=lM_bin_max;
     lM_bin_max=lM_min+((double)i_bin)*dlM;
   }
   bin[i_bin]=lM_bin_max;
   SID_log("Done.",SID_LOG_CLOSE);

   // Write mass function
   FILE *fp_out;
   if((fp_out=fopen(filename_out,"w"))==NULL){
     fprintf(stderr,"Error opening output file {%s}.\n",filename_out);
     SID_free(SID_FARG data);
     return(1);
   }
   SID_log("Writing results to {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_out);
   double box_volume=box_size*box_size*box_size;
   fprintf(fp_out,"# Mass function for column %d in {%s}\n",M_column,filename_in);
   fprintf(fp_out,"# Column (1): M_lo     [source units]\n");
   fprintf(fp_out,"#        (2): M_median [source units]\n");
   fprintf(fp_out,"#        (3): M_hi     [source units]\n");
   fprintf(fp_out,"#        (4): No. in bin\n");
   fprintf(fp_out,"#        (5): MFn (per unit volume, per dlogM)\n");
   fprintf(fp_out,"#        (6): +/- MFn\n");
   fprintf(fp_out,"#        (7): Sheth & Tormen MFn\n");
   for(int i=0;i<n_bins;i++){
     double dn_dlogM_theory=mass_function(bin_median[i]*M_SOL/h_Hubble,
                                          redshift,
                                          &cosmo,
                                          MF_ST)*pow(M_PER_MPC/h_Hubble,3.0);
     fprintf(fp_out,"%11.4le %11.4le %11.4le %6d %11.4le %11.4le %10.4le\n",
             bin[i],
             bin_median[i],
             bin[i+1],
             hist[i],
             (double)(hist[i])/(box_volume*dlM),
             sqrt((double)(hist[i]))/(box_volume*dlM),
             dn_dlogM_theory);
   }
   fclose(fp_out);
   SID_log("Done.",SID_LOG_CLOSE);

   // Free allocated memory
   SID_free(SID_FARG data);
   SID_free(SID_FARG bin);
   SID_free(SID_FARG bin_median);
   SID_free(SID_FARG hist);
  
   SID_log("Done.",SID_LOG_CLOSE);
   SID_exit(ERROR_NONE);
}

