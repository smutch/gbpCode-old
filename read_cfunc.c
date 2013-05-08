#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>

void read_cfunc(cfunc_info *cfunc,char *filename_in_root,int j_run){
   // Now that all 4 runs are done, let's write the results
   SID_log("Reading correlation functions {%s;run=%d}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in_root,j_run);

   if(j_run<0 || j_run>3)
      SID_trap_error("Invalid run iteration specified (%d).",ERROR_LOGIC,j_run);

   // Set output filenames
   char filename_in_1D[MAX_FILENAME_LENGTH];
   char filename_in_2D[MAX_FILENAME_LENGTH];
   sprintf(filename_in_1D,"%s_1D_cfunc.dat",filename_in_root);
   sprintf(filename_in_2D,"%s_2D_cfunc.dat",filename_in_root);

   // Read 1-D Results to an ascii file
   FILE *fp_in;
   char  grouping_name[MAX_FILENAME_LENGTH];
   strcpy(grouping_name,filename_in_root);
   strip_path(grouping_name);

   // Read 1-D Results to a binary file
   int i_jack;
   int i_run;
   if((fp_in=fopen(filename_in_1D,"r"))==NULL)
      SID_trap_error("Error opening {%s}.",ERROR_IO_OPEN,filename_in_1D);
   for(i_run=0;i_run<=j_run;i_run++){
      if(i_run==0){
         fread(&(cfunc->n_1D),     sizeof(int),   1,fp_in);
         fread(&(cfunc->n_jack),   sizeof(int),   1,fp_in);
         fread(&(cfunc->r_min_l1D),sizeof(double),1,fp_in);
         fread(&(cfunc->r_max_1D), sizeof(double),1,fp_in);
         fread(&(cfunc->n_data),   sizeof(int),   1,fp_in);
         fread(&(cfunc->n_random), sizeof(int),   1,fp_in);
      }
      if(i_run==j_run){
         fread(cfunc->CFUNC_1D[i_run], sizeof(double),cfunc->n_1D,fp_in);
         fread(cfunc->dCFUNC_1D[i_run],sizeof(double),cfunc->n_1D,fp_in);
         fread(cfunc->COVMTX_1D[i_run],sizeof(double),cfunc->n_1D*cfunc->n_1D,fp_in);
         for(i_jack=0;i_jack<=cfunc->n_jack_total;i_jack++){
            fread(cfunc->DD_1D[i_jack],sizeof(long long),cfunc->n_1D,fp_in);
            fread(cfunc->DR_1D[i_jack],sizeof(long long),cfunc->n_1D,fp_in);
            fread(cfunc->RR_1D[i_jack],sizeof(long long),cfunc->n_1D,fp_in);
         }
         fread(cfunc->CFUNC_l1D[i_run], sizeof(double),cfunc->n_1D,fp_in);
         fread(cfunc->dCFUNC_l1D[i_run],sizeof(double),cfunc->n_1D,fp_in);
         fread(cfunc->COVMTX_l1D[i_run],sizeof(double),cfunc->n_1D*cfunc->n_1D,fp_in);
         for(i_jack=0;i_jack<=cfunc->n_jack_total;i_jack++){
            fread(cfunc->DD_l1D[i_jack],sizeof(long long),cfunc->n_1D,fp_in);
            fread(cfunc->DR_l1D[i_jack],sizeof(long long),cfunc->n_1D,fp_in);
            fread(cfunc->RR_l1D[i_jack],sizeof(long long),cfunc->n_1D,fp_in);
         }
      }
      else{
         fseeko(fp_in,(off_t)(sizeof(double)*cfunc->n_1D),            SEEK_CUR);
         fseeko(fp_in,(off_t)(sizeof(double)*cfunc->n_1D),            SEEK_CUR);
         fseeko(fp_in,(off_t)(sizeof(double)*cfunc->n_1D*cfunc->n_1D),SEEK_CUR);
         for(i_jack=0;i_jack<=cfunc->n_jack_total;i_jack++){
            fseeko(fp_in,(off_t)(sizeof(long long)*cfunc->n_1D),SEEK_CUR);
            fseeko(fp_in,(off_t)(sizeof(long long)*cfunc->n_1D),SEEK_CUR);
            fseeko(fp_in,(off_t)(sizeof(long long)*cfunc->n_1D),SEEK_CUR);
         }
         fseeko(fp_in,(off_t)(sizeof(double)*cfunc->n_1D),            SEEK_CUR);
         fseeko(fp_in,(off_t)(sizeof(double)*cfunc->n_1D),            SEEK_CUR);
         fseeko(fp_in,(off_t)(sizeof(double)*cfunc->n_1D*cfunc->n_1D),SEEK_CUR);
         for(i_jack=0;i_jack<=cfunc->n_jack_total;i_jack++){
            fseeko(fp_in,(off_t)(sizeof(long long)*cfunc->n_1D),SEEK_CUR);
            fseeko(fp_in,(off_t)(sizeof(long long)*cfunc->n_1D),SEEK_CUR);
            fseeko(fp_in,(off_t)(sizeof(long long)*cfunc->n_1D),SEEK_CUR);
         }
      }
   }
   fclose(fp_in);

   // Read 2D correlation functions to a binary file
   if((fp_in=fopen(filename_in_2D,"r"))==NULL)
      SID_trap_error("Error opening {%s}.",ERROR_IO_OPEN,filename_in_2D);
   for(i_run=0;i_run<=j_run;i_run++){
      if(i_run==0){
         fread(&(cfunc->n_2D),     sizeof(int),   1,fp_in);
         fread(&(cfunc->n_jack),   sizeof(int),   1,fp_in);
         fread(&(cfunc->r_min_2D), sizeof(double),1,fp_in);
         fread(&(cfunc->r_max_2D), sizeof(double),1,fp_in);
         fread(&(cfunc->n_data),   sizeof(int),   1,fp_in);
         fread(&(cfunc->n_random), sizeof(int),   1,fp_in);
      }
      if(i_run==j_run){
         fread(cfunc->CFUNC_2D[i_run], sizeof(double),cfunc->n_2D_total,fp_in);
         fread(cfunc->dCFUNC_2D[i_run],sizeof(double),cfunc->n_2D_total,fp_in);
         fread(cfunc->COVMTX_2D[i_run],sizeof(double),cfunc->n_2D_total*cfunc->n_2D_total,fp_in);
         for(i_jack=0;i_jack<=cfunc->n_jack_total;i_jack++){
            fread(cfunc->DD_2D[i_jack],sizeof(long long),cfunc->n_2D_total,fp_in);
            fread(cfunc->DR_2D[i_jack],sizeof(long long),cfunc->n_2D_total,fp_in);
            fread(cfunc->RR_2D[i_jack],sizeof(long long),cfunc->n_2D_total,fp_in);
         }
      }
      else{
         fseeko(fp_in,(off_t)(sizeof(double)*cfunc->n_2D_total),                  SEEK_CUR);
         fseeko(fp_in,(off_t)(sizeof(double)*cfunc->n_2D_total),                  SEEK_CUR);
         fseeko(fp_in,(off_t)(sizeof(double)*cfunc->n_2D_total*cfunc->n_2D_total),SEEK_CUR);
         for(i_jack=0;i_jack<=cfunc->n_jack_total;i_jack++){
            fseeko(fp_in,(off_t)(sizeof(long long)*cfunc->n_2D_total),SEEK_CUR);
            fseeko(fp_in,(off_t)(sizeof(long long)*cfunc->n_2D_total),SEEK_CUR);
            fseeko(fp_in,(off_t)(sizeof(long long)*cfunc->n_2D_total),SEEK_CUR);
         }
      }
   }
   fclose(fp_in);
   SID_log("Done.",SID_LOG_CLOSE);
}

