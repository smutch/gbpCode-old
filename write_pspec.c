#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>

void write_pspec(pspec_info *pspec,char *filename_out_root,plist_info *plist,char *species_name){
   // Now that all 4 runs are done, let's write the results
   SID_log("Writing power spectra...",SID_LOG_OPEN|SID_LOG_TIMER);
   if(SID.I_am_Master){
      // Set output filenames
      char filename_out_1D[MAX_FILENAME_LENGTH];
      char filename_out_2D[MAX_FILENAME_LENGTH];
      sprintf(filename_out_1D,"%s_1D_pspec.dat",filename_out_root);
      sprintf(filename_out_2D,"%s_2D_pspec.dat",filename_out_root);

      // Fetch the number of objects that contributed to the calculation
      size_t n_species;
      n_species=((size_t *)ADaPS_fetch(plist->data,"n_all_%s",species_name))[0];

      // Write 1-D Results
      FILE *fp_out;
      int   i_k;
      char  grouping_name[MAX_FILENAME_LENGTH];
      strcpy(grouping_name,filename_out_root);
      strip_path(grouping_name);
      fp_out=fopen(filename_out_1D,"w");
      fprintf(fp_out,"# 1D Power spectra for dataset {%s}\n",grouping_name);
      fprintf(fp_out,"#\n");
      switch(pspec->mass_assignment_scheme){
         case MAP2GRID_DIST_DWT20:
            fprintf(fp_out,"# Mass assignment scheme: DWT20\n");
            break;
         case MAP2GRID_DIST_DWT12:
            fprintf(fp_out,"# Mass assignment scheme: DWT12\n");
            break;
         case MAP2GRID_DIST_NGP:
            fprintf(fp_out,"# Mass assignment scheme: NGP\n");
            break;
         case MAP2GRID_DIST_TSC:
            fprintf(fp_out,"# Mass assignment scheme: TSC\n");
            break;
      }
      fprintf(fp_out,"# No. of objects:         %zd\n",           n_species);
      fprintf(fp_out,"# Redshift:               %5.3lf\n",        pspec->redshift);
      fprintf(fp_out,"# Box size:               %9.3le [Mpc/h]\n",pspec->box_size);
      fprintf(fp_out,"# Grid size:              %d^3\n",          pspec->grid_size);
      int i_column=1;
      int i_run;
      fprintf(fp_out,"#\n");
      fprintf(fp_out,"# Column: (%02d) k_minimum [h/Mpc]\n",i_column++);
      fprintf(fp_out,"#         (%02d) k_average [h/Mpc]\n",i_column++);
      fprintf(fp_out,"#         (%02d) k_maximum [h/Mpc]\n",i_column++);
      fprintf(fp_out,"#         (%02d) # of Fourier modes in this bin\n",i_column++);
      for(i_run=0;i_run<4;i_run++){
         char run_name[128];
         switch(i_run){
            case 0:
               sprintf(run_name,"real");
               break;
            case 1:
               sprintf(run_name,"x-projected z");
               break;
            case 2:
               sprintf(run_name,"y-projected z");
               break;
            case 3:
               sprintf(run_name,"z-projected z");
               break;
         }
         if(pspec->flag_processed[i_run]){
            fprintf(fp_out,"#         (%02d) P(k)  [(Mpc/h)^3];  %s-space\n",i_column++,run_name);
            fprintf(fp_out,"#         (%02d) dP(k) [(Mpc/h)^3];  %s-space\n",i_column++,run_name);
         }
      }
      fprintf(fp_out,"#\n");
      double k_min;
      double k_max;
      k_min=pspec->k_min_1D;
      k_max=k_min+pspec->dk_1D;
      for(i_k=0;i_k<pspec->n_k_1D;i_k++){
         if(pspec->n_modes_1D[i_k]>0){
            if(i_k==(pspec->n_k_1D-1))
               k_max=pspec->k_max_1D;
            fprintf(fp_out,"%9.3le %9.3le %9.3le %8d ",k_min,pspec->k_1D[i_k],k_max,pspec->n_modes_1D[i_k]);
            for(i_run=0;i_run<4;i_run++){
               if(pspec->flag_processed[i_run])
                  fprintf(fp_out,"  %le %le",pspec->P_k_1D[i_run][i_k],pspec->dP_k_1D[i_run][i_k]);
            }
            fprintf(fp_out,"\n");
            k_min=k_max;
         }
         k_max+=pspec->dk_1D;
      }
      fclose(fp_out);

      // Write 2D power spectra
      fp_out=fopen(filename_out_2D,"w");
      fwrite(&(pspec->n_k_2D),  sizeof(int),   1,fp_out);
      fwrite(&(pspec->k_min_2D),sizeof(double),1,fp_out);
      fwrite(&(pspec->k_max_2D),sizeof(double),1,fp_out);
      fwrite(&(pspec->dk_2D),   sizeof(double),1,fp_out);
      fwrite(pspec->n_modes_2D, sizeof(int),   pspec->n_k_2D*pspec->n_k_2D,fp_out);
      for(i_run=0;i_run<4;i_run++){
         fwrite(&(pspec->flag_processed),sizeof(double),pspec->n_k_2D*pspec->n_k_2D,fp_out);
         if(pspec->flag_processed){
            fwrite(pspec->P_k_2D[i_run],    sizeof(double),pspec->n_k_2D*pspec->n_k_2D,fp_out);
            fwrite(pspec->dP_k_2D[i_run],   sizeof(double),pspec->n_k_2D*pspec->n_k_2D,fp_out);
         }
      }
      fclose(fp_out);
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

