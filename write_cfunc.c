#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>

void write_cfunc(cfunc_info *cfunc,char *filename_out_root,plist_info *plist,char *species_name,char *randoms_name){
   // Now that all 4 runs are done, let's write the results
   SID_log("Writing correlation functions...",SID_LOG_OPEN|SID_LOG_TIMER);
   if(SID.I_am_Master){
      // Set output filenames
      char filename_out_1D[MAX_FILENAME_LENGTH];
      char filename_out_2D[MAX_FILENAME_LENGTH];
      sprintf(filename_out_1D,"%s_1D_cfunc.dat",filename_out_root);
      sprintf(filename_out_2D,"%s_2D_cfunc.dat",filename_out_root);

      // Fetch the number of objects that contributed to the calculation
      size_t n_species;
      size_t n_randoms;
      n_species=((size_t *)ADaPS_fetch(plist->data,"n_all_%s",species_name))[0];
      n_randoms=((size_t *)ADaPS_fetch(plist->data,"n_all_%s",randoms_name))[0];

      // Write 1-D Results
      FILE *fp_out;
      int   i_r;
      char  grouping_name[MAX_FILENAME_LENGTH];
      strcpy(grouping_name,filename_out_root);
      strip_path(grouping_name);
      fp_out=fopen(filename_out_1D,"w");
      fprintf(fp_out,"# 1D Correlation functions for grouping {%s}\n",grouping_name);
      fprintf(fp_out,"#\n");
      fprintf(fp_out,"# No. of objects:         %zd\n",           n_species);
      fprintf(fp_out,"# No. of randoms:         %zd\n",           n_randoms);
      fprintf(fp_out,"# Redshift:               %5.3lf\n",        cfunc->redshift);
      fprintf(fp_out,"# Box size:               %9.3le [Mpc/h]\n",cfunc->box_size);
      fprintf(fp_out,"# No. of jack-knifes:     %d^3\n",          cfunc->n_jack);
      fprintf(fp_out,"#\n");
      fprintf(fp_out,"# Start of first log-bin: %9.3le [Mpc/h]\n",cfunc->r_min_l1D);
      fprintf(fp_out,"# Size of linear 1D bins: %9.3le [Mpc/h]\n",cfunc->dr_1D);
      fprintf(fp_out,"# Size of linear 2D bins: %9.3le [Mpc/h]\n",cfunc->dr_2D);
      int i_column=1;
      int i_run;
      fprintf(fp_out,"#\n");
      fprintf(fp_out,"# Column: (%02d) r_log [Mpc/h]; Logarithmically spaced bins\n",i_column++);
      fprintf(fp_out,"#         (%02d) r_lin [Mpc/h]; Linearly        spaced bins\n",i_column++);
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
         fprintf(fp_out,"#         (%02d) xi(r_log);  %s-space\n",i_column++,run_name);
         fprintf(fp_out,"#         (%02d) dxi(r_log); %s-space\n",i_column++,run_name);
         fprintf(fp_out,"#         (%02d) xi(r_lin);  %s-space\n",i_column++,run_name);
         fprintf(fp_out,"#         (%02d) dxi(r_lin); %s-space\n",i_column++,run_name);
      }
      fprintf(fp_out,"#\n");
      double r_log;
      double r_lin;
      r_log=take_log10(cfunc->r_min_l1D)+0.5*cfunc->dr_l1D;
      r_lin=0.5*cfunc->dr_1D;
      for(i_r=0;i_r<cfunc->n_1D;i_r++,r_log+=cfunc->dr_l1D,r_lin+=cfunc->dr_1D){
         fprintf(fp_out,"%9.4le %9.3le ",take_alog10(r_log),r_lin);
         for(i_run=0;i_run<4;i_run++)
            fprintf(fp_out,"  %le %le %le %le",
                           cfunc->CFUNC_l1D[i_run][i_r],
                           cfunc->dCFUNC_l1D[i_run][i_r],
                           cfunc->CFUNC_1D[i_run][i_r],
                           cfunc->dCFUNC_1D[i_run][i_r]);
         fprintf(fp_out,"\n");
      }
      fclose(fp_out);

      // Write 2D correlation functions
      /*
      fp_out=fopen(filename_out_2D,"w");
      fwrite(&(cfunc->n_r_2D),  sizeof(int),   1,fp_out);
      fwrite(&(cfunc->k_min_2D),sizeof(double),1,fp_out);
      fwrite(&(cfunc->k_max_2D),sizeof(double),1,fp_out);
      fwrite(&(cfunc->dk_2D),   sizeof(double),1,fp_out);
      fwrite(cfunc->n_modes_2D, sizeof(int),   cfunc->n_r_2D*cfunc->n_r_2D,fp_out);
      for(i_run=0;i_run<4;i_run++){
         fwrite(cfunc->P_k_2D[i_run], sizeof(double),cfunc->n_r_2D*cfunc->n_r_2D,fp_out);
         fwrite(cfunc->dP_k_2D[i_run],sizeof(double),cfunc->n_r_2D*cfunc->n_r_2D,fp_out);
      }
      fclose(fp_out);
      */
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

