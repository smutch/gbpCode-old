#include <stdio.h>
#include <gbpLib.h>
#include <gbpHist.h>

void write_trend_ascii(trend_info *trend,const char *filename_output_root){

   // Set filename and open file
   char filename[MAX_FILENAME_LENGTH];
   sprintf(filename,"%s_%s.txt",filename_output_root,trend->ordinate->name);
   FILE *fp_out;
   SID_log("Writing trend to {%s}...",SID_LOG_OPEN,filename);
   SID_Barrier(SID.COMM_WORLD);
   if(SID.I_am_Master){
      fp_out=fopen(filename,"w");

      // Write header
      int i_column=1;
      fprintf(fp_out,"#  Column (%03d): Snapshot\n",   i_column++);
      fprintf(fp_out,"#         (%03d): %s bin - lo\n",i_column++,trend->ordinate->name);
      fprintf(fp_out,"#         (%03d): %s bin - hi\n",i_column++,trend->ordinate->name);
      fprintf(fp_out,"#         (%03d): n_halos_all\n",i_column++);
      // Loop over the properties
      trend_property_info *current_property=trend->coordinate_first;
      while(current_property!=NULL){
         fprintf(fp_out,"#         (%03d): n_halos_hist (%s)\n",      i_column++,current_property->name);
         fprintf(fp_out,"#         (%03d): %s\n",                     i_column++,current_property->name);
         fprintf(fp_out,"#         (%03d): %s - 68%% confidence lo\n",i_column++,current_property->name);
         fprintf(fp_out,"#         (%03d): %s - 68%% confidence hi\n",i_column++,current_property->name);
         fprintf(fp_out,"#         (%03d): %s - 95%% confidence lo\n",i_column++,current_property->name);
         fprintf(fp_out,"#         (%03d): %s - 95%% confidence hi\n",i_column++,current_property->name);
         current_property=current_property->next;
      }
   }

   // Finalize the snapshot histograms and write results
   hist_info *hist_ordinate=trend->ordinate->hist;
   finalize_histogram(hist_ordinate);
   for(int i_bin=0;i_bin<hist_ordinate->n_bins;i_bin++){
      if(SID.I_am_Master)
         fprintf(fp_out,"%3d %le %le %d",
                        i_bin,
                        histogram_bin_x_lo(hist_ordinate,i_bin),
                        histogram_bin_x_hi(hist_ordinate,i_bin),
                        hist_ordinate->bin_count[i_bin]);

      // Loop over the properties
      trend_property_info *current_property=trend->coordinate_first;
      while(current_property!=NULL){
         double x_peak;
         double x_68_lo;
         double x_68_hi;
         double x_95_lo;
         double x_95_hi;
         hist_info *hist_i=&(current_property->hist[i_bin]);
         finalize_histogram(hist_i);
         compute_histogram_range(hist_i,68.,GBP_HISTOGRAM_RANGE_HIST,&x_peak,&x_68_lo,&x_68_hi);
         compute_histogram_range(hist_i,95.,GBP_HISTOGRAM_RANGE_HIST,&x_peak,&x_95_lo,&x_95_hi);
         if(SID.I_am_Master)
            fprintf(fp_out," %d %le %le %le %le %le",
                           hist_i->count_hist,
                           x_peak,
                           x_68_lo,
                           x_68_hi,
                           x_95_lo,
                           x_95_hi);
         current_property=current_property->next;
      }
      if(SID.I_am_Master) fprintf(fp_out,"\n");
   }
   if(SID.I_am_Master) fclose(fp_out);

   SID_log("Done.",SID_LOG_CLOSE);
}

