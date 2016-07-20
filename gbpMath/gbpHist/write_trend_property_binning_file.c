#include <stdio.h>
#include <gbpLib.h>
#include <gbpHist.h>

void write_trend_property_binning_file(trend_property_info *property,const char *filename_output_root){
   char  filename_out[MAX_FILENAME_LENGTH];
   sprintf(filename_out,"%s_%s_bins.txt",filename_output_root,property->name);
   SID_log("Writing binning description file to {%s}...",SID_LOG_OPEN,filename_out);
   if(SID.I_am_Master){
      FILE *fp_out=fopen(filename_out,"w");
      int i_column=1;
      fprintf(fp_out,"#  Ordinate %s-binning for {%s}\n",property->name,filename_output_root);
      fprintf(fp_out,"#  Column (%03d): Bin index\n",i_column++);
      fprintf(fp_out,"#         (%03d): %s - lo\n",  i_column++,property->name);
      fprintf(fp_out,"#         (%03d): %s - hi\n",  i_column++,property->name);
      hist_info *hist=property->hist;
      for(int i_bin=0;i_bin<hist->n_bins;i_bin++)
         fprintf(fp_out,"%03d %le %le\n",i_bin,histogram_bin_x_lo(hist,i_bin),histogram_bin_x_hi(hist,i_bin));
      fclose(fp_out);
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

