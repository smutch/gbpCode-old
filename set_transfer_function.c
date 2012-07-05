#include <gbpLib.h>
#include <gbpRender.h>

int set_transfer_function(char *line,int i_word,interp_info **return_interp){
   int n_transfer;
   n_transfer=count_words(line)-i_word+1;
   if(n_transfer>2){
     double *transfer_array_x;
     double *transfer_array_y;
     n_transfer+=2; // We need to add low/hi interpolation anchors
     transfer_array_x=(double *)SID_malloc(sizeof(double)*n_transfer);
     transfer_array_y=(double *)SID_malloc(sizeof(double)*n_transfer);
     int  j_word;
     char temp_word[1024];
     for(j_word=1;j_word<n_transfer-1;j_word++){
       grab_word(line,i_word++,temp_word);
       search_and_replace(temp_word,","," ");
       if(count_words(temp_word)!=2)
         SID_trap_error("Error in formatting of transfer array {%s}{%s}",ERROR_LOGIC,line,temp_word);
       grab_double(temp_word,1,&(transfer_array_x[j_word]));
       grab_double(temp_word,2,&(transfer_array_y[j_word]));
     }
     // Create low/hi interpolation anchors
     transfer_array_y[0]           =transfer_array_y[1];
     transfer_array_y[n_transfer-1]=transfer_array_y[n_transfer-2];
     if(transfer_array_x[1]<0.)
       transfer_array_x[0]=10.*transfer_array_x[1];
     else
       transfer_array_x[0]=0.1*transfer_array_x[1];
     if(transfer_array_x[n_transfer-2]<0.)
       transfer_array_x[n_transfer-1]= 0.1*transfer_array_x[n_transfer-2];
     else
       transfer_array_x[n_transfer-1]=10.0*transfer_array_x[n_transfer-2];
     // Create interpolation array
     init_interpolate(transfer_array_x,transfer_array_y,n_transfer,gsl_interp_cspline,return_interp);
     // Free temporary arrays
     SID_free(SID_FARG transfer_array_x);
     SID_free(SID_FARG transfer_array_y);
   }
   else
     SID_trap_error("Transfer arrays must be >2 elements long.",ERROR_LOGIC);
   return(i_word);
}

