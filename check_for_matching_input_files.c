#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int check_for_matching_input_files(const char *filename_root_in,int i_read);
int check_for_matching_input_files(const char *filename_root_in,int i_read){
   int   flag_all_inputs_present=TRUE;
   char  filename_test[MAX_FILENAME_LENGTH];
   FILE *fp_test;
   // Test from groups file
   sprintf(filename_test,"%s_%03d.catalog_groups",   filename_root_in,i_read);
   if((fp_test=fopen(filename_test,"r"))==NULL)
      flag_all_inputs_present=FALSE;
   else
      fclose(fp_test);
   // Test from subgroups file
   sprintf(filename_test,"%s_%03d.catalog_subgroups",filename_root_in,i_read);
   if((fp_test=fopen(filename_test,"r"))==NULL)
      flag_all_inputs_present=FALSE;
   else
      fclose(fp_test);
   // Test from particles file
   sprintf(filename_test,"%s_%03d.catalog_particles",filename_root_in,i_read);
   if((fp_test=fopen(filename_test,"r"))==NULL)
      flag_all_inputs_present=FALSE;
   else
      fclose(fp_test);
   // Test from particles file
   sprintf(filename_test,"%s_%03d.catalog_PHKs",filename_root_in,i_read);
   if((fp_test=fopen(filename_test,"r"))==NULL)
      flag_all_inputs_present=FALSE;
   else
      fclose(fp_test);
   return(flag_all_inputs_present);
}

