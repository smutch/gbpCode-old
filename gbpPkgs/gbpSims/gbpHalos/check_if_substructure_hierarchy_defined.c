#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHalos.h>

int check_if_substructure_hierarchy_defined(const char *filename_SSimPL_root,const char *filename_halos_version,int i_snap){
  // Perform check
  int flag_read_sub_pointers=FALSE;
  if(SID.I_am_Master){
     char  filename_in[MAX_FILENAME_LENGTH];
     sprintf(filename_in,"%s/halos/%s_%03d.catalog_subgroups",filename_SSimPL_root,filename_halos_version,i_snap);
     FILE *fp_test=fopen(filename_in,"r");
     int n_subgroups_in;fread_verify(&n_subgroups_in,sizeof(int),1,fp_test);
     int offset_size;   fread_verify(&offset_size,   sizeof(int),1,fp_test);
     fseeko(fp_test,n_subgroups_in*(sizeof(int)+offset_size),SEEK_CUR);
     int test;          fread(&offset_size,sizeof(int),1,fp_test); // make sure not to verify here!
     flag_read_sub_pointers=!feof(fp_test);
     fclose(fp_test);
  }
  SID_Bcast(&flag_read_sub_pointers,sizeof(int),MASTER_RANK,SID.COMM_WORLD);
  if(flag_read_sub_pointers)
     SID_log("Substructure hierarchy pointers present and will be used.",SID_LOG_COMMENT);
  else
     SID_log("Substructure hierarchy pointers not present and will not be used.",SID_LOG_COMMENT);
  return(flag_read_sub_pointers);
}
