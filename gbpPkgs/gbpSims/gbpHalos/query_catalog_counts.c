#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>

int main(int argc, char *argv[]){
  SID_init(&argc,&argv,NULL,NULL);

  char filename_SSimPL[MAX_FILENAME_LENGTH];
  char filename_halo_type[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL,   argv[1]);
  strcpy(filename_halo_type,argv[2]);

  if(SID.I_am_Master){
    FILE *fp_out=stdout;

    // Count the number of snapshots and read their expansion factors
    char filename_a_list[MAX_FILENAME_LENGTH];
    sprintf(filename_a_list,"%s/run/a_list.txt",filename_SSimPL);
    FILE  *fp_a_list  =fopen(filename_a_list,"r");
    int    n_snaps    =count_lines_data(fp_a_list);
    float *a_list     =(float *)SID_malloc(sizeof(float *)*n_snaps);
    char  *line       =NULL;
    size_t line_length=0;
    for(int i_snap=0;i_snap<n_snaps;i_snap++){
       grab_next_line_data(fp_a_list,&line,&line_length);
       grab_float(line,1,&(a_list[i_snap]));
    }
    SID_free(SID_FARG line);
    fclose(fp_a_list);

    for(int i_snap=0;i_snap<n_snaps;i_snap++){
       char filename_group_properties[MAX_FILENAME_LENGTH];
       char filename_subgroup_properties[MAX_FILENAME_LENGTH];
       sprintf(filename_group_properties,   "%s/catalogs/%s",filename_SSimPL,filename_halo_type,i_snap);
       sprintf(filename_subgroup_properties,"%s/catalogs/%s",filename_SSimPL,filename_halo_type,i_snap);
       fp_catalog_info fp_group_properties;
       fp_catalog_info fp_subgroup_properties;
       fopen_catalog(filename_group_properties,
                     i_snap,
                     READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES,
                     &fp_group_properties);
       fopen_catalog(filename_subgroup_properties,
                     i_snap,
                     READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES,
                     &fp_subgroup_properties);

       // Write header
       if(i_snap==0){
          int i_column=1;
          fprintf(fp_out,"# Column (%02d): Snapshot\n",        i_column++);
          fprintf(fp_out,"#        (%02d): Expansion factor\n",i_column++);
          fprintf(fp_out,"#        (%02d): Redshift\n",        i_column++);
          fprintf(fp_out,"#        (%02d): No. of groups\n",   i_column++);
          fprintf(fp_out,"#        (%02d): No. of subgroups\n",i_column++);
       }
       fprintf(fp_out,"%03d %le %7.4f %d %d\n",i_snap,a_list[i_snap],z_of_a(a_list[i_snap]),fp_group_properties.n_halos_total,fp_subgroup_properties.n_halos_total);

       // Close file
       fclose_catalog(&fp_group_properties);
       fclose_catalog(&fp_subgroup_properties);
    }
  }  

  SID_exit(ERROR_NONE);
}
