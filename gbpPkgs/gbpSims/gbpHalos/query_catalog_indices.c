#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>

int main(int argc, char *argv[]){
  char    filename_root[256];
  char    filename_properties[256];
  char    filename_indices[256];
  char    filename_out[256];
  char    prefix_text[5];
  FILE   *fp_properties=NULL;
  FILE   *fp_profiles  =NULL;
  FILE   *fp_out       =NULL;
  int     i_file;
  int     n_files;
  int     n_groups_all;
  int     i_group;
  int     i_group_selected;
  int     i_profile;
  int     flag_process_group;
  int     snap_number;
  int     n_groups_properties;
  int     n_groups_profiles;  
  halo_properties_info properties;
  halo_profile_info    profile;
  float  lambda,v_c;
  float  offset_COM;
  float  r_min,r_max;

  SID_init(&argc,&argv,NULL,NULL);

  strcpy(filename_root, argv[1]);
  snap_number     =atoi(argv[2]);
  i_group_selected=atoi(argv[3]);

  if(i_group_selected<0){
    flag_process_group=TRUE;
    i_group_selected*=-1;
    sprintf(prefix_text,"");
  }
  else{
    flag_process_group=FALSE;
    sprintf(prefix_text,"sub");
  }

  if(SID.I_am_Master){
    sprintf(filename_indices,"%s_%03d.catalog_%sgroups_indices",filename_root,snap_number,prefix_text);

    FILE *fp_in;
    if((fp_in=fopen(filename_indices,"r"))==NULL)
       SID_trap_error("Could not open file {%s}.",ERROR_IO_OPEN,filename_indices);
    int i_file,n_files,n_groups,n_groups_all;
    fread(&i_file,      sizeof(int),1,fp_in);
    fread(&n_files,     sizeof(int),1,fp_in);
    fread(&n_groups,    sizeof(int),1,fp_in);
    fread(&n_groups_all,sizeof(int),1,fp_in);
    int n_particles_i;
    for(i_group=0;i_group<i_group_selected;i_group++){
       fread(&n_particles_i,sizeof(int),1,fp_in);
       fseeko(fp_in,(off_t)(sizeof(size_t)*n_particles_i),SEEK_CUR);
    }
    fread(&n_particles_i,sizeof(int),1,fp_in);
    fprintf(stderr,"n_particles=%d\n",n_particles_i);
    int i_particle;
    for(i_particle=0;i_particle<n_particles_i;i_particle++){
       size_t index_i;
       fread(&index_i,sizeof(size_t),1,fp_in);
       fprintf(stderr,"%4d %lld\n",i_particle,index_i);
    }
   
    fclose(fp_in);

  }  

  SID_exit(ERROR_NONE);
}
