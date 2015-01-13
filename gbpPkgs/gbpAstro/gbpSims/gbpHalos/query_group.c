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
  char    prefix_text[5];
  int     i_file;
  int     n_files;
  int     n_groups_all;
  int     i_group;
  int     i_group_selected;
  int     i_profile;
  int     flag_process_group;
  int     snap_number;
  char   *group_type_in;

  SID_init(&argc,&argv,NULL,NULL);

  strcpy(filename_root, argv[1]);
  snap_number     =atoi(argv[2]);
  group_type_in   =     argv[3];
  i_group_selected=atoi(argv[4]);

  if(!strcmp(group_type_in,"sub")||
     !strcmp(group_type_in,"subgroup")){
     flag_process_group=FALSE;
     sprintf(prefix_text,"sub");
  }
  else if(!strcmp(group_type_in,"group")){
     flag_process_group=TRUE;
    sprintf(prefix_text,"");
  }
  else
     SID_trap_error("Incorrect syntax.",ERROR_SYNTAX);

  char   filename_subgroups[MAX_FILENAME_LENGTH];
  char   filename_groups[MAX_FILENAME_LENGTH];
  char   filename_particles[MAX_FILENAME_LENGTH];
  FILE  *fp_groups;
  FILE  *fp_particles;
  int    n_particles_i;
  size_t particle_offset_i;
  int    id_byte_size;
  int    n_particles_int;
  size_t n_particles; 
  int    particle_offset;
  int    MBP_ID_int;
  size_t MBP_ID;
  int    n_sub_i;

  if(SID.I_am_Master){
    sprintf(filename_subgroups,"%s_%03d.catalog_subgroups",filename_root,snap_number,prefix_text);
    sprintf(filename_groups,   "%s_%03d.catalog_groups",   filename_root,snap_number,prefix_text);
    sprintf(filename_particles,"%s_%03d.catalog_particles",filename_root,snap_number,prefix_text);
    SID_log("Selected subgroup  file:{%s}",SID_LOG_COMMENT,filename_subgroups);
    SID_log("Selected group     file:{%s}",SID_LOG_COMMENT,filename_groups);
    SID_log("Selected particles file:{%s}",SID_LOG_COMMENT,filename_particles);
    SID_log("", SID_LOG_COMMENT);

    // Open file and read header
    int n_bytes_groups;
    if(flag_process_group)
       fp_groups=fopen(filename_groups,"r");
    else
       fp_groups=fopen(filename_subgroups,"r");
    fread(&n_groups_all,  sizeof(int),1,fp_groups);
    fread(&n_bytes_groups,sizeof(int),1,fp_groups);
    if(n_bytes_groups!=4 && n_bytes_groups!=8)
       SID_trap_error("Invalid group offset byte length (%d).",ERROR_LOGIC,n_bytes_groups);

    // Sanity check
    if(i_group_selected<0 || i_group_selected>=n_groups_all)
       SID_trap_error("Invalid group selection {%d;n_groups=%d}.",ERROR_LOGIC,i_group_selected,n_groups_all);

    // Find the group we want and get the needed info
    int particle_offset_i_int;
    fseeko(fp_groups,sizeof(int)*i_group_selected,SEEK_CUR);
    fread(&n_particles_i,sizeof(int),1,fp_groups);
    fseeko(fp_groups,sizeof(int)*(n_groups_all-i_group_selected-1),SEEK_CUR);
    fseeko(fp_groups,n_bytes_groups*i_group_selected,SEEK_CUR);
    if(n_bytes_groups==sizeof(int)){
       fread(&particle_offset_i_int,sizeof(int),1,fp_groups);
       particle_offset_i=(size_t)particle_offset_i_int;
    }
    else
       fread(&particle_offset_i,sizeof(size_t),1,fp_groups);
    fseeko(fp_groups,n_bytes_groups*(n_groups_all-i_group_selected-1),SEEK_CUR);
    if(flag_process_group){
       fseeko(fp_groups,sizeof(int)*i_group_selected,SEEK_CUR);
       fread(&n_sub_i,sizeof(int),1,fp_groups);
       fseeko(fp_groups,sizeof(int)*(n_groups_all-i_group_selected-1),SEEK_CUR);
    }
    else
       n_sub_i=0;
    fclose(fp_groups);

    // Fetch the MBP ID
    fp_particles=fopen(filename_particles,"r");
    fread(&id_byte_size,sizeof(int),1,fp_particles);
    if(id_byte_size==sizeof(int)){
       fread(&n_particles_int,sizeof(int),1,fp_particles);
       n_particles=(size_t)n_particles_int;
    }
    else
       fread(&n_particles,sizeof(size_t),1,fp_particles);
    if(((size_t)particle_offset_i)>=n_particles)
       SID_trap_error("Invalid particle offset {%lld;n_particles=%lld}.",ERROR_LOGIC,particle_offset_i,n_particles);
    fseeko(fp_particles,id_byte_size*particle_offset_i,SEEK_CUR);
    if(id_byte_size==sizeof(int)){
       fread(&MBP_ID_int,sizeof(int),1,fp_particles);
       MBP_ID=(size_t)MBP_ID_int;
    }
    else
       fread(&MBP_ID,sizeof(size_t),1,fp_particles);
    fclose(fp_particles);

    SID_log("snapshot no.    = %d",SID_LOG_COMMENT,snap_number);
    SID_log("halo            = %d of %d",    SID_LOG_COMMENT,i_group_selected,n_groups_all);
    SID_log("id byte size    = %d",SID_LOG_COMMENT,id_byte_size);
    SID_log("particle offset = %lld",SID_LOG_COMMENT,particle_offset_i);
    SID_log("no. of particles= %lld of %lld",SID_LOG_COMMENT,n_particles_i,n_particles);
    SID_log("MBP ID          = %d",SID_LOG_COMMENT,MBP_ID);
    if(flag_process_group)
       SID_log("no. of subhalos = %d",SID_LOG_COMMENT,n_sub_i);

  }  

  SID_exit(ERROR_NONE);
}
