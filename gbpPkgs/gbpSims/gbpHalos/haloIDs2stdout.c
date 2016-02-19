#define _MAIN
#include <gbpLib.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc,char *argv[]){

   char filename_root[MAX_FILENAME_LENGTH];
   char filename_groups[MAX_FILENAME_LENGTH];
   char filename_subgroups[MAX_FILENAME_LENGTH];
   char filename_ids[MAX_FILENAME_LENGTH];
   sprintf(filename_root,"%s",argv[1]);
   int i_snapshot=atoi(argv[2]);
   sprintf(filename_groups,   "%s_%03d.catalog_groups",   filename_root,i_snapshot);
   sprintf(filename_subgroups,"%s_%03d.catalog_subgroups",filename_root,i_snapshot);
   sprintf(filename_ids,      "%s_%03d.catalog_particles",filename_root,i_snapshot);
   
   FILE *fp_groups   =fopen(filename_groups,   "r");
   FILE *fp_subgroups=fopen(filename_subgroups,"r");
   FILE *fp_ids      =fopen(filename_ids,      "r");

   int    n_groups,group_offset_byte_size,n_subgroups,subgroup_offset_byte_size,n_ids_i,id_byte_size;
   size_t n_ids;
   fread_verify(&n_groups,                 sizeof(int),1,fp_groups);
   fread_verify(&group_offset_byte_size,   sizeof(int),1,fp_groups);
   fread_verify(&n_subgroups,              sizeof(int),1,fp_subgroups);
   fread_verify(&subgroup_offset_byte_size,sizeof(int),1,fp_subgroups);
   fread_verify(&id_byte_size,sizeof(int),1,fp_ids);
   if(id_byte_size==sizeof(int)){
      fread_verify(&n_ids_i,sizeof(int),1,fp_ids);
      n_ids=(size_t)n_ids_i;
   }
   else
      fread_verify(&n_ids,sizeof(size_t),1,fp_ids);

   size_t i_group   =0;
   size_t i_subgroup=0;
   for(size_t i_particle=0;i_particle<n_ids;i_particle++){
      long int id_i;
      if(id_byte_size==sizeof(int)){
         int id_i_int;
         fread_verify(&id_i_int,sizeof(int),1,fp_ids);
         id_i=(long int)id_i_int;
      }
      else
         fread_verify(&id_i,sizeof(long int),1,fp_ids);
      fprintf(stdout,"%zd %zd %zd %lld\n",i_particle,i_group,i_subgroup,id_i);
   }
   fclose(fp_groups);
   fclose(fp_subgroups);
   fclose(fp_ids);
}

