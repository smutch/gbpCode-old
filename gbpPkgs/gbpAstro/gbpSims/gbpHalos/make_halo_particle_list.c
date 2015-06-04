#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpSPH.h>

int main(int argc, char *argv[]){
  SID_init(&argc,&argv,NULL,NULL);

  // Parse command line
  select_gadget_ids_params_info params;
  int  snapshot;
  int  halo_index;
  int  halo_type;
  int  n_files_out;
  int  select_mode;
  char filename_SSimPL_root[MAX_FILENAME_LENGTH];
  char filename_halo_version[MAX_FILENAME_LENGTH];
  char filename_in_root[MAX_FILENAME_LENGTH];
  char filename_out_root[MAX_FILENAME_LENGTH];
  GBPREAL cen_select[3];
  GBPREAL select_size;
  strcpy(filename_SSimPL_root,      argv[1]);
  strcpy(filename_halo_version,     argv[2]);
  snapshot                  =  atoi(argv[3]);
  halo_index                =  atoi(argv[4]);
  halo_type                 =  atoi(argv[5]);
  strcpy(filename_out_root,         argv[6]);
  n_files_out               =  atoi(argv[7]);

  SID_log("Writing halo particles to ascii file {%s;snapshot=%d;halo_index=%d}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_SSimPL_root,snapshot,halo_index);

  // Initialize the plist data structure 
  plist_info plist;
  params.plist=&plist;
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);

  // Read the halo ID list.  Generate sort indicies and copy the 
  //    list to a duplicate array.  Make sure this is the one added
  //    to params.  Pass this to write_gadget_ascii below, so that
  //    the particles get written in the same order that they are
  //    in the halo catalog.

  char   filename_halos[MAX_FILENAME_LENGTH];
  if(halo_type==0)
     sprintf(filename_halos,"%s/halos/%s_%03d.catalog_groups",   filename_SSimPL_root,filename_halo_version,snapshot,filename_halo_version);
  else
     sprintf(filename_halos,"%s/halos/%s_%03d.catalog_subgroups",filename_SSimPL_root,filename_halo_version,snapshot,filename_halo_version);
  FILE   *fp_groups=fopen(filename_halos,"r");
  int     n_groups;
  int     offset_size;
  int     halo_length;
  size_t  halo_offset;
  fread(&n_groups,   sizeof(int),1,fp_groups);
  fread(&offset_size,sizeof(int),1,fp_groups);
  fseeko(fp_groups,(off_t)(2*sizeof(int)+halo_index*sizeof(int)),SEEK_SET);
  fread(&halo_length,sizeof(int),1,fp_groups);
  fseeko(fp_groups,(off_t)(2*sizeof(int)+n_groups*sizeof(int)+halo_index*offset_size),SEEK_SET);
  if(offset_size==sizeof(int)){
     int halo_offset_i;
     fread(&halo_offset_i,offset_size,1,fp_groups);
     halo_offset=(size_t)halo_offset_i;
  }
  else
     fread(&halo_offset,offset_size,1,fp_groups);
  fclose(fp_groups);

  char filename_ids[MAX_FILENAME_LENGTH];
  sprintf(filename_ids,"%s/halos/%s_%03d.catalog_particles",filename_SSimPL_root,filename_halo_version,snapshot,filename_halo_version);
  FILE  *fp_ids=fopen(filename_ids,"r");
  int    id_byte_size;
  size_t n_ids;
  fread(&id_byte_size,sizeof(int),1,fp_ids);
  SID_log("%d %d-byte IDs to be read (offset=%d)",SID_LOG_COMMENT,halo_length,id_byte_size,halo_offset);
  if(id_byte_size==sizeof(int)){
     int n_ids_i;
     fread(&n_ids_i,sizeof(int),1,fp_ids);
     n_ids=(size_t)n_ids_i;
  }
  else
     fread(&n_ids,sizeof(size_t),1,fp_ids);
  fseeko(fp_ids,(off_t)(sizeof(int)+id_byte_size+halo_offset*id_byte_size),SEEK_SET);
  params.n_ids  =halo_length;
  params.id_list          =(size_t *)SID_malloc(sizeof(size_t)*halo_length);
  size_t *id_list_unsorted=(size_t *)SID_malloc(sizeof(size_t)*halo_length);
  int flag_long_ids=TRUE;
  if(id_byte_size==sizeof(int)){
     flag_long_ids=FALSE;
     int *id_list_i=(int *)SID_malloc(sizeof(int)*halo_length);
     fread(id_list_i,id_byte_size,halo_length,fp_ids);
     for(int i_p=0;i_p<halo_length;i_p++)
        id_list_unsorted[i_p]=(size_t)id_list_i[i_p];
     SID_free(SID_FARG id_list_i);
  }
  else
     fread(id_list_unsorted,id_byte_size,halo_length,fp_ids);
  fclose(fp_ids);
  memcpy(params.id_list,id_list_unsorted,sizeof(size_t)*halo_length);
  merge_sort(params.id_list,halo_length,NULL,SID_SIZE_T,SORT_INPLACE_ONLY,SORT_COMPUTE_INPLACE);

  // Count the particles
  size_t n_particles_type_local[N_GADGET_TYPE];
  size_t n_particles_type[N_GADGET_TYPE];
  int    flag_long_IDs;
  sprintf(filename_in_root,"%s/snapshots/snapshot",filename_SSimPL_root);
  process_gadget_file("Counting particles in selection...",
                      filename_in_root,
                      snapshot,
                      select_gadget_all,
                      process_gadget_file_fctn_null,
                      &params,
                      n_particles_type_local,
                      n_particles_type,
                      &flag_long_IDs,
                      PROCESS_GADGET_BINARY_DEFAULT);

  // Allocate RAM for the particles
  allocate_gadget_particles(&plist,n_particles_type_local,n_particles_type,flag_long_IDs);

  // Read the particles
  process_gadget_file("Performing read...",
                      filename_in_root,
                      snapshot,
                      select_gadget_all,
                      store_gadget_particles,
                      &params,
                      NULL,
                      NULL,
                      &flag_long_IDs,
                      PROCESS_GADGET_BINARY_DEFAULT);

  // Write the snapshot
  if(SID.I_am_Master){
     char filename_out[MAX_FILENAME_LENGTH];
     sprintf(filename_out,"%s_%03d_%08d.ascii",filename_out_root,snapshot,halo_index);
     FILE *fp=fopen(filename_out,"w");
     fprintf(fp,"#Columns:\n");
     fprintf(fp,"#  1) Gadget particle type\n");
     fprintf(fp,"#  2) x   [Mpc/h]\n");
     fprintf(fp,"#  3) y   [Mpc/h]\n");
     fprintf(fp,"#  4) z   [Mpc/h]\n");
     fprintf(fp,"#  5) v_x [km/s]\n");
     fprintf(fp,"#  6) v_y [km/s]\n");
     fprintf(fp,"#  7) v_z [km/s]\n");
     fprintf(fp,"#  8) id\n");
     int i_species=GADGET_TYPE_DARK;
     size_t  n_p=((size_t *)ADaPS_fetch( plist.data,"n_%s", plist.species[i_species]))[0];
     GBPREAL *x =(GBPREAL  *)ADaPS_fetch(plist.data,"x_%s", plist.species[i_species]);
     GBPREAL *y =(GBPREAL  *)ADaPS_fetch(plist.data,"y_%s", plist.species[i_species]);
     GBPREAL *z =(GBPREAL  *)ADaPS_fetch(plist.data,"z_%s", plist.species[i_species]);
     GBPREAL *vx=(GBPREAL  *)ADaPS_fetch(plist.data,"vx_%s",plist.species[i_species]);
     GBPREAL *vy=(GBPREAL  *)ADaPS_fetch(plist.data,"vy_%s",plist.species[i_species]);
     GBPREAL *vz=(GBPREAL  *)ADaPS_fetch(plist.data,"vz_%s",plist.species[i_species]);
     size_t  *id=(size_t   *)ADaPS_fetch(plist.data,"id_%s",plist.species[i_species]);
     size_t  *id_indices=NULL;
     SID_log("Sorting IDs...",SID_LOG_OPEN);
     merge_sort(id,n_p,&id_indices,SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
     SID_log("Done.",SID_LOG_CLOSE);
     SID_log("Writing particles...",SID_LOG_OPEN);
     pcounter_info pcounter;
     SID_init_pcounter(&pcounter,halo_length,10);
     int n_unfound=0;
     for(int i_p=0;i_p<halo_length;i_p++){
       size_t k_p=id_indices[find_index(id,id_list_unsorted[i_p],n_p,id_indices)];
       if(id[k_p]!=id_list_unsorted[i_p])
          n_unfound++;
       else
          fprintf(fp,"%1d %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %7zd\n",
                  i_species,
                  (double)(x[k_p]),
                  (double)(y[k_p]),
                  (double)(z[k_p]),
                  (double)(vx[k_p]),
                  (double)(vy[k_p]),
                  (double)(vz[k_p]),
                  id[k_p]);
       SID_check_pcounter(&pcounter,i_p);
     }
     fclose(fp);
     if(n_unfound!=0)
        SID_log("(%d unfound)...",SID_LOG_CONTINUE,n_unfound);
     SID_log("Done.",SID_LOG_CLOSE);
     SID_free(SID_FARG id_indices);
  }

  // Clean-up 
  SID_free(SID_FARG id_list_unsorted);
  SID_free(SID_FARG params.id_list);
  free_plist(&plist);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

