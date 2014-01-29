#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL);

  SID_log("Constructing match catalog...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Parse arguments
  char filename_catalog_root[MAX_FILENAME_LENGTH];
  char filename_matches_root[MAX_FILENAME_LENGTH];
  char filename_out[MAX_FILENAME_LENGTH];
  char prefix_text[32];
  int  flag_matches_type;
  int  i_read;
  int  j_read;
  int  catalog_read_mode;
  int  matches_read_mode;
  strcpy(filename_catalog_root,argv[1]);
  strcpy(filename_matches_root,argv[2]);
  flag_matches_type      =atoi(argv[3]);
  strcpy(prefix_text,          argv[4]);
  i_read                 =atoi(argv[5]);
  j_read                 =atoi(argv[6]);
  strcpy(filename_out,         argv[7]);
  if(!strcpy(prefix_text,"subgroup") ||
     !strcpy(prefix_text,"subgroups") ||
     !strcpy(prefix_text,"sub")){
     sprintf(prefix_text,"sub");
     catalog_read_mode=READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES;
     matches_read_mode=MATCH_SUBGROUPS;
  }
  else if(!strcpy(prefix_text,"group") ||
          !strcpy(prefix_text,"groups")){
     sprintf(prefix_text,"");
     catalog_read_mode=READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES;
     matches_read_mode=MATCH_GROUPS;
  }
  else
     SID_trap_error("Invalid catalog type (%s).",ERROR_SYNTAX,prefix_text);

  // Set filenames 
  char filename_cat1[MAX_FILENAME_LENGTH];
  char filename_cat2[MAX_FILENAME_LENGTH];
  sprintf(filename_cat1,"%s_%03d.catalog_%sgroups_properties",filename_catalog_root,prefix_text,i_read);
  sprintf(filename_cat2,"%s_%03d.catalog_%sgroups_properties",filename_catalog_root,prefix_text,j_read);

  // Contents of the halo properties structure
  //struct halo_properties_info{
  //  long long id_MBP;                    // ID of most bound particle in structure
  //  int       n_particles;               // Number of particles in the structure
  //  float     position_COM[3];           // Centre-of-mass position      [Mpc/h]
  //  float     position_MBP[3];           // Most bound particle position [Mpc/h]
  //  float     velocity_COM[3];           // Centre-of-mass velocity      [km/s]
  //  float     velocity_MBP[3];           // Most bound particle velocity [km/s]
  //  double    M_vir;                     // Bryan & Norman (ApJ 495, 80, 1998) virial mass [M_sol/h]
  //  float     R_vir;                     // Virial radius [Mpc/h]
  //  float     R_halo;                    // Distance of last halo particle from MBP [Mpc/h]
  //  float     R_max;                     // Radius of maximum circular velocity     [Mpc/h]
  //  float     V_max;                     // Maximum circular velocity               [km/s]
  //  float     sigma_v;                   // Total 3D velocity dispersion            [km/s]
  //  float     spin[3];                   // Specific angular momentum vector        [Mpc/h*km/s]
  //  float     q_triaxial;                // Triaxial shape parameter q=b/a
  //  float     s_triaxial;                // Triaxial shape parameter s=c/a
  //  float     shape_eigen_vectors[3][3]; // Normalized triaxial shape eigenvectors
  //};

  // Read matches
  int    *n_subgroups  =NULL;
  int    *n_groups     =NULL;
  int    *n_particles_i=NULL;
  int    *n_particles_j=NULL;
  int    *match_ids    =NULL;
  float  *match_score  =NULL;
  size_t *match_index  =NULL;
  int     n_halos_i;
  int     n_halos_j;
  int     n_files;
  int     n_subgroups_max;
  int     n_groups_max;
  int     n_halos_max;
  read_matches_header(filename_matches_root,
                      0,
                      MAX(i_read,j_read),
                      1,
                      &n_files,
                      &n_subgroups,
                      &n_groups,
                      &n_subgroups_max,
                      &n_groups_max,
                      &n_halos_max);
  read_matches(filename_matches_root,
               i_read,
               j_read,
               n_halos_max,
               matches_read_mode,
               &n_halos_i,
               &n_halos_j,
               n_particles_i,
               n_particles_j,
               NULL,
               NULL,
               match_ids,
               match_score,
               match_index);

  // Create a storage array mapping the indices of the second catalog
  //    to those of the halos they are matched to in the first catalog
  int  i_halo;
  int  j_halo;
  int *storage_index;
  storage_index=(int *)SID_malloc(sizeof(int)*n_halos_j);
  for(j_halo=0;j_halo<n_halos_j;j_halo++)
     storage_index[j_halo]=-1;
  for(i_halo=0,j_halo=0;j_halo<n_halos_j && i_halo<n_halos_i;j_halo++){
     while(match_ids[match_index[i_halo]]<j_halo && i_halo<(n_halos_i-1)) i_halo++;
     if(match_ids[match_index[i_halo]]<j_halo)                            i_halo++;
     if(match_ids[match_index[i_halo]]==j_halo)
        storage_index[j_halo]=match_index[i_halo];
  }

  // Open catalog files
  fp_catalog_info fp_properties_i;
  fp_catalog_info fp_properties_j;
  fopen_catalog(filename_catalog_root,
                i_read,
                catalog_read_mode,
                &fp_properties_i);
  fopen_catalog(filename_catalog_root,
                j_read,
                catalog_read_mode,
                &fp_properties_j);

  // Read catalogs
  halo_properties_info *properties_i;
  halo_properties_info *properties_j;
  properties_i    =(halo_properties_info *)SID_malloc(sizeof(halo_properties_info)*n_halos_i);
  properties_j    =(halo_properties_info *)SID_malloc(sizeof(halo_properties_info)*n_halos_i);
  for(i_halo=0;i_halo<n_halos_i;i_halo++)
     fread_catalog_file(&fp_properties_i,NULL,&(properties_i[i_halo]),NULL,i_halo);
  for(j_halo=0;j_halo<n_halos_j;j_halo++){
     if(storage_index[j_halo]>=0)
        fread_catalog_file(&fp_properties_j,NULL,&(properties_j[storage_index[j_halo]]),NULL,j_halo);
  }
  fclose_catalog(&fp_properties_i);
  fclose_catalog(&fp_properties_j);

  // Write results
  int   i_column=0;
  FILE *fp_out;
  fp_out=fopen(filename_out,"w");
  fprintf(fp_out,"# Catalog for matches {root %s} and catalog {root %s}; snap No. %d to %d.\n",
                 filename_matches_root,
                 filename_catalog_root,
                 i_read,j_read);
  fprintf(fp_out,"# Columns:(%02d) id (catalog No. 1)\n",          i_column++);
  fprintf(fp_out,"#         (%02d) id (catalog No. 2)\n",          i_column++);
  fprintf(fp_out,"#         (%02d)  M (catalog No. 1) [M_sol/h]\n",i_column++);
  fprintf(fp_out,"#         (%02d)  M (catalog No. 2) [M_sol/h]\n",i_column++);
  fprintf(fp_out,"#         (%02d)  x (catalog No. 1) [Mpc/h]\n",  i_column++);
  fprintf(fp_out,"#         (%02d)  y (catalog No. 1) [Mpc/h]\n",  i_column++);
  fprintf(fp_out,"#         (%02d)  z (catalog No. 1) [Mpc/h]\n",  i_column++);
  fprintf(fp_out,"#         (%02d)  x (catalog No. 2) [Mpc/h]\n",  i_column++);
  fprintf(fp_out,"#         (%02d)  y (catalog No. 2) [Mpc/h]\n",  i_column++);
  fprintf(fp_out,"#         (%02d)  z (catalog No. 2) [Mpc/h]\n",  i_column++);
  for(i_halo=0;i_halo<n_halos_i;i_halo++){
     fprintf(fp_out,"%10d %10d %10.4le %10.4le %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                    i_halo,
                    match_ids[i_halo],
                    properties_i[i_halo].M_vir,
                    properties_j[i_halo].M_vir,
                    properties_i[i_halo].position_COM[0],
                    properties_i[i_halo].position_COM[1],
                    properties_i[i_halo].position_COM[2],
                    properties_j[i_halo].position_COM[0],
                    properties_j[i_halo].position_COM[1],
                    properties_j[i_halo].position_COM[2]);
  }
  fclose(fp_out);

  // Clean-up
  SID_free(SID_FARG n_subgroups);
  SID_free(SID_FARG n_groups);
  SID_free(SID_FARG n_particles_i);
  SID_free(SID_FARG n_particles_j);
  SID_free(SID_FARG properties_i);
  SID_free(SID_FARG properties_j);
  SID_free(SID_FARG match_ids);
  SID_free(SID_FARG match_score);
  SID_free(SID_FARG match_index);
  SID_free(SID_FARG storage_index);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

