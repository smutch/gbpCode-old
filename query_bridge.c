#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){
  int     i_halo;
  char    group_text_prefix[4];
  int     n_files;
  int     k_read;
  int     max_n_groups;
  int     l_read;
  int     n_groups;
  int    *n_particles_i;
  int    *n_particles_j;
  int     mode;
  int     n_groups_i;
  int     n_groups_j;
  int    *match_ids;
  float  *match_score;
  size_t *match_index;
  int     j_halo;
  int     match;
  int     i_read;
  int     i_read_in;
  int     i_read_start;
  int     i_read_stop;
  SID_fp  fp_in;

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  char    filename_SSimPL_root[MAX_FILENAME_LENGTH];
  char    filename_halos_root[MAX_FILENAME_LENGTH];
  char    filename_trees_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_root,argv[1]);
  strcpy(filename_halos_root, argv[2]);
  strcpy(filename_trees_root, argv[3]);
  if(!strcmp(argv[4],"groups") || !strcmp(argv[4],"group"))
     mode=MATCH_GROUPS;
  else if(!strcmp(argv[4],"subgroups") || !strcmp(argv[4],"subgroup"))
     mode=MATCH_SUBGROUPS;
  else{
     SID_log("Invalid mode selection {%s}.  Should be 'group' or 'subgroup'.",SID_LOG_COMMENT,argv[4]);
     SID_exit(ERROR_SYNTAX);
  }
  i_read_in=atoi(argv[5]);
  i_halo   =atoi(argv[6]);
  SID_log("Querying bridge matches for halo #%d in file #%d from {%s}...",SID_LOG_OPEN,i_halo,i_read_in,filename_trees_root);

  // Convert filename_root to filename
  switch(mode){
     case MATCH_SUBGROUPS:
     sprintf(group_text_prefix,"sub");
     break;
     case MATCH_GROUPS:
     sprintf(group_text_prefix,"");
     break;
  }

  // Read tree header information
  tree_info *trees;
  char       filename_file_root[MAX_FILENAME_LENGTH];
  sprintf(filename_file_root,"%s/trees/%s",filename_SSimPL_root,filename_trees_root);
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
  init_trees_read(filename_file_root,TREE_READ_HEADER_ONLY,&trees);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

  // Check that the given snapshot is valid
  int i_file=0;
  int j_read=trees->snap_list[i_file];
  while(i_file<(trees->n_snaps-1) && j_read<i_read_in){
    i_file++;
    j_read=trees->snap_list[i_file];
  }
  if(i_read_in!=trees->snap_list[i_file])
     SID_trap_error("Invalid snapshot specified {%d}.",ERROR_LOGIC,i_read_in);

  // Open output file
  char filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_out,"%s_%d_%d_bridges.txt",filename_trees_root,i_read_in,i_halo);
  FILE *fp_out=fopen(filename_out,"w");

  // Write header
  int i_column=1;
  fprintf(fp_out,"# Column (%02d): Halo snapshot\n",              i_column++);
  fprintf(fp_out,"#        (%02d): Halo index\n",                 i_column++);
  fprintf(fp_out,"#        (%02d): Halo ID\n",                    i_column++);
  fprintf(fp_out,"#        (%02d): Halo log10(M_vir [M_sol/h])\n",i_column++);
  fprintf(fp_out,"#        (%02d): Halo tree ID\n",               i_column++);
  fprintf(fp_out,"#        (%02d): Descendant file offset\n",     i_column++);
  fprintf(fp_out,"#        (%02d): Descendant snapshot\n",        i_column++);
  fprintf(fp_out,"#        (%02d): Descendant index\n",           i_column++);
  fprintf(fp_out,"#        (%02d): Descendant ID\n",              i_column++);
  fprintf(fp_out,"#        (%02d): Bridge match snapshot\n",      i_column++);
  fprintf(fp_out,"#        (%02d): Bridge match index\n",         i_column++);
  fprintf(fp_out,"#        (%02d): Halo type\n",                  i_column++);
  fprintf(fp_out,"#        (%02d): Halo type string\n",           i_column++);

  int halo_file_offset =1;
  int flag_write_header=TRUE;
  for(i_read=MAX(0,i_read_in-trees->n_search*trees->i_read_step),i_file=MAX(0,i_file-trees->n_search);
      i_read<=MIN(trees->i_read_stop,i_read_in+trees->n_search*trees->i_read_step);
      i_read+=(trees->i_read_step),i_file++){
    //if(i_read!=i_read_in){
    if(TRUE){
       SID_log("Processing snapshot %03d (%03d of %03d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read,i_file+1,trees->n_snaps);
       int halo_id            =0;
       int halo_type          =0;
       int halo_descendant_id =0;
       int halo_tree_id       =0;
       int halo_index         =0;
       int bridge_match_id    =0;
       int bridge_match_file  =0;
       int bridge_match_index =0;

       // Open files
       fp_catalog_info      fp_properties;
       halo_properties_info properties;
       char                 filename_catalog_root[MAX_FILENAME_LENGTH];
       int                  read_props_mode;
       sprintf(filename_catalog_root,"%s/catalogs/%s",filename_SSimPL_root,filename_halos_root);
       if(mode==MATCH_GROUPS)
          read_props_mode=READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES;
       else
          read_props_mode=READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES;
       fopen_catalog(filename_catalog_root,
                     i_read,
                     read_props_mode,
                     &fp_properties);
       char filename_in[MAX_FILENAME_LENGTH];
       SID_fp fp_in_trees;
       SID_fp fp_in_bridge_match;
       sprintf(filename_in,"%s/trees/%s/horizontal/trees/horizontal_trees_%03d.dat",filename_SSimPL_root,filename_trees_root,i_read);
       SID_fopen(filename_in,"r",&fp_in_trees);
       if(i_read<i_read_in)
          sprintf(filename_in,"%s/trees/%s/horizontal/trees/horizontal_trees_bridge_forematch_pointers_%03d.dat",
                              filename_SSimPL_root,filename_trees_root,i_read);
       else
          sprintf(filename_in,"%s/trees/%s/horizontal/trees/horizontal_trees_bridge_backmatch_pointers_%03d.dat",
                              filename_SSimPL_root,filename_trees_root,i_read);
       SID_fopen(filename_in,"r",&fp_in_bridge_match);

       // Perform search
       int n_groups;
       int n_subgroups;
       int n_step_in;
       int n_search_in;
       int n_groups_max_in;
       int n_subgroups_max_in;
       int n_trees_group_in;
       int n_trees_subgroup_in;
       SID_fread_all(&n_step_in,          sizeof(int),1,&fp_in_trees);
       SID_fread_all(&n_search_in,        sizeof(int),1,&fp_in_trees);
       SID_fread_all(&n_groups,           sizeof(int),1,&fp_in_trees);
       SID_fread_all(&n_subgroups,        sizeof(int),1,&fp_in_trees);
       SID_fread_all(&n_groups_max_in,    sizeof(int),1,&fp_in_trees);
       SID_fread_all(&n_subgroups_max_in, sizeof(int),1,&fp_in_trees);
       SID_fread_all(&n_trees_subgroup_in,sizeof(int),1,&fp_in_trees);
       SID_fread_all(&n_trees_group_in,   sizeof(int),1,&fp_in_trees);
       if(i_read==i_read_in){
          if(mode==MATCH_GROUPS         && i_halo>=n_groups)    SID_trap_error("Invalid group halo index (ie. %d>=%d).",   ERROR_LOGIC,i_halo,n_groups);
          else if(mode==MATCH_SUBGROUPS && i_halo>=n_subgroups) SID_trap_error("Invalid subgroup halo index (ie. %d>=%d).",ERROR_LOGIC,i_halo,n_subgroups);
       }
       SID_fskip(sizeof(int),8,&fp_in_bridge_match);
       int i_group;
       int i_subgroup;
       for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
          // Read groups
          int n_subgroups_group;
          if(mode==MATCH_GROUPS){
             SID_fread_all(&halo_id,           sizeof(int),1,&fp_in_trees);
             SID_fread_all(&halo_type,         sizeof(int),1,&fp_in_trees);
             SID_fread_all(&halo_descendant_id,sizeof(int),1,&fp_in_trees);
             SID_fread_all(&halo_tree_id,      sizeof(int),1,&fp_in_trees);
             SID_fread_all(&halo_file_offset,  sizeof(int),1,&fp_in_trees);
             SID_fread_all(&halo_index,        sizeof(int),1,&fp_in_trees);
             SID_fread_all(&bridge_match_id,   sizeof(int),1,&fp_in_bridge_match);
             SID_fread_all(&bridge_match_file, sizeof(int),1,&fp_in_bridge_match);
             SID_fread_all(&bridge_match_index,sizeof(int),1,&fp_in_bridge_match);
             // Write match
             int flag_is_halo=(i_read==i_read_in && i_halo==i_group);
             if(trees->snap_list[bridge_match_file]==i_read_in && bridge_match_index==i_halo || flag_is_halo){
                char *halo_type_string=NULL;
                tree_case_flags_text(halo_type,"+",&halo_type_string);
                fread_catalog_file(&fp_properties,NULL,&properties,NULL,i_subgroup);
                int descendant_snap;
                int bridge_match_snap;
                if(halo_file_offset>0)
                   descendant_snap=trees->snap_list[i_file+halo_file_offset];
                else
                   descendant_snap=-1;
                if(bridge_match_file>=0)
                   bridge_match_snap=trees->snap_list[bridge_match_file];
                else
                   bridge_match_snap=-1;
                fprintf(fp_out,"%3d %7d %7d %5.2lf %7d %3d %3d %7d %7d %3d %7d %7d %s\n",
                               i_read,
                               i_group,
                               halo_id,
                               take_log10(properties.M_vir),
                               halo_tree_id,
                               halo_file_offset,
                               descendant_snap,
                               halo_index,
                               halo_descendant_id,
                               bridge_match_snap,
                               bridge_match_index,
                               halo_type,
                               halo_type_string);
                SID_free(SID_FARG halo_type_string);
             }
          }
          else{
             SID_fskip(sizeof(int),6,&fp_in_trees);
             SID_fskip(sizeof(int),4,&fp_in_bridge_match); // 3+1 'casue we're skipping n_subgroups_group too
          }
          SID_fread_all(&n_subgroups_group,sizeof(int),1,&fp_in_trees);
          if(mode==MATCH_SUBGROUPS){
             int j_subgroup;
             for(j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
                // Read subgroups
                SID_fread_all(&halo_id,           sizeof(int),1,&fp_in_trees);
                SID_fread_all(&halo_type,         sizeof(int),1,&fp_in_trees);
                SID_fread_all(&halo_descendant_id,sizeof(int),1,&fp_in_trees);
                SID_fread_all(&halo_tree_id,      sizeof(int),1,&fp_in_trees);
                SID_fread_all(&halo_file_offset,  sizeof(int),1,&fp_in_trees);
                SID_fread_all(&halo_index,        sizeof(int),1,&fp_in_trees);
                SID_fread_all(&bridge_match_id,   sizeof(int),1,&fp_in_bridge_match);
                SID_fread_all(&bridge_match_file, sizeof(int),1,&fp_in_bridge_match);
                SID_fread_all(&bridge_match_index,sizeof(int),1,&fp_in_bridge_match);
                // Write match
                int flag_is_halo=(i_read==i_read_in && i_halo==i_subgroup);
                if(trees->snap_list[bridge_match_file]==i_read_in && bridge_match_index==i_halo || flag_is_halo){
                   char *halo_type_string=NULL;
                   tree_case_flags_text(halo_type,"+",&halo_type_string);
                   fread_catalog_file(&fp_properties,NULL,&properties,NULL,i_subgroup);
                   int descendant_snap;
                   int bridge_match_snap;
                   if(halo_file_offset>0)
                      descendant_snap=trees->snap_list[i_file+halo_file_offset];
                   else
                      descendant_snap=-1;
                   if(bridge_match_file>=0)
                      bridge_match_snap=trees->snap_list[bridge_match_file];
                   else
                      bridge_match_snap=-1;
                   fprintf(fp_out,"%3d %7d %7d %5.2lf %7d %3d %3d %7d %7d %3d %7d %7d %s\n",
                                  i_read,
                                  i_subgroup,
                                  halo_id,
                                  take_log10(properties.M_vir),
                                  halo_tree_id,
                                  halo_file_offset,
                                  descendant_snap,
                                  halo_index,
                                  halo_descendant_id,
                                  bridge_match_snap,
                                  bridge_match_index,
                                  halo_type,
                                  halo_type_string);
                   SID_free(SID_FARG halo_type_string);
                }
             }
          }
          else{
             SID_fskip(sizeof(int),6*n_subgroups_group,&fp_in_trees);
             SID_fskip(sizeof(int),3*n_subgroups_group,&fp_in_bridge_match); 
          }
       }
       SID_fclose(&fp_in_trees);
       SID_fclose(&fp_in_bridge_match);
       fclose_catalog(&fp_properties);

       SID_log("Done.",SID_LOG_CLOSE);
    }
  }
  if(fp_out!=stderr)
    fclose(fp_out);

  SID_log("Output written to {%s}",SID_LOG_COMMENT,filename_out);

  // Clean-up
  free_trees(&trees);
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

