#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){
  int     n_search;
  int     i_halo;
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
  int     i_read_start;
  int     i_read_stop;
  SID_fp  fp_in;

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  if(argc!=7 && argc!=6)
     SID_trap_error("Invalid Syntax.",ERROR_SYNTAX);
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
  int halo_id_find;
  int find_mode;
  char filename_out[MAX_FILENAME_LENGTH];
  if(argc==6){
     halo_id_find=atoi(argv[5]);
     i_read      =0;
     i_halo      =0;
     find_mode   =0;
     SID_log("Querying trees for halo ID #%d from {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,halo_id_find,filename_trees_root);
     sprintf(filename_out,"%s_%d.txt",filename_trees_root,halo_id_find);
  }
  else{
     i_read   =atoi(argv[5]);
     i_halo   =atoi(argv[6]);
     find_mode=1;
     SID_log("Querying trees for halo #%d in file #%d from {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,i_halo,i_read,filename_trees_root);
     sprintf(filename_out,"%s_%d_%d.txt",filename_trees_root,i_read,i_halo);
  }

  // Read tree header information
  tree_info *trees;
  char       filename_file_root[MAX_FILENAME_LENGTH];
  sprintf(filename_file_root,"%s/trees/%s",filename_SSimPL_root,filename_trees_root);
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
  init_trees_read(filename_SSimPL_root,filename_trees_root,TREE_READ_HEADER_ONLY,&trees);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

  // If we're searching for an ID, set starting snapshot to be the first ...
  int i_file=0;
  if(find_mode==0)
     i_read=trees->snap_list[0];
  // ... else check that the given snapshot is valid
  else{
     int j_read=trees->snap_list[i_file];
     while(i_file<(trees->n_snaps-1) && j_read<i_read){
       i_file++;
       j_read=trees->snap_list[i_file];
     }
     if(i_read!=trees->snap_list[i_file])
        SID_trap_error("Invalid snapshot specified {%d}.",ERROR_LOGIC,i_read);
  }

  // Open output file
  FILE *fp_out=fopen(filename_out,"w");

  int halo_file_offset=1;
  int flag_write_header=TRUE;
  for(;
      i_read<=trees->i_read_stop && ((find_mode==1 && i_halo>=0) || find_mode==0);
      i_read+=(halo_file_offset*trees->i_read_step),i_file+=halo_file_offset){
    SID_log("Processing snapshot %03d (%03d of %03d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read,i_file+1,trees->n_snaps);
    int halo_id               =0;
    int halo_type             =0;
    int halo_descendant_id    =0;
    int halo_tree_id          =0;
    int halo_index            =0;
    int bridge_forematch_id   =0;
    int bridge_forematch_file =0;
    int bridge_forematch_index=0;
    int bridge_backmatch_id   =0;
    int bridge_backmatch_file =0;
    int bridge_backmatch_index=0;

    // Read tree entry
    char filename_in[MAX_FILENAME_LENGTH];
    SID_fp fp_in_trees;
    SID_fp fp_in_bridge_forematch;
    SID_fp fp_in_bridge_backmatch;
    sprintf(filename_in,"%s/trees/%s/horizontal/trees/horizontal_trees_%03d.dat",filename_SSimPL_root,filename_trees_root,i_read);
    SID_fopen(filename_in,"r",&fp_in_trees);
    sprintf(filename_in,"%s/trees/%s/horizontal/trees/horizontal_trees_bridge_forematch_pointers_%03d.dat",
                        filename_SSimPL_root,filename_trees_root,i_read);
    SID_fopen(filename_in,"r",&fp_in_bridge_forematch);
    sprintf(filename_in,"%s/trees/%s/horizontal/trees/horizontal_trees_bridge_backmatch_pointers_%03d.dat",
                        filename_SSimPL_root,filename_trees_root,i_read);
    SID_fopen(filename_in,"r",&fp_in_bridge_backmatch);
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
    if(find_mode==1){
       if(mode==MATCH_GROUPS         && i_halo>=n_groups)    SID_trap_error("Invalid group halo index (ie. %d>=%d).",   ERROR_LOGIC,i_halo,n_groups);
       else if(mode==MATCH_SUBGROUPS && i_halo>=n_subgroups) SID_trap_error("Invalid subgroup halo index (ie. %d>=%d).",ERROR_LOGIC,i_halo,n_subgroups);
    }
    SID_fskip(sizeof(int),8,&fp_in_bridge_forematch);
    SID_fskip(sizeof(int),8,&fp_in_bridge_backmatch);
    int flag_done=FALSE;
    int i_group;
    int i_subgroup;
    for(i_group=0,i_subgroup=0;i_group<n_groups && !flag_done;i_group++){
       // Read groups
       int n_subgroups_group;
       if(mode==MATCH_GROUPS){
          SID_fread_all(&halo_id,               sizeof(int),1,&fp_in_trees);
          SID_fread_all(&halo_type,             sizeof(int),1,&fp_in_trees);
          SID_fread_all(&halo_descendant_id,    sizeof(int),1,&fp_in_trees);
          SID_fread_all(&halo_tree_id,          sizeof(int),1,&fp_in_trees);
          SID_fread_all(&halo_file_offset,      sizeof(int),1,&fp_in_trees);
          SID_fread_all(&halo_index,            sizeof(int),1,&fp_in_trees);
          SID_fread_all(&bridge_forematch_id,   sizeof(int),1,&fp_in_bridge_forematch);
          SID_fread_all(&bridge_forematch_file, sizeof(int),1,&fp_in_bridge_forematch);
          SID_fread_all(&bridge_forematch_index,sizeof(int),1,&fp_in_bridge_forematch);
          SID_fread_all(&bridge_backmatch_id,   sizeof(int),1,&fp_in_bridge_backmatch);
          SID_fread_all(&bridge_backmatch_file, sizeof(int),1,&fp_in_bridge_backmatch);
          SID_fread_all(&bridge_backmatch_index,sizeof(int),1,&fp_in_bridge_backmatch);
          if((find_mode==0 && halo_id==halo_id_find) || (find_mode==1 && i_halo==i_group)){
             if(find_mode==0) i_halo=i_group;
             flag_done=TRUE;
          }
       }
       else{
          SID_fskip(sizeof(int),6,&fp_in_trees);
          SID_fskip(sizeof(int),4,&fp_in_bridge_forematch); // 3+1 'casue we're skipping n_subgroups_group too
          SID_fskip(sizeof(int),4,&fp_in_bridge_backmatch); // 3+1 'casue we're skipping n_subgroups_group too
       }
       SID_fread_all(&n_subgroups_group,sizeof(int),1,&fp_in_trees);
       if(mode==MATCH_SUBGROUPS){
          int j_subgroup;
          for(j_subgroup=0;j_subgroup<n_subgroups_group && !flag_done;i_subgroup++,j_subgroup++){
             // Read subgroups
             SID_fread_all(&halo_id,               sizeof(int),1,&fp_in_trees);
             SID_fread_all(&halo_type,             sizeof(int),1,&fp_in_trees);
             SID_fread_all(&halo_descendant_id,    sizeof(int),1,&fp_in_trees);
             SID_fread_all(&halo_tree_id,          sizeof(int),1,&fp_in_trees);
             SID_fread_all(&halo_file_offset,      sizeof(int),1,&fp_in_trees);
             SID_fread_all(&halo_index,            sizeof(int),1,&fp_in_trees);
             SID_fread_all(&bridge_forematch_id,   sizeof(int),1,&fp_in_bridge_forematch);
             SID_fread_all(&bridge_forematch_file, sizeof(int),1,&fp_in_bridge_forematch);
             SID_fread_all(&bridge_forematch_index,sizeof(int),1,&fp_in_bridge_forematch);
             SID_fread_all(&bridge_backmatch_id,   sizeof(int),1,&fp_in_bridge_backmatch);
             SID_fread_all(&bridge_backmatch_file, sizeof(int),1,&fp_in_bridge_backmatch);
             SID_fread_all(&bridge_backmatch_index,sizeof(int),1,&fp_in_bridge_backmatch);
             if((find_mode==0 && halo_id==halo_id_find) || (find_mode==1 && i_halo==i_subgroup)){
                if(find_mode==0) i_halo=i_subgroup;
                flag_done=TRUE;
             }
          }
       }
       else{
          SID_fskip(sizeof(int),6*n_subgroups_group,&fp_in_trees);
          SID_fskip(sizeof(int),3*n_subgroups_group,&fp_in_bridge_forematch); 
          SID_fskip(sizeof(int),3*n_subgroups_group,&fp_in_bridge_backmatch); 
       }
    }
    SID_fclose(&fp_in_trees);
    SID_fclose(&fp_in_bridge_forematch);
    SID_fclose(&fp_in_bridge_backmatch);
    if(find_mode==1 && !flag_done)
       SID_trap_error("Halo {%d} could not be processed.",ERROR_LOGIC,i_halo);

    if(flag_done){
       // Read properties
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
       fread_catalog_file(&fp_properties,NULL,&properties,NULL,i_halo);
       fclose_catalog(&fp_properties);

       // Write results
       char *halo_type_string=NULL;
       tree_case_flags_text(halo_type,"+",&halo_type_string);
       if(flag_write_header){
          int i_column=1;
          fprintf(fp_out,"# Column (%02d): Halo snapshot\n",               i_column++);
          fprintf(fp_out,"#        (%02d): Halo expansion factor\n",       i_column++);
          fprintf(fp_out,"#        (%02d): Halo redshift\n",               i_column++);
          fprintf(fp_out,"#        (%02d): Halo index\n",                  i_column++);
          fprintf(fp_out,"#        (%02d): Halo ID\n",                     i_column++);
          fprintf(fp_out,"#        (%02d): Halo log10(M_vir [M_sol/h])\n", i_column++);
          fprintf(fp_out,"#        (%02d): Halo x [Mpc/h])\n",             i_column++);
          fprintf(fp_out,"#        (%02d): Halo y [Mpc/h])\n",             i_column++);
          fprintf(fp_out,"#        (%02d): Halo z [Mpc/h])\n",             i_column++);
          fprintf(fp_out,"#        (%02d): Halo tree ID\n",                i_column++);
          fprintf(fp_out,"#        (%02d): Descendant file offset\n",      i_column++);
          fprintf(fp_out,"#        (%02d): Descendant snapshot\n",         i_column++);
          fprintf(fp_out,"#        (%02d): Descendant index\n",            i_column++);
          fprintf(fp_out,"#        (%02d): Descendant ID\n",               i_column++);
          fprintf(fp_out,"#        (%02d): Bridge forematch snapshot\n",   i_column++);
          fprintf(fp_out,"#        (%02d): Bridge forematch index\n",      i_column++);
          fprintf(fp_out,"#        (%02d): Bridge backmatch snapshot\n",   i_column++);
          fprintf(fp_out,"#        (%02d): Bridge backmatch index\n",      i_column++);
          fprintf(fp_out,"#        (%02d): Halo type\n",                   i_column++);
          fprintf(fp_out,"#        (%02d): Halo type string\n",            i_column++);
          flag_write_header=FALSE;
       }
       int descendant_snap;
       int bridge_forematch_snap;
       int bridge_backmatch_snap;
       if(halo_file_offset>0)
          descendant_snap=trees->snap_list[i_file+halo_file_offset];
       else
          descendant_snap=-1;
       if(bridge_forematch_file>=0)
          bridge_forematch_snap=trees->snap_list[bridge_forematch_file];
       else
          bridge_forematch_snap=-1;
       if(bridge_backmatch_file>=0)
          bridge_backmatch_snap=trees->snap_list[bridge_backmatch_file];
       else
          bridge_backmatch_snap=-1;
       fprintf(fp_out,"%3d %le %7.3lf %7d %7d %5.2lf %5.2lf %5.2lf %5.2lf %7d %3d %3d %7d %7d %3d %7d %3d %7d %7d %s\n",
                      i_read,
                      a_of_z(trees->z_list[i_read]),
                      trees->z_list[i_read],
                      i_halo,
                      halo_id,
                      take_log10(properties.M_vir),
                      properties.position_COM[0],
                      properties.position_COM[1],
                      properties.position_COM[2],
                      halo_tree_id,
                      halo_file_offset,
                      descendant_snap,
                      halo_index,
                      halo_descendant_id,
                      bridge_forematch_snap,
                      bridge_forematch_index,
                      bridge_backmatch_snap,
                      bridge_backmatch_index,
                      halo_type,
                      halo_type_string);
       SID_free(SID_FARG halo_type_string);

       // Set descendant
       i_halo=halo_index;
    }

    // This ends the loop if we've reached the end
    //   of a descendant line (where halo_file_offset is
    //   set to -1 and we end-up in an infinite loop
    //   if we don't do something)    
    if(find_mode==0){
       if(flag_done && halo_file_offset<0)
          i_read=trees->snap_list[trees->n_snaps-1];
       halo_file_offset=1;
    }

    SID_log("Done.",SID_LOG_CLOSE);
  }
  if(fp_out!=stderr)
    fclose(fp_out);

  SID_log("Output written to {%s}",SID_LOG_COMMENT,filename_out);

  // Clean-up
  free_trees(&trees);
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

