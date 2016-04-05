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
  if(argc!=13)
     SID_trap_error("Invalid Syntax.",ERROR_SYNTAX);
  char    filename_SSimPL_root[MAX_FILENAME_LENGTH];
  char    filename_halos_root[MAX_FILENAME_LENGTH];
  char    filename_trees_root[MAX_FILENAME_LENGTH];
  char    filename_out_root[MAX_FILENAME_LENGTH];
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
  double x_cen   =(double)atof(argv[5]);
  double y_cen   =(double)atof(argv[6]);
  double z_cen   =(double)atof(argv[7]);
  double radius  =(double)atof(argv[8]);
  double z_min_in=(double)atof(argv[9]);
  double z_max_in=(double)atof(argv[10]);
  double M_min   =(double)atof(argv[11]);
  strcpy(filename_out_root,    argv[12]);
  double radius2=radius*radius;

  SID_log("Query trees for %s in sphere (x,y,z,r)=(%.2lf,%.2lf,%.2lf,%.2lf) between z=%.2lf and z=%.2lf...",SID_LOG_OPEN,argv[4],x_cen,y_cen,z_cen,radius,z_min_in,z_max_in);

  char                 filename_catalog_root[MAX_FILENAME_LENGTH];
  sprintf(filename_catalog_root,"%s/catalogs/%s",filename_SSimPL_root,filename_halos_root);

  int read_props_mode;
  if(mode==MATCH_GROUPS){
     strcat(filename_out_root,"_group");
     read_props_mode=READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES;
  }
  else{
     strcat(filename_out_root,"_subgroup");
     read_props_mode=READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES;
  }

  // Read tree header information
  tree_info *trees;
  char       filename_file_root[MAX_FILENAME_LENGTH];
  sprintf(filename_file_root,"%s/trees/%s",filename_SSimPL_root,filename_trees_root);
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
  init_trees_read(filename_SSimPL_root,filename_halos_root,filename_trees_root,TREE_READ_HEADER_ONLY,&trees);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

  // Turn given redshift range into snapshot range
  int i_snap_min_z=find_treesnap_z(trees,z_max_in);
  int i_snap_max_z=find_treesnap_z(trees,z_min_in);
  SID_log("z=%.2lf -> snapshot=%d",SID_LOG_COMMENT,z_min_in,trees->snap_list[i_snap_max_z]);
  SID_log("z=%.2lf -> snapshot=%d",SID_LOG_COMMENT,z_max_in,trees->snap_list[i_snap_min_z]);

  // Find the halos we want to query
  int  n_list=0;
  int *halo_list=NULL;
  for(int i_pass=0;i_pass<3;i_pass++){
     if(i_pass==0)      SID_log("Counting halos to be queried...",SID_LOG_OPEN|SID_LOG_TIMER);
     else if(i_pass==1) SID_log("Identifying halos to be queried...",SID_LOG_OPEN|SID_LOG_TIMER);
     else if(i_pass==2) SID_log("Performing query...",SID_LOG_OPEN|SID_LOG_TIMER);
     // Write headers
     if(i_pass==2){
        for(int i_list=0;i_list<n_list;i_list++){
           int   i_column=1;
           char  filename_out[MAX_FILENAME_LENGTH];
           sprintf(filename_out,"%s_%09d.txt",filename_out_root,halo_list[i_list]);
           FILE *fp_out=fopen(filename_out,"w");
           fprintf(fp_out,"# Column (%02d): Halo expansion factor\n",       i_column++);
           fprintf(fp_out,"#        (%02d): Halo redshift\n",               i_column++);
           fprintf(fp_out,"#        (%02d): Halo snapshot\n",               i_column++);
           fprintf(fp_out,"#        (%02d): Halo index\n",                  i_column++);
           fprintf(fp_out,"#        (%02d): Halo ID\n",                     i_column++);
           if(mode==MATCH_SUBGROUPS){
              fprintf(fp_out,"#        (%02d): Group index\n",               i_column++);
              fprintf(fp_out,"#        (%02d): Group ID\n",                  i_column++);
              fprintf(fp_out,"#        (%02d): Subgroup index\n",              i_column++);
           }
           fprintf(fp_out,"#        (%02d): Halo log10(M_vir [M_sol/h])\n", i_column++);
           fprintf(fp_out,"#        (%02d): Halo n_particles\n",            i_column++);
           fprintf(fp_out,"#        (%02d): Halo n_particles_peak\n",       i_column++);
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
           fclose(fp_out);
        }
     }

     int i_snap_start;
     int i_snap_stop;
     if(i_pass<2){
        i_snap_start=i_snap_min_z;
        i_snap_stop =i_snap_max_z;
     }
     else{
        i_snap_start=0;
        i_snap_stop =trees->n_snaps-1;
     }

     int i_list=0;
     for(int i_snap=i_snap_start;i_snap<=i_snap_stop;i_snap++){
        // Get the snapshot
        int snapshot=trees->snap_list[i_snap];
        if(i_pass==2) SID_log("Processing snapshot %03d...",SID_LOG_OPEN,snapshot);

        fp_catalog_info      fp_properties;
        halo_properties_info properties;
        fopen_catalog(filename_catalog_root,
                      snapshot,
                      read_props_mode,
                      &fp_properties);

        SID_fp fp_in_trees;
        SID_fp fp_in_bridge_forematch;
        SID_fp fp_in_bridge_backmatch;
        char filename_in[MAX_FILENAME_LENGTH];
        sprintf(filename_in,"%s/trees/%s/horizontal/trees/horizontal_trees_%03d.dat",filename_SSimPL_root,filename_trees_root,snapshot);
        SID_fopen(filename_in,"r",&fp_in_trees);
        sprintf(filename_in,"%s/trees/%s/horizontal/trees/horizontal_trees_forematch_pointers_%03d.dat",
                            filename_SSimPL_root,filename_trees_root,snapshot);
        SID_fopen(filename_in,"r",&fp_in_bridge_forematch);
        sprintf(filename_in,"%s/trees/%s/horizontal/trees/horizontal_trees_backmatch_pointers_%03d.dat",
                            filename_SSimPL_root,filename_trees_root,snapshot);
        SID_fopen(filename_in,"r",&fp_in_bridge_backmatch);

        int n_step_in;
        int n_search_in;
        int n_groups;
        int n_subgroups;
        int n_groups_max_in;
        int n_subgroups_max_in;
        int n_trees_subgroup_in;
        int n_trees_group_in;
        SID_fread_all(&n_step_in,          sizeof(int),1,&fp_in_trees);
        SID_fread_all(&n_search_in,        sizeof(int),1,&fp_in_trees);
        SID_fread_all(&n_groups,           sizeof(int),1,&fp_in_trees);
        SID_fread_all(&n_subgroups,        sizeof(int),1,&fp_in_trees);
        SID_fread_all(&n_groups_max_in,    sizeof(int),1,&fp_in_trees);
        SID_fread_all(&n_subgroups_max_in, sizeof(int),1,&fp_in_trees);
        SID_fread_all(&n_trees_subgroup_in,sizeof(int),1,&fp_in_trees);
        SID_fread_all(&n_trees_group_in,   sizeof(int),1,&fp_in_trees);
        SID_fskip(sizeof(int),8,&fp_in_bridge_forematch);
        SID_fskip(sizeof(int),8,&fp_in_bridge_backmatch);

        int   i_subgroup=0;
        char *halo_type_string=NULL;
        for(int i_group=0;i_group<n_groups;i_group++){
           // Read trees
           int   halo_id;
           int   halo_type;
           int   halo_descendant_id;
           int   halo_tree_id;
           int   halo_file_offset;
           int   halo_index;
           int   halo_n_particles_peak;
           int   group_id;
           int   group_type;
           int   group_descendant_id;
           int   group_tree_id;
           int   group_file_offset;
           int   group_index;
           int   group_n_particles_peak;
           int   n_subgroups_group;
           int   bridge_forematch_id;
           int   bridge_forematch_first_file;
           int   bridge_forematch_first_index;
           float bridge_forematch_first_score;
           int   bridge_forematch_default_file;
           int   bridge_forematch_default_index;
           float bridge_forematch_default_score;
           int   bridge_forematch_best_file;
           int   bridge_forematch_best_index;
           float bridge_forematch_best_score;
           float bridge_forematch_score_prog;
           int   bridge_backmatch_id;
           int   bridge_backmatch_file;
           int   bridge_backmatch_index;
           float bridge_backmatch_score;
           float bridge_backmatch_score_prog;
           SID_fread_all(&group_id,                      sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_type,                    sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_descendant_id,           sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_tree_id,                 sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_file_offset,             sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_index,                   sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_n_particles_peak,        sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&n_subgroups_group,             sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&bridge_forematch_id,           sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&bridge_forematch_first_file,   sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&bridge_forematch_first_index,  sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&bridge_forematch_first_score,  sizeof(float),1,&fp_in_bridge_forematch);
           SID_fread_all(&bridge_forematch_default_file, sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&bridge_forematch_default_index,sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&bridge_forematch_default_score,sizeof(float),1,&fp_in_bridge_forematch);
           SID_fread_all(&bridge_forematch_best_file,    sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&bridge_forematch_best_index,   sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&bridge_forematch_best_score,   sizeof(float),1,&fp_in_bridge_forematch);
           SID_fread_all(&bridge_forematch_score_prog,   sizeof(float),1,&fp_in_bridge_forematch);
           SID_fread_all(&bridge_backmatch_id,           sizeof(int),  1,&fp_in_bridge_backmatch);
           SID_fread_all(&bridge_backmatch_file,         sizeof(int),  1,&fp_in_bridge_backmatch);
           SID_fread_all(&bridge_backmatch_index,        sizeof(int),  1,&fp_in_bridge_backmatch);
           SID_fread_all(&bridge_backmatch_score,        sizeof(float),1,&fp_in_bridge_backmatch);
           SID_fread_all(&bridge_backmatch_score_prog,   sizeof(float),1,&fp_in_bridge_backmatch);
           SID_fskip(sizeof(int),1,&fp_in_bridge_forematch); 
           SID_fskip(sizeof(int),1,&fp_in_bridge_backmatch); 

           // Read group catalogs
           int flag_process=FALSE;
           if(mode==MATCH_GROUPS){
              fread_catalog_file(&fp_properties,NULL,NULL,&properties,NULL,i_group);
              flag_process=TRUE;
              i_halo=i_group;
           }

           for(int j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
              // Read subgroups
              if(mode==MATCH_SUBGROUPS){
                 SID_fread_all(&halo_id,                       sizeof(int),  1,&fp_in_trees);
                 SID_fread_all(&halo_type,                     sizeof(int),  1,&fp_in_trees);
                 SID_fread_all(&halo_descendant_id,            sizeof(int),  1,&fp_in_trees);
                 SID_fread_all(&halo_tree_id,                  sizeof(int),  1,&fp_in_trees);
                 SID_fread_all(&halo_file_offset,              sizeof(int),  1,&fp_in_trees);
                 SID_fread_all(&halo_index,                    sizeof(int),  1,&fp_in_trees);
                 SID_fread_all(&halo_n_particles_peak,         sizeof(int),  1,&fp_in_trees);
                 SID_fread_all(&bridge_forematch_id,           sizeof(int),  1,&fp_in_bridge_forematch);
                 SID_fread_all(&bridge_forematch_first_file,   sizeof(int),  1,&fp_in_bridge_forematch);
                 SID_fread_all(&bridge_forematch_first_index,  sizeof(int),  1,&fp_in_bridge_forematch);
                 SID_fread_all(&bridge_forematch_first_score,  sizeof(float),1,&fp_in_bridge_forematch);
                 SID_fread_all(&bridge_forematch_default_file, sizeof(int),  1,&fp_in_bridge_forematch);
                 SID_fread_all(&bridge_forematch_default_index,sizeof(int),  1,&fp_in_bridge_forematch);
                 SID_fread_all(&bridge_forematch_default_score,sizeof(float),1,&fp_in_bridge_forematch);
                 SID_fread_all(&bridge_forematch_best_file,    sizeof(int),  1,&fp_in_bridge_forematch);
                 SID_fread_all(&bridge_forematch_best_index,   sizeof(int),  1,&fp_in_bridge_forematch);
                 SID_fread_all(&bridge_forematch_best_score,   sizeof(float),1,&fp_in_bridge_forematch);
                 SID_fread_all(&bridge_forematch_score_prog,   sizeof(float),1,&fp_in_bridge_forematch);
                 SID_fread_all(&bridge_backmatch_id,           sizeof(int),  1,&fp_in_bridge_backmatch);
                 SID_fread_all(&bridge_backmatch_file,         sizeof(int),  1,&fp_in_bridge_backmatch);
                 SID_fread_all(&bridge_backmatch_index,        sizeof(int),  1,&fp_in_bridge_backmatch);
                 SID_fread_all(&bridge_backmatch_score,        sizeof(float),1,&fp_in_bridge_backmatch);
                 SID_fread_all(&bridge_backmatch_score_prog,   sizeof(float),1,&fp_in_bridge_backmatch);
                 fread_catalog_file(&fp_properties,NULL,NULL,&properties,NULL,i_subgroup);
                 i_halo=i_subgroup;
                 flag_process=TRUE;
              }
              else{
                 halo_id              =group_id;
                 halo_type            =group_type;
                 halo_descendant_id   =group_descendant_id;
                 halo_tree_id         =group_tree_id;
                 halo_file_offset     =group_file_offset;
                 halo_index           =group_index;
                 halo_n_particles_peak=group_n_particles_peak;
                 SID_fskip(sizeof(int),  7,&fp_in_trees); 
                 SID_fskip(sizeof(int),  7,&fp_in_bridge_forematch); 
                 SID_fskip(sizeof(int),  3,&fp_in_bridge_backmatch); 
                 SID_fskip(sizeof(float),4,&fp_in_bridge_forematch); 
                 SID_fskip(sizeof(float),2,&fp_in_bridge_backmatch); 
              }

              if(flag_process){
                 if(i_pass==2){
                    int flag_keep=FALSE;
                    for(int j_list=0;j_list<n_list && !flag_keep;j_list++) if(halo_id==halo_list[j_list]) flag_keep=TRUE;
                    if(flag_keep){
                       char filename_out[MAX_FILENAME_LENGTH];
                       sprintf(filename_out,"%s_%09d.txt",filename_out_root,halo_id);
                       FILE *fp_out=fopen(filename_out,"a");
                       int descendant_snap      =-1;
                       int bridge_forematch_first_snap=-1;
                       int bridge_backmatch_snap=-1;
                       if(halo_file_offset>0)
                          descendant_snap=trees->snap_list[i_snap+halo_file_offset];
                       if(bridge_forematch_first_file>=0)
                          bridge_forematch_first_snap=trees->snap_list[bridge_forematch_first_file];
                       if(bridge_backmatch_file>=0)
                          bridge_backmatch_snap=trees->snap_list[bridge_backmatch_file];
                       tree_case_flags_text(halo_type,"+",&halo_type_string);
                       if(mode==MATCH_SUBGROUPS)
                          fprintf(fp_out,"%le %7.3lf %3d %7d %7d %7d %7d %4d %5.2lf %6d %6d %5.2lf %5.2lf %5.2lf %7d %3d %3d %7d %7d %3d %7d %3d %7d %7d %s\n",
                                         trees->a_list[i_snap],
                                         trees->z_list[i_snap],
                                         trees->snap_list[i_snap],
                                         i_halo,
                                         halo_id,
                                         i_group,
                                         group_id,
                                         j_subgroup,
                                         take_log10(properties.M_vir),
                                         properties.n_particles,
                                         halo_n_particles_peak,
                                         properties.position_MBP[0],
                                         properties.position_MBP[1],
                                         properties.position_MBP[2],
                                         halo_tree_id,
                                         halo_file_offset,
                                         descendant_snap,
                                         halo_index,
                                         halo_descendant_id,
                                         bridge_forematch_first_snap,
                                         bridge_forematch_first_index,
                                         bridge_backmatch_snap,
                                         bridge_backmatch_index,
                                         halo_type,
                                         halo_type_string);
                       else
                          fprintf(fp_out,"%le %7.3lf %3d %7d %7d %5.2lf %6d %6d %5.2lf %5.2lf %5.2lf %7d %3d %3d %7d %7d %3d %7d %3d %7d %7d %s\n",
                                         trees->a_list[i_snap],
                                         trees->z_list[i_snap],
                                         trees->snap_list[i_snap],
                                         i_halo,
                                         halo_id,
                                         take_log10(properties.M_vir),
                                         properties.n_particles,
                                         halo_n_particles_peak,
                                         properties.position_MBP[0],
                                         properties.position_MBP[1],
                                         properties.position_MBP[2],
                                         halo_tree_id,
                                         halo_file_offset,
                                         descendant_snap,
                                         halo_index,
                                         halo_descendant_id,
                                         bridge_forematch_first_snap,
                                         bridge_forematch_first_index,
                                         bridge_backmatch_snap,
                                         bridge_backmatch_index,
                                         halo_type,
                                         halo_type_string);
                       fclose(fp_out);
                    }
                 }
                 else{
                    double dx_i=d_periodic(((double)properties.position_MBP[0]-x_cen),trees->box_size);
                    double dy_i=d_periodic(((double)properties.position_MBP[1]-y_cen),trees->box_size);
                    double dz_i=d_periodic(((double)properties.position_MBP[2]-z_cen),trees->box_size);
                    double r2_i=dx_i*dx_i+dy_i*dy_i+dz_i*dz_i;
                    if(r2_i<radius2 && properties.M_vir>=M_min){
                       if(i_pass==0)
                          n_list++;
                       else{
                          // Check if we have added this halo ID yet
                          int flag_continue=TRUE;
                          for(int i_scan=0;i_scan<i_list && flag_continue;i_scan++){
                             if(halo_list[i_scan]==halo_id)
                                flag_continue=FALSE;
                          }
                          // Add halo ID uf it isn't in the list
                          if(flag_continue)
                             halo_list[i_list++]=halo_id;
                       }
                    }
                 }
              }
              // This is needed so we print results for groups only once
              flag_process=FALSE;
           }
        }
        SID_free(SID_FARG halo_type_string);
        fclose_catalog(&fp_properties);
        SID_fclose(&fp_in_trees);
        SID_fclose(&fp_in_bridge_forematch);
        SID_fclose(&fp_in_bridge_backmatch);
        if(i_pass==2) 
           SID_log("Done.",SID_LOG_CLOSE);
     } // i_snap
     if(i_pass==0)
        halo_list=(int *)SID_malloc(sizeof(int)*n_list);
     else if(i_pass==1){
        n_list=i_list;
        SID_log("(%d found)...",SID_LOG_CONTINUE,n_list);
     }
     SID_log("Done.",SID_LOG_CLOSE);
  } 

  // Clean-up
  SID_free(SID_FARG halo_list);
  free_trees(&trees);
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

