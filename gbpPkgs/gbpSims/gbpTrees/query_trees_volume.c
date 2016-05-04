#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int process_local(tree_info            *trees,
                  int                   i_pass,
                  char                 *filename_out_root,
                  double                radius2,
                  double                M_min,
                  int                   i_snap,
                  int                   i_halo,
                  int                   j_subgroup,
                  int                   i_group,
                  int                   halo_id,
                  int                   halo_file_offset,
                  int                   halo_type,
                  int                   halo_n_particles_peak,
                  int                   halo_tree_id,
                  int                   halo_index,
                  int                   halo_descendant_id,
                  int                   bridge_forematch_first_index,
                  int                   bridge_forematch_first_file,
                  int                   bridge_backmatch_index,
                  int                   bridge_backmatch_file,
                  int                   group_id,
                  halo_properties_info *properties,
                  halo_properties_info *properties_group,
                  int                  *n_list,
                  int                  *halo_list,
                  double                x_cen,
                  double                y_cen,
                  double                z_cen);
int process_local(tree_info            *trees,
                  int                   i_pass,
                  char                 *filename_out_root,
                  double                radius2,
                  double                M_min,
                  int                   i_snap,
                  int                   i_halo,
                  int                   j_subgroup,
                  int                   i_group,
                  int                   halo_id,
                  int                   halo_file_offset,
                  int                   halo_type,
                  int                   halo_n_particles_peak,
                  int                   halo_tree_id,
                  int                   halo_index,
                  int                   halo_descendant_id,
                  int                   bridge_forematch_first_index,
                  int                   bridge_forematch_first_file,
                  int                   bridge_backmatch_index,
                  int                   bridge_backmatch_file,
                  int                   group_id,
                  halo_properties_info *properties,
                  halo_properties_info *properties_group,
                  int                  *n_list,
                  int                  *halo_list,
                  double                x_cen,
                  double                y_cen,
                  double                z_cen){
   int flag_found=FALSE;
   int i_type=!(properties_group==NULL);
   // Perform search
   if(i_pass<2){
      double dx_i=d_periodic(((double)properties->position_MBP[0]-x_cen),trees->box_size);
      double dy_i=d_periodic(((double)properties->position_MBP[1]-y_cen),trees->box_size);
      double dz_i=d_periodic(((double)properties->position_MBP[2]-z_cen),trees->box_size);
      double r2_i=dx_i*dx_i+dy_i*dy_i+dz_i*dz_i;
      if(r2_i<radius2 && properties->M_vir>=M_min){
         flag_found=TRUE;
         if(i_pass==0)
            (*n_list)++;
         else{
            // Check if we have added this halo ID yet
            for(int i_scan=0;i_scan<(*n_list) && flag_found;i_scan++){
               if(halo_list[i_scan]==halo_id)
                  flag_found=FALSE;
            }
            // Add halo ID if it isn't in the list
            if(flag_found)
               halo_list[(*n_list)++]=halo_id;
         }
      }
   }
   // Perform write
   else if(i_pass==2){
      int flag_keep=FALSE;
      for(int j_list=0;j_list<(*n_list) && !flag_keep;j_list++) if(halo_id==halo_list[j_list]) flag_keep=TRUE;
      if(flag_keep){
         char filename_out[MAX_FILENAME_LENGTH];
         if(i_type==0)
            sprintf(filename_out,"%s_group_%09d.txt",filename_out_root,halo_id);
         else
            sprintf(filename_out,"%s_subgroup_%09d.txt",filename_out_root,halo_id);
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
         
         char *halo_type_string=NULL;
         tree_case_flags_text(halo_type,"+",&halo_type_string);
         fprintf(fp_out,"%le %7.3lf %3d %7d %7d %5.2lf %6d %6d %10.3le %10.3le %10.3le %10.3le %7d %2d %3d %7d %7d %3d %7d %3d %7d %7d %s",
                        trees->a_list[i_snap],
                        trees->z_list[i_snap],
                        trees->snap_list[i_snap],
                        i_halo,
                        halo_id,
                        take_log10(properties->M_vir),
                        properties->n_particles,
                        halo_n_particles_peak,
                        properties->position_MBP[0],
                        properties->position_MBP[1],
                        properties->position_MBP[2],
                        properties->R_vir,
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
         SID_free(SID_FARG halo_type_string);
         if(i_type==1)
            fprintf(fp_out," %7d %5d %7d %10.3le %10.3le %10.3le\n",
                           i_group,
                           group_id,
                           j_subgroup,
                           properties->position_MBP[0]-properties_group->position_MBP[0],
                           properties->position_MBP[1]-properties_group->position_MBP[1],
                           properties->position_MBP[2]-properties_group->position_MBP[2]);
         else
            fprintf(fp_out,"\n");
         // Write subgroup-specific stuff
         if(i_type==1){

         }
         fclose(fp_out);
      }
   }
   // Sanity check
   else
      SID_trap_error("Invalid mode passed to process_local.",ERROR_LOGIC);
}

int main(int argc, char *argv[]){
  int     n_search;
  int     n_files;
  int     k_read;
  int     max_n_groups;
  int     l_read;
  int    *n_particles_i;
  int    *n_particles_j;
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
  if(argc!=12)
     SID_trap_error("Invalid Syntax.",ERROR_SYNTAX);
  char    filename_SSimPL_root[MAX_FILENAME_LENGTH];
  char    filename_halos_root[MAX_FILENAME_LENGTH];
  char    filename_trees_root[MAX_FILENAME_LENGTH];
  char    filename_out_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_root,argv[1]);
  strcpy(filename_halos_root, argv[2]);
  strcpy(filename_trees_root, argv[3]);
  double x_cen   =(double)atof(argv[4]);
  double y_cen   =(double)atof(argv[5]);
  double z_cen   =(double)atof(argv[6]);
  double radius  =(double)atof(argv[7]);
  double z_min_in=(double)atof(argv[8]);
  double z_max_in=(double)atof(argv[9]);
  double M_min   =(double)atof(argv[10]);
  strcpy(filename_out_root,    argv[11]);
  double radius2=radius*radius;

  SID_log("Query trees for sphere (x,y,z,r)=(%.2lf,%.2lf,%.2lf,%.2lf) between z=%.2lf and z=%.2lf...",SID_LOG_OPEN,x_cen,y_cen,z_cen,radius,z_min_in,z_max_in);

  char                 filename_catalog_root[MAX_FILENAME_LENGTH];
  sprintf(filename_catalog_root,"%s/catalogs/%s",filename_SSimPL_root,filename_halos_root);

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

  // Perform query
  int  n_groups_list   =0;
  int  n_subgroups_list=0;
  int *group_list      =NULL;
  int *subgroup_list   =NULL;
  for(int i_pass=0;i_pass<3;i_pass++){
     if(i_pass==0)      SID_log("Counting halos to be queried...",SID_LOG_OPEN|SID_LOG_TIMER);
     else if(i_pass==1) SID_log("Identifying halos to be queried...",SID_LOG_OPEN|SID_LOG_TIMER);
     else if(i_pass==2) SID_log("Performing query...",SID_LOG_OPEN|SID_LOG_TIMER);
     // Write headers
     if(i_pass==2){
        for(int i_type=0;i_type<2;i_type++){
           int  n_list   =0;
           int *halo_list=NULL;
           if(i_type==0){
              n_list   =n_groups_list;
              halo_list=group_list;
           }
           else{
              n_list   =n_subgroups_list;
              halo_list=subgroup_list;
           }
           for(int i_list=0;i_list<n_list;i_list++){
              int   i_column=1;
              char  filename_out[MAX_FILENAME_LENGTH];
              if(i_type==0) sprintf(filename_out,"%s_group_%09d.txt",   filename_out_root,halo_list[i_list]);
              else          sprintf(filename_out,"%s_subgroup_%09d.txt",filename_out_root,halo_list[i_list]);
              FILE *fp_out=fopen(filename_out,"w");
              fprintf(fp_out,"# Column (%02d): Halo expansion factor\n",       i_column++);
              fprintf(fp_out,"#        (%02d): Halo redshift\n",               i_column++);
              fprintf(fp_out,"#        (%02d): Halo snapshot\n",               i_column++);
              fprintf(fp_out,"#        (%02d): Halo index\n",                  i_column++);
              fprintf(fp_out,"#        (%02d): Halo ID\n",                     i_column++);
              fprintf(fp_out,"#        (%02d): Halo log10(M_vir [M_sol/h])\n", i_column++);
              fprintf(fp_out,"#        (%02d): Halo n_particles\n",            i_column++);
              fprintf(fp_out,"#        (%02d): Halo n_particles_peak\n",       i_column++);
              fprintf(fp_out,"#        (%02d): Halo x [Mpc/h])\n",             i_column++);
              fprintf(fp_out,"#        (%02d): Halo y [Mpc/h])\n",             i_column++);
              fprintf(fp_out,"#        (%02d): Halo z [Mpc/h])\n",             i_column++);
              fprintf(fp_out,"#        (%02d): Halo radius [Mpc/h])\n",        i_column++);
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
              if(i_type==1){
                 fprintf(fp_out,"#        (%02d): Group index\n",              i_column++);
                 fprintf(fp_out,"#        (%02d): Group ID\n",                 i_column++);
                 fprintf(fp_out,"#        (%02d): Subgroup index\n",           i_column++);
                 fprintf(fp_out,"#        (%02d): FoF Centre dx [Mpc/h])\n",   i_column++);
                 fprintf(fp_out,"#        (%02d): FoF Centre dy [Mpc/h])\n",   i_column++);
                 fprintf(fp_out,"#        (%02d): FoF Centre dz [Mpc/h])\n",   i_column++);
              }
              fclose(fp_out);
           }
        }
     }

     // Set the snapshot range we need to scan
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

     // Loop over the range of snapshots
     int i_list=0;
     for(int i_snap=i_snap_start;i_snap<=i_snap_stop;i_snap++){
        // Get the snapshot
        int snapshot=trees->snap_list[i_snap];
        if(i_pass==2) SID_log("Processing snapshot %03d...",SID_LOG_OPEN,snapshot);

        // Open properties for this snapshot
        fp_catalog_info      fp_properties_groups;
        fp_catalog_info      fp_properties_subgroups;
        halo_properties_info properties_groups;
        halo_properties_info properties_subgroups;
        fopen_catalog(filename_catalog_root,
                      snapshot,
                      READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES,
                      &fp_properties_groups);
        fopen_catalog(filename_catalog_root,
                      snapshot,
                      READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES,
                      &fp_properties_subgroups);

        // Open horizontal trees for this snapshot
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

        // Read trees header
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

        // Scan through the trees
        int   i_subgroup=0;
        for(int i_group=0;i_group<n_groups;i_group++){
           // Read group
           int   group_id;
           int   group_type;
           int   group_descendant_id;
           int   group_tree_id;
           int   group_file_offset;
           int   group_index;
           int   group_n_particles_peak;
           int   n_subgroups_group;
           int   group_forematch_id;
           int   group_forematch_first_file;
           int   group_forematch_first_index;
           float group_forematch_first_score;
           int   group_forematch_default_file;
           int   group_forematch_default_index;
           float group_forematch_default_score;
           int   group_forematch_best_file;
           int   group_forematch_best_index;
           float group_forematch_best_score;
           float group_forematch_score_prog;
           int   group_backmatch_id;
           int   group_backmatch_file;
           int   group_backmatch_index;
           float group_backmatch_score;
           float group_backmatch_score_prog;
           SID_fread_all(&group_id,                     sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_type,                   sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_descendant_id,          sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_tree_id,                sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_file_offset,            sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_index,                  sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_n_particles_peak,       sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&n_subgroups_group,            sizeof(int),  1,&fp_in_trees);
           SID_fread_all(&group_forematch_id,           sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&group_forematch_first_file,   sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&group_forematch_first_index,  sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&group_forematch_first_score,  sizeof(float),1,&fp_in_bridge_forematch);
           SID_fread_all(&group_forematch_default_file, sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&group_forematch_default_index,sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&group_forematch_default_score,sizeof(float),1,&fp_in_bridge_forematch);
           SID_fread_all(&group_forematch_best_file,    sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&group_forematch_best_index,   sizeof(int),  1,&fp_in_bridge_forematch);
           SID_fread_all(&group_forematch_best_score,   sizeof(float),1,&fp_in_bridge_forematch);
           SID_fread_all(&group_forematch_score_prog,   sizeof(float),1,&fp_in_bridge_forematch);
           SID_fread_all(&group_backmatch_id,           sizeof(int),  1,&fp_in_bridge_backmatch);
           SID_fread_all(&group_backmatch_file,         sizeof(int),  1,&fp_in_bridge_backmatch);
           SID_fread_all(&group_backmatch_index,        sizeof(int),  1,&fp_in_bridge_backmatch);
           SID_fread_all(&group_backmatch_score,        sizeof(float),1,&fp_in_bridge_backmatch);
           SID_fread_all(&group_backmatch_score_prog,   sizeof(float),1,&fp_in_bridge_backmatch);
           SID_fskip(sizeof(int),1,&fp_in_bridge_forematch); // skip subhalo count 
           SID_fskip(sizeof(int),1,&fp_in_bridge_backmatch); // skip subhalo count 
           fread_catalog_file(&fp_properties_groups,NULL,NULL,&properties_groups,NULL,i_group);

           // Process group
           process_local(trees,
                         i_pass,
                         filename_out_root,
                         radius2,
                         M_min,
                         i_snap,
                         i_group,
                         0,
                         i_group,
                         group_id,
                         group_file_offset,
                         group_type,
                         group_n_particles_peak,
                         group_tree_id,
                         group_index,
                         group_descendant_id,
                         group_forematch_first_index,
                         group_forematch_first_file,
                         group_backmatch_index,
                         group_backmatch_file,
                         group_id,
                         &properties_groups,
                         NULL,
                         &n_groups_list,
                         group_list,
                         x_cen,
                         y_cen,
                         z_cen);

           for(int j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
              // Read subgroup
              int   subgroup_id;
              int   subgroup_type;
              int   subgroup_descendant_id;
              int   subgroup_tree_id;
              int   subgroup_file_offset;
              int   subgroup_index;
              int   subgroup_n_particles_peak;
              int   subgroup_forematch_id;
              int   subgroup_forematch_first_file;
              int   subgroup_forematch_first_index;
              float subgroup_forematch_first_score;
              int   subgroup_forematch_default_file;
              int   subgroup_forematch_default_index;
              float subgroup_forematch_default_score;
              int   subgroup_forematch_best_file;
              int   subgroup_forematch_best_index;
              float subgroup_forematch_best_score;
              float subgroup_forematch_score_prog;
              int   subgroup_backmatch_id;
              int   subgroup_backmatch_file;
              int   subgroup_backmatch_index;
              float subgroup_backmatch_score;
              float subgroup_backmatch_score_prog;
              SID_fread_all(&subgroup_id,                     sizeof(int),  1,&fp_in_trees);
              SID_fread_all(&subgroup_type,                   sizeof(int),  1,&fp_in_trees);
              SID_fread_all(&subgroup_descendant_id,          sizeof(int),  1,&fp_in_trees);
              SID_fread_all(&subgroup_tree_id,                sizeof(int),  1,&fp_in_trees);
              SID_fread_all(&subgroup_file_offset,            sizeof(int),  1,&fp_in_trees);
              SID_fread_all(&subgroup_index,                  sizeof(int),  1,&fp_in_trees);
              SID_fread_all(&subgroup_n_particles_peak,       sizeof(int),  1,&fp_in_trees);
              SID_fread_all(&subgroup_forematch_id,           sizeof(int),  1,&fp_in_bridge_forematch);
              SID_fread_all(&subgroup_forematch_first_file,   sizeof(int),  1,&fp_in_bridge_forematch);
              SID_fread_all(&subgroup_forematch_first_index,  sizeof(int),  1,&fp_in_bridge_forematch);
              SID_fread_all(&subgroup_forematch_first_score,  sizeof(float),1,&fp_in_bridge_forematch);
              SID_fread_all(&subgroup_forematch_default_file, sizeof(int),  1,&fp_in_bridge_forematch);
              SID_fread_all(&subgroup_forematch_default_index,sizeof(int),  1,&fp_in_bridge_forematch);
              SID_fread_all(&subgroup_forematch_default_score,sizeof(float),1,&fp_in_bridge_forematch);
              SID_fread_all(&subgroup_forematch_best_file,    sizeof(int),  1,&fp_in_bridge_forematch);
              SID_fread_all(&subgroup_forematch_best_index,   sizeof(int),  1,&fp_in_bridge_forematch);
              SID_fread_all(&subgroup_forematch_best_score,   sizeof(float),1,&fp_in_bridge_forematch);
              SID_fread_all(&subgroup_forematch_score_prog,   sizeof(float),1,&fp_in_bridge_forematch);
              SID_fread_all(&subgroup_backmatch_id,           sizeof(int),  1,&fp_in_bridge_backmatch);
              SID_fread_all(&subgroup_backmatch_file,         sizeof(int),  1,&fp_in_bridge_backmatch);
              SID_fread_all(&subgroup_backmatch_index,        sizeof(int),  1,&fp_in_bridge_backmatch);
              SID_fread_all(&subgroup_backmatch_score,        sizeof(float),1,&fp_in_bridge_backmatch);
              SID_fread_all(&subgroup_backmatch_score_prog,   sizeof(float),1,&fp_in_bridge_backmatch);
              fread_catalog_file(&fp_properties_subgroups,NULL,NULL,&properties_subgroups,NULL,i_subgroup);

              // Process subgroup
              process_local(trees,
                            i_pass,
                            filename_out_root,
                            radius2,
                            M_min,
                            i_snap,
                            i_subgroup,
                            j_subgroup,
                            i_group,
                            subgroup_id,
                            subgroup_file_offset,
                            subgroup_type,
                            subgroup_n_particles_peak,
                            subgroup_tree_id,
                            subgroup_index,
                            subgroup_descendant_id,
                            subgroup_forematch_first_index,
                            subgroup_forematch_first_file,
                            subgroup_backmatch_index,
                            subgroup_backmatch_file,
                            group_id,
                            &properties_subgroups,
                            &properties_groups,
                            &n_subgroups_list,
                            subgroup_list,
                            x_cen,
                            y_cen,
                            z_cen);
           }
        }
        fclose_catalog(&fp_properties_groups);
        fclose_catalog(&fp_properties_subgroups);
        SID_fclose(&fp_in_trees);
        SID_fclose(&fp_in_bridge_forematch);
        SID_fclose(&fp_in_bridge_backmatch);
        if(i_pass==2) 
           SID_log("Done.",SID_LOG_CLOSE);
     } // i_snap
     if(i_pass==0){
        group_list   =(int *)SID_malloc(sizeof(int)*n_groups_list);
        subgroup_list=(int *)SID_malloc(sizeof(int)*n_subgroups_list);
        n_groups_list   =0;
        n_subgroups_list=0;
     }
     else if(i_pass==1){
        SID_log("(%d groups and %d subgroups found)...",SID_LOG_CONTINUE,n_groups_list,n_subgroups_list);
        // Write list files
        for(int i_type=0;i_type<2;i_type++){
           char  filename_out[MAX_FILENAME_LENGTH];
           int   n_list   =0;
           int  *halo_list=NULL;
           if(i_type==0){
              sprintf(filename_out,"%s_group_list.txt",filename_out_root);
              n_list   =n_groups_list;
              halo_list=group_list;
           }
           else{
              sprintf(filename_out,"%s_subgroup_list.txt",filename_out_root);
              n_list   =n_subgroups_list;
              halo_list=subgroup_list;
           }
           FILE *fp_out=fopen(filename_out,"w");
           fprintf(fp_out,"# Halo IDs in list of tree tracks with base {%s}\n",filename_out_root);
           for(int i_list=0;i_list<n_list;i_list++)
              fprintf(fp_out,"%d\n",halo_list[i_list]);
           fclose(fp_out);
        }
     }
     SID_log("Done.",SID_LOG_CLOSE);
  } 

  // Clean-up
  SID_free(SID_FARG group_list);
  SID_free(SID_FARG subgroup_list);
  free_trees(&trees);
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

