#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  char filename_SSimPL_dir[MAX_FILENAME_LENGTH];
  char filename_halo_version_root[MAX_FILENAME_LENGTH];
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_trees_name[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  char filename_out_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_dir,       argv[1]);
  strcpy(filename_halo_version_root,argv[2]);
  strcpy(filename_trees_name,       argv[3]);
  double z_lo         =(double)atof(argv[4]);
  double z_hi         =(double)atof(argv[5]);
  double M_cut_min    =(double)atof(argv[6]);
  strcpy(filename_out_root,         argv[7]);

  // Set some filenames
  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_dir,filename_trees_name);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_dir,filename_halo_version_root);

  SID_log("Creating a catalog matched across redshifts z_lo~%lf and z_hi~%lf...",SID_LOG_OPEN|SID_LOG_TIMER,z_lo,z_hi);

  // Perform analysis
  tree_info *trees;
  read_trees(filename_SSimPL_dir,
             filename_halo_version_root,
             filename_trees_name,
             TREE_MODE_DEFAULT|TREE_READ_EXTENDED_POINTERS,
             &trees);

  // Read ancillary data
  read_trees_catalogs(trees,
                      filename_SSimPL_dir,
                      filename_halo_version_root,
                      READ_TREES_CATALOGS_GROUPS|READ_TREES_CATALOGS_SUBGROUPS);

  // Determine which snapshots best match the given redshifts
  int i_z_lo=find_treesnap_z(trees,z_lo);
  int i_z_hi=find_treesnap_z(trees,z_hi);

  // Generate two catalogs
  SID_log("Compiling catalogs...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(int i_cat=0;i_cat<2;i_cat++){
     char filename_out[MAX_FILENAME_LENGTH];
     switch(i_cat){
        case 0:
           SID_log("Processing group catalog...",SID_LOG_OPEN);
           sprintf(filename_out,"%s_groups.txt",filename_out_root);
           break;
        case 1:
           SID_log("Processing subgroup catalog...",SID_LOG_OPEN);
           sprintf(filename_out,"%s_subgroups.txt",filename_out_root);
           break;
     }

     // Open file and write header
     FILE *fp_out  =fopen(filename_out,"w");
     int   i_column=1;
     fprintf(fp_out,"# Catalog matched from z_hi=%lf to z_lo=%lf\n",trees->z_list[i_z_hi],trees->z_list[i_z_lo]);
     fprintf(fp_out,"#   SSiMPL_dir  ={%s}\n",filename_SSimPL_dir);
     fprintf(fp_out,"#   halo_version={%s}\n",filename_halo_version_root);
     fprintf(fp_out,"#   tree_version={%s}\n",filename_trees_name);
     fprintf(fp_out,"#\n");
     if(i_cat==0){
        fprintf(fp_out,"# Column (%02d): FoF      index   (z_hi)\n",i_column++);
        fprintf(fp_out,"#        (%02d): FoF      index   (z_lo)\n",i_column++);
     }
     else{
        fprintf(fp_out,"# Column (%02d): subgroup index   (z_hi)\n",i_column++);
        fprintf(fp_out,"#        (%02d): subgroup index   (z_lo)\n",i_column++);
        fprintf(fp_out,"#        (%02d): FoF      index   (z_hi)\n",i_column++);
        fprintf(fp_out,"#        (%02d): FoF      index   (z_lo)\n",i_column++);
     }
     fprintf(fp_out,"#        (%02d): n_particles      (z_hi)\n",i_column++);
     if(i_cat==0){
        fprintf(fp_out,"#        (%02d): n_subgroups      (z_hi)\n",i_column++);
        fprintf(fp_out,"#        (%02d): n_subgroups      (z_lo)\n",i_column++);
        fprintf(fp_out,"#        (%02d): sig_v_1D  [km/s] (z_hi)\n",i_column++);
        fprintf(fp_out,"#        (%02d): sig_v_1D  [km/s] (z_lo)\n",i_column++);
        fprintf(fp_out,"#        (%02d): sig_v_1Dp [km/s] (z_hi)\n",i_column++);
        fprintf(fp_out,"#        (%02d): sig_v_1Dp [km/s] (z_lo)\n",i_column++);
     }
     else{
        fprintf(fp_out,"#        (%02d): M_vir_sub      (z_hi)\n",i_column++);
        fprintf(fp_out,"#        (%02d): M_vir_sub      (z_lo)\n",i_column++);
     }
     fprintf(fp_out,"#        (%02d): M_vir_FoF        (z_hi)\n",i_column++);
     fprintf(fp_out,"#        (%02d): M_vir_FoF        (z_lo)\n",i_column++);
     fprintf(fp_out,"#        (%02d): x                (z_hi)\n",i_column++);
     fprintf(fp_out,"#        (%02d): y                (z_hi)\n",i_column++);
     fprintf(fp_out,"#        (%02d): z                (z_hi)\n",i_column++);
     fprintf(fp_out,"#        (%02d): x                (z_lo)\n",i_column++);
     fprintf(fp_out,"#        (%02d): y                (z_lo)\n",i_column++);
     fprintf(fp_out,"#        (%02d): z                (z_lo)\n",i_column++);
     fprintf(fp_out,"#        (%02d): v_x              (z_hi)\n",i_column++);
     fprintf(fp_out,"#        (%02d): v_y              (z_hi)\n",i_column++);
     fprintf(fp_out,"#        (%02d): v_z              (z_hi)\n",i_column++);
     fprintf(fp_out,"#        (%02d): v_x              (z_lo)\n",i_column++);
     fprintf(fp_out,"#        (%02d): v_y              (z_lo)\n",i_column++);
     fprintf(fp_out,"#        (%02d): v_z              (z_lo)\n",i_column++);

     // Write catalog
     tree_node_info        *current;
     halo_properties_info **group_properties   =(halo_properties_info **)ADaPS_fetch(trees->data,"properties_groups");
     halo_properties_info **subgroup_properties=(halo_properties_info **)ADaPS_fetch(trees->data,"properties_subgroups");
     if(i_cat==0)
        current=trees->first_neighbour_groups[i_z_hi];
     else
        current=trees->first_neighbour_subgroups[i_z_hi];
     while(current!=NULL){
        tree_node_info       *current_subgroup;
        tree_node_info       *current_group;
        halo_properties_info *current_properties;
        halo_properties_info *current_group_properties;
        halo_properties_info *current_subgroup_properties;
        if(i_cat==0){
           current_subgroup           =NULL;
           current_group              =current;
           current_group_properties   =&(group_properties[current_group->snap_tree][current_group->neighbour_index]);
           current_subgroup_properties=&(subgroup_properties[current->snap_tree][current->neighbour_index]);
           current_properties         =current_group_properties;
        }
        else{
           current_subgroup           =current;
           current_group              =current->parent;
           current_group_properties   =&(group_properties[current_group->snap_tree][current_group->neighbour_index]);
           current_subgroup_properties=&(subgroup_properties[current->snap_tree][current->neighbour_index]);
           current_properties         =current_subgroup_properties;
        }
        if(current_properties->M_vir>=M_cut_min){
           tree_node_info *matched;
           int             flag_exact=find_treenode_snap_equals_given(trees,current,i_z_lo,&matched,TREE_PROGENITOR_ORDER_N_PARTICLES_PEAK);
           if(matched!=NULL && flag_exact){
              int n_sub_lo;
              int n_sub_hi;
              tree_node_info *matched_group;
              tree_node_info *matched_subgroup;
              halo_properties_info *matched_group_properties;
              halo_properties_info *matched_subgroup_properties;
              double current_group_sigma_v=0.;
              double matched_group_sigma_v=0.;
              if(i_cat==0){
                 double v_mean;
                 matched_subgroup           =NULL;
                 matched_group              =matched;
                 matched_subgroup_properties=NULL;
                 matched_group_properties   =&(group_properties[matched_group->snap_tree][matched_group->neighbour_index]);
                 tree_node_info *current_j;

                 // Compute 1D velocity dispersion for current group                
                 current_j=current_group->substructure_first;
                 v_mean   =0.;
                 n_sub_hi =0;
                 while(current_j!=NULL){
                    double M_i=subgroup_properties[current_j->snap_tree][current_j->neighbour_index].M_vir;
                    if(M_i>M_cut_min){
                       float  v_i=subgroup_properties[current_j->snap_tree][current_j->neighbour_index].velocity_COM[0];
                       v_mean+=v_i;
                       n_sub_hi++;
                    }
                    current_j=current_j->substructure_next;
                 }
                 current_j=current_group->substructure_first;
                 current_group_sigma_v=0.;
                 if(n_sub_hi>1){
                    v_mean/=(double)n_sub_hi;
                    while(current_j!=NULL){
                       double M_i=subgroup_properties[current_j->snap_tree][current_j->neighbour_index].M_vir;
                       if(M_i>M_cut_min){
                          float v_i=subgroup_properties[current_j->snap_tree][current_j->neighbour_index].velocity_COM[0];
                          current_group_sigma_v+=pow(v_i-v_mean,2.);
                       }
                       current_j=current_j->substructure_next;
                    }
                    current_group_sigma_v=sqrt(current_group_sigma_v/(double)(n_sub_hi-1));
                 }

                 // Compute 1D velocity dispersion for matched group                
                 current_j=matched_group->substructure_first;
                 v_mean   =0.;
                 n_sub_lo =0;
                 while(current_j!=NULL){
                    double M_i=subgroup_properties[current_j->snap_tree][current_j->neighbour_index].M_vir;
                    if(M_i>M_cut_min){
                       float v_i=subgroup_properties[current_j->snap_tree][current_j->neighbour_index].velocity_COM[0];
                       v_mean+=v_i;
                       n_sub_lo++;
                    }
                    current_j=current_j->substructure_next;
                 }
                 current_j=matched_group->substructure_first;
                 matched_group_sigma_v=0.;
                 if(n_sub_lo>1){
                    v_mean/=(double)n_sub_lo;
                    while(current_j!=NULL){
                       double M_i=subgroup_properties[current_j->snap_tree][current_j->neighbour_index].M_vir;
                       if(M_i>M_cut_min){
                          float v_i=subgroup_properties[current_j->snap_tree][current_j->neighbour_index].velocity_COM[0];
                          matched_group_sigma_v+=pow(v_i-v_mean,2.);
                       }
                       current_j=current_j->substructure_next;
                    }
                    matched_group_sigma_v=sqrt(matched_group_sigma_v/(double)(n_sub_lo-1));
                 }

              }
              else{
                 matched_subgroup           =matched;
                 matched_group              =matched_subgroup->parent;
                 matched_subgroup_properties=&(subgroup_properties[matched_subgroup->snap_tree][matched_subgroup->neighbour_index]);
                 matched_group_properties   =&(group_properties[matched_group->snap_tree][matched_group->neighbour_index]);
              }
              if(i_cat==0)
                 fprintf(fp_out,"%7d %7d %6d %5d %5d %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le\n",
                    current_group->neighbour_index,
                    matched_group->neighbour_index,
                    current_group_properties->n_particles,
                    n_sub_hi,
                    n_sub_lo,
                    current_group_sigma_v,
                    matched_group_sigma_v,
                    current_group_properties->sigma_v,
                    matched_group_properties->sigma_v,
                    current_group_properties->M_vir,
                    matched_group_properties->M_vir,
                    current_group_properties->position_MBP[0],
                    current_group_properties->position_MBP[1],
                    current_group_properties->position_MBP[2],
                    matched_group_properties->position_MBP[0],
                    matched_group_properties->position_MBP[1],
                    matched_group_properties->position_MBP[2],
                    current_group_properties->velocity_COM[0],
                    current_group_properties->velocity_COM[1],
                    current_group_properties->velocity_COM[2],
                    matched_group_properties->velocity_COM[0],
                    matched_group_properties->velocity_COM[1],
                    matched_group_properties->velocity_COM[2]);
              else
                 fprintf(fp_out,"%7d %7d %7d %7d %6d %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le %10.3le\n",
                    current_subgroup->neighbour_index,
                    matched_subgroup->neighbour_index,
                    current_group->neighbour_index,
                    matched_group->neighbour_index,
                    current_subgroup_properties->n_particles,
                    current_subgroup_properties->M_vir,
                    matched_subgroup_properties->M_vir,
                    current_group_properties->M_vir,
                    matched_group_properties->M_vir,
                    current_subgroup_properties->position_MBP[0],
                    current_subgroup_properties->position_MBP[1],
                    current_subgroup_properties->position_MBP[2],
                    matched_subgroup_properties->position_MBP[0],
                    matched_subgroup_properties->position_MBP[1],
                    matched_subgroup_properties->position_MBP[2],
                    current_subgroup_properties->velocity_COM[0],
                    current_subgroup_properties->velocity_COM[1],
                    current_subgroup_properties->velocity_COM[2],
                    matched_subgroup_properties->velocity_COM[0],
                    matched_subgroup_properties->velocity_COM[1],
                    matched_subgroup_properties->velocity_COM[2]);
           }
        }
        current=current->next_neighbour;
     }
     fclose(fp_out);
     SID_log("Done.",SID_LOG_CLOSE);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

