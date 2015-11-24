#define _MAIN
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

#define TAU_FORM_RECOVERED   2.0
#define TAU_MERGER_RECOVERED 3.0
#define XOFF_RELAXED         0.07
#define FSUB_RELAXED         0.1
#define VIR_RELAXED          1.35

// Structure that will carry the needed information to the select-and-analyze function
typedef struct process_trees_params_local process_trees_params_local;
struct process_trees_params_local{
   char  filename_output_root[MAX_FILENAME_LENGTH];
   FILE *fp_out;
};

void process_trees_fctn_init_snap_local(tree_info *trees,void *params_in,int mode,int i_type,int flag_init,int i_snap);
void process_trees_fctn_init_snap_local(tree_info *trees,void *params_in,int mode,int i_type,int flag_init,int i_snap){
   process_trees_params_local *params=(process_trees_params_local *)params_in;
   char filename_out[MAX_FILENAME_LENGTH];
   if(i_type==0)
      sprintf(filename_out,"%s_%03d_groups.txt",params->filename_output_root,trees->snap_list[i_snap]);
   else
      sprintf(filename_out,"%s_%03d_subgroups.txt",params->filename_output_root,trees->snap_list[i_snap]);
   if(SID.n_proc>1)
      sprintf(filename_out,"%s.%d",filename_out,SID.My_rank);
   params->fp_out=fopen(filename_out,"w");
   int i_column=1;
   if(i_type==0){
      fprintf(params->fp_out,"# ASCII group catalog of snap #%d of %s\n",trees->snap_list[i_snap],params->filename_output_root);
      fprintf(params->fp_out,"# File columns: (%02d)     group number\n",      i_column++);
   }
   else{
      fprintf(params->fp_out,"# ASCII subgroup catalog of snap #%d of %s\n",trees->snap_list[i_snap],params->filename_output_root);
      fprintf(params->fp_out,"# File columns: (%02d)     subgroup number\n",      i_column++);
   }
   fprintf(params->fp_out,"#               (%02d)     # of particles\n",      i_column++);
   fprintf(params->fp_out,"#               (%02d)     id_MBP\n",              i_column++);
   fprintf(params->fp_out,"#               (%02d)     x_COM      [Mpc/h]\n",  i_column++);
   fprintf(params->fp_out,"#               (%02d)     y_COM      [Mpc/h]\n",  i_column++);
   fprintf(params->fp_out,"#               (%02d)     z_COM      [Mpc/h]\n",  i_column++);
   fprintf(params->fp_out,"#               (%02d)     v_x_COM    [km/s]\n",   i_column++);
   fprintf(params->fp_out,"#               (%02d)     v_y_COM    [km/s]\n",   i_column++);
   fprintf(params->fp_out,"#               (%02d)     v_z_COM    [km/s]\n",   i_column++);
   fprintf(params->fp_out,"#               (%02d)     x_MBP      [Mpc/h]\n",  i_column++);
   fprintf(params->fp_out,"#               (%02d)     y_MBP      [Mpc/h]\n",  i_column++);
   fprintf(params->fp_out,"#               (%02d)     z_MBP      [Mpc/h]\n",  i_column++);
   fprintf(params->fp_out,"#               (%02d)     v_x_MBP    [km/s]\n",   i_column++);
   fprintf(params->fp_out,"#               (%02d)     v_y_MBP    [km/s]\n",   i_column++);
   fprintf(params->fp_out,"#               (%02d)     v_z_MBP    [km/s]\n",   i_column++);
   fprintf(params->fp_out,"#               (%02d)     M_vir      [M_sol/h]\n",i_column++);
   fprintf(params->fp_out,"#               (%02d)     R_vir      [Mpc/h]\n",  i_column++);
   fprintf(params->fp_out,"#               (%02d)     R_halo     [Mpc/h]\n",  i_column++);
   fprintf(params->fp_out,"#               (%02d)     R_max      [Mpc/h]\n",  i_column++);
   fprintf(params->fp_out,"#               (%02d)     V_max      [km/s]\n",   i_column++);
   fprintf(params->fp_out,"#               (%02d)     V_vir      [km/s]\n",   i_column++);
   fprintf(params->fp_out,"#               (%02d)     sigma_v    [km/s]\n",   i_column++);
   fprintf(params->fp_out,"#               (%02d,%02d,%02d) spin       [Mpc/h km/s]\n",i_column,i_column+1,i_column+2);i_column+=3;
   fprintf(params->fp_out,"#               (%02d)     lambda (spin parameter)\n",i_column++);
   fprintf(params->fp_out,"#               (%02d)     q_triaxial\n",             i_column++);
   fprintf(params->fp_out,"#               (%02d)     s_triaxial\n",             i_column++);
   fprintf(params->fp_out,"#               (%02d)     offset_COM [R_vir]\n",     i_column++);
   if(i_type==0){
      fprintf(params->fp_out,"#               (%02d)     n_substructures\n",     i_column++);
   }
   else{
      fprintf(params->fp_out,"#               (%02d)     group index\n",           i_column++);
   }
   fprintf(params->fp_out,"#               (%02d)     f_sub\n",i_column++);
   fprintf(params->fp_out,"#               (%02d)     M_peak [M_sol/h]\n",         i_column++);
   fprintf(params->fp_out,"#               (%02d)     tau_form\n",           i_column++);
   fprintf(params->fp_out,"#               (%02d)     tau_3to1\n",           i_column++);
   fprintf(params->fp_out,"#               (%02d)     tau_10to1\n",          i_column++);
}

void process_trees_fctn_fin_snap_local(tree_info *trees,void *params_in,int mode,int i_type,int flag_init,int i_snap);
void process_trees_fctn_fin_snap_local(tree_info *trees,void *params_in,int mode,int i_type,int flag_init,int i_snap){
   process_trees_params_local *params=(process_trees_params_local *)params_in;
   fclose(params->fp_out);
}

// ** Define the calculation here **
void process_trees_fctn_init_local(tree_info *trees,void *params_in,int mode,int i_type);
void process_trees_fctn_init_local(tree_info *trees,void *params_in,int mode,int i_type){
  // Create an alias for the passed void pointer
  process_trees_params_local *params=(process_trees_params_local *)params_in;
  // Initialize markers (just one-at-a-time to save RAM)
  if(i_type==0)
     precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
  else
     precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
}

// ** Perform halo-by-halo processing here **
void process_trees_fctn_analyze_local(tree_info *trees,void *params_in,int mode,int i_type,int flag_init,tree_node_info *halo);
void process_trees_fctn_analyze_local(tree_info *trees,void *params_in,int mode,int i_type,int flag_init,tree_node_info *halo){
   process_trees_params_local *params=(process_trees_params_local *)params_in;

   int n_subgroups       =halo->n_substructures;
   int flag_process_group=(i_type==0);
   halo_properties_info *properties      =fetch_treenode_properties(trees,halo);
   halo_properties_info *properties_group=fetch_treenode_properties(trees,halo->parent);
   double expansion_factor=a_of_z(trees->z_list[halo->snap_tree]);
   double box_size=(double)trees->box_size;
   int    i_group=halo->file_index;

   // Create a few properties
   double dx,dy,dz,v_c,lambda;
   float  offset_COM;
   float  r_min,r_max;
   v_c       =sqrt(G_NEWTON*properties->M_vir*M_SOL/(properties->R_vir*M_PER_MPC))*1e-3;
   lambda    =sqrt(properties->spin[0]*properties->spin[0]+properties->spin[1]*properties->spin[1]+properties->spin[2]*properties->spin[2])/(sqrt(2.)*properties->R_vir*v_c);
   dx        =d_periodic(properties->position_COM[0]-properties->position_MBP[0],box_size);
   dy        =d_periodic(properties->position_COM[1]-properties->position_MBP[1],box_size);
   dz        =d_periodic(properties->position_COM[2]-properties->position_MBP[2],box_size);
   offset_COM=sqrt(dx*dx+dy*dy+dz*dz)/(properties->R_vir/expansion_factor);

   // Perform write
   fprintf(params->fp_out,"%9d %9d %9lld  %12.5e %12.5e %12.5e %11.5f %11.5f %11.5f  %12.5e %12.5e %12.5e %11.5f %11.5f %11.5f  %10.5le %11.5f %11.5f %11.5f  %11.5f %11.5f %11.5f  %12.5e %12.5e %12.5e %11.5f  %11.5f %11.5f  %11.5f",
           i_group,properties->n_particles,properties->id_MBP,
           properties->position_COM[0],properties->position_COM[1],properties->position_COM[2],
           properties->velocity_COM[0],properties->velocity_COM[1],properties->velocity_COM[2],
           properties->position_MBP[0],properties->position_MBP[1],properties->position_MBP[2],
           properties->velocity_MBP[0],properties->velocity_MBP[1],properties->velocity_MBP[2],
           properties->M_vir,
           properties->R_vir,
           properties->R_halo,
           properties->R_max,
           properties->V_max,
           v_c,
           properties->sigma_v,
           properties->spin[0],properties->spin[1],properties->spin[2],
           lambda,
           properties->q_triaxial,
           properties->s_triaxial,
           offset_COM); // converts to kpc/h
    if(!flag_process_group){
       float f_sub=(float)properties->n_particles/(float)properties_group->n_particles;
       fprintf(params->fp_out,"  %9d %11.5f",halo->parent->file_index,f_sub);
    }
    else{
       double f_sub=fetch_treenode_SSFctn(trees,halo);
       fprintf(params->fp_out,"  %3d %11.5lf",n_subgroups,f_sub);
    }
    double tau_form =fetch_treenode_tau_form (trees,halo);
    double tau_3to1 =fetch_treenode_tau_3to1 (trees,halo);
    double tau_10to1=fetch_treenode_tau_10to1(trees,halo);
    double M_peak   =fetch_treenode_M_peak(trees,halo);
    fprintf(params->fp_out,"  %10.5le %11.5lf  %11.5lf  %11.5lf",M_peak,tau_form,tau_3to1,tau_10to1);
    fprintf(params->fp_out,"\n");
}

// ** Write the results here **
void process_trees_fctn_fin_local(tree_info *trees,void *params_in,int mode,int i_type);
void process_trees_fctn_fin_local(tree_info *trees,void *params_in,int mode,int i_type){
   // Create an alias for the passed void pointer
   process_trees_params_local *params=(process_trees_params_local *)params_in;
   // Clean-up markers
   if(i_type==0)
      free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
   else
      free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
}

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  char filename_SSimPL_dir[MAX_FILENAME_LENGTH];
  char filename_halo_version_root[MAX_FILENAME_LENGTH];
  char filename_trees_name[MAX_FILENAME_LENGTH];
  char filename_output_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_dir,       argv[1]);
  strcpy(filename_halo_version_root,argv[2]);
  strcpy(filename_trees_name,       argv[3]);
  int  snap_number_start      =atoi(argv[4]);
  int  snap_number_stop       =atoi(argv[5]);
  strcpy(filename_output_root,      argv[6]);

  // Set the halo and tree filename roots
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_dir,filename_trees_name);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_dir,filename_halo_version_root);

  SID_log("Generating ascii version of catalogs with tree markers...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Perform analysis
  tree_info *trees;
  read_trees(filename_SSimPL_dir,
             filename_halo_version_root,
             filename_trees_name,
             TREE_MODE_DEFAULT,
             &trees);

  // Convert given snapshot ranges to snap_tree ranges
  snap_number_start=find_treesnap_snap(trees,snap_number_start);
  snap_number_stop =find_treesnap_snap(trees,snap_number_stop);

  // Read catalogs
  read_trees_catalogs(trees,
                      filename_SSimPL_dir,
                      filename_halo_version_root,
                      READ_TREES_CATALOGS_BOTH);

  // ** PERFORM the calculation here **a
  process_trees_params_local params;
  strcpy(params.filename_output_root,filename_output_root);
  process_trees_by_snap(trees,&params,PROCESS_TREES_BOTH,snap_number_start,snap_number_stop-snap_number_start+1,
                        process_trees_fctn_init_local,
                        process_trees_fctn_init_snap_local,
                        process_trees_fctn_select_null,
                        process_trees_fctn_analyze_local,
                        process_trees_fctn_fin_snap_local,
                        process_trees_fctn_fin_local);

  // Clean-up
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

