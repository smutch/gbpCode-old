#define _MAIN
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

#define TAU_FORM_RECOVERED   1.5
#define TAU_MERGER_RECOVERED 2.0
#define XOFF_RELAXED         0.07
#define FSUB_RELAXED         0.1
#define VIR_RELAXED          1.35

// Structure that will carry the needed information to the select-and-analyze function
typedef struct process_trees_params_local process_trees_params_local;
struct process_trees_params_local{
   int          flag_calc_Vir;
   char         filename_output_root[MAX_FILENAME_LENGTH];
   int          snap_lo_tau_trends;
   int          snap_hi_tau_trends;
   trend_info  *mass_binning;
   trend_info **trends_z;
   trend_info **trends_tau_form;
   trend_info **trends_tau_3to1;
   trend_info **trends_tau_10to1;
   int        **n_halos_z;
   int        **n_halos_z_recovered_form;
   int        **n_halos_z_recovered_3to1;
   int        **n_halos_z_recovered_10to1;
   int        **n_halos_z_relaxed_x_off;
   int        **n_halos_z_relaxed_f_sub;
   int        **n_halos_z_relaxed_Vir_ratio;
   int        **n_halos_z_relaxed_all;
};

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
  // Create the trend which will describe the mass binning
  init_treenode_trend(trees,&(params->mass_binning),"logM_course");
  // Allocate the trends
  int n_M=params->mass_binning->ordinate->hist->n_bins;
  params->trends_z        =(trend_info **)SID_malloc(sizeof(trend_info *)*n_M);
  params->trends_tau_form =(trend_info **)SID_malloc(sizeof(trend_info *)*n_M);
  params->trends_tau_3to1 =(trend_info **)SID_malloc(sizeof(trend_info *)*n_M);
  params->trends_tau_10to1=(trend_info **)SID_malloc(sizeof(trend_info *)*n_M);
  // Initialize the trend(s) to be populated
  for(int i_M=0;i_M<n_M;i_M++){
     // f(z) trends
     init_treenode_trend           (trees,&(params->trends_z[i_M]),"z");
     init_treenode_trend_coordinate(trees, (params->trends_z[i_M]),"tau_form");
     init_treenode_trend_coordinate(trees, (params->trends_z[i_M]),"tau_3to1");
     init_treenode_trend_coordinate(trees, (params->trends_z[i_M]),"tau_10to1");
     init_treenode_trend_coordinate(trees, (params->trends_z[i_M]),"xoff");
     if(params->flag_calc_Vir)
        init_treenode_trend_coordinate(trees, (params->trends_z[i_M]),"Vir_ratio");
     if(i_type==0)
        init_treenode_trend_coordinate(trees,(params->trends_z[i_M]),"SSFctn");
     // f(tau_form) trends
     init_treenode_trend           (trees,&(params->trends_tau_form[i_M]),"tau_form");
     init_treenode_trend_coordinate(trees, (params->trends_tau_form[i_M]),"xoff");
     if(params->flag_calc_Vir)
        init_treenode_trend_coordinate(trees, (params->trends_tau_form[i_M]),"Vir_ratio");
     if(i_type==0)
        init_treenode_trend_coordinate(trees, (params->trends_tau_form[i_M]),"SSFctn");
     // f(tau_3to1) trends
     init_treenode_trend           (trees,&(params->trends_tau_3to1[i_M]),"tau_3to1");
     init_treenode_trend_coordinate(trees, (params->trends_tau_3to1[i_M]),"xoff");
     if(params->flag_calc_Vir)
        init_treenode_trend_coordinate(trees, (params->trends_tau_3to1[i_M]),"Vir_ratio");
     if(i_type==0)
        init_treenode_trend_coordinate(trees, (params->trends_tau_3to1[i_M]), "SSFctn");
     // f(tau_10to1) trends
     init_treenode_trend           (trees,&(params->trends_tau_10to1[i_M]),"tau_10to1");
     init_treenode_trend_coordinate(trees, (params->trends_tau_10to1[i_M]),"xoff");
     if(params->flag_calc_Vir)
        init_treenode_trend_coordinate(trees, (params->trends_tau_10to1[i_M]),"Vir_ratio");
     if(i_type==0)
        init_treenode_trend_coordinate(trees, (params->trends_tau_10to1[i_M]),"SSFctn");
  }
  // Initialize recovery/relaxation arrays
  params->n_halos_z                  =(int **)SID_calloc(sizeof(int *)*n_M);
  params->n_halos_z_recovered_form   =(int **)SID_calloc(sizeof(int *)*n_M);
  params->n_halos_z_recovered_3to1   =(int **)SID_calloc(sizeof(int *)*n_M);
  params->n_halos_z_recovered_10to1  =(int **)SID_calloc(sizeof(int *)*n_M);
  params->n_halos_z_relaxed_x_off    =(int **)SID_calloc(sizeof(int *)*n_M);
  params->n_halos_z_relaxed_f_sub    =(int **)SID_calloc(sizeof(int *)*n_M);
  params->n_halos_z_relaxed_Vir_ratio=(int **)SID_calloc(sizeof(int *)*n_M);
  params->n_halos_z_relaxed_all      =(int **)SID_calloc(sizeof(int *)*n_M);
  for(int i_M=0;i_M<n_M;i_M++){
     int n_z=params->trends_z[i_M]->ordinate->hist->n_bins;
     params->n_halos_z[i_M]                  =(int *)SID_calloc(sizeof(int)*n_z);
     params->n_halos_z_recovered_form[i_M]   =(int *)SID_calloc(sizeof(int)*n_z);
     params->n_halos_z_recovered_3to1[i_M]   =(int *)SID_calloc(sizeof(int)*n_z);
     params->n_halos_z_recovered_10to1[i_M]  =(int *)SID_calloc(sizeof(int)*n_z);
     params->n_halos_z_relaxed_x_off[i_M]    =(int *)SID_calloc(sizeof(int)*n_z);
     params->n_halos_z_relaxed_f_sub[i_M]    =(int *)SID_calloc(sizeof(int)*n_z);
     params->n_halos_z_relaxed_Vir_ratio[i_M]=(int *)SID_calloc(sizeof(int)*n_z);
     params->n_halos_z_relaxed_all[i_M]      =(int *)SID_calloc(sizeof(int)*n_z);
  }
}

// ** Perform the calculation here **
void process_trees_fctn_analyze_local(tree_info *trees,void *params_in,int mode,int i_type,int flag_init,tree_node_info *halo);
void process_trees_fctn_analyze_local(tree_info *trees,void *params_in,int mode,int i_type,int flag_init,tree_node_info *halo){
   process_trees_params_local *params=(process_trees_params_local *)params_in;
   int i_M;
   int i_snap   =halo->snap_tree;
   int i_snap_lo=((process_trees_params_local *)params)->snap_lo_tau_trends;
   int i_snap_hi=((process_trees_params_local *)params)->snap_hi_tau_trends;
   if((i_M=add_item_to_trend(((process_trees_params_local *)params)->mass_binning,GBP_ADD_ITEM_TO_TREND_DEFAULT,halo))>=0){
      // Add halo to redshift trend
      int i_z=add_item_to_trend(((process_trees_params_local *)params)->trends_z[i_M],GBP_ADD_ITEM_TO_TREND_DEFAULT,halo);
      // Add halo to redshift arrays
      if(i_z>=0 && i_z<((process_trees_params_local *)params)->trends_z[i_M]->ordinate->hist->n_bins){
         double tau_form_i =fetch_treenode_tau_form (trees,halo);
         double tau_3to1_i =fetch_treenode_tau_3to1 (trees,halo);
         double tau_10to1_i=fetch_treenode_tau_10to1(trees,halo);
         double x_off_i    =fetch_treenode_x_off    (trees,halo);
         double f_sub_i    =fetch_treenode_SSFctn   (trees,halo);
         double vir_i;
         if(params->flag_calc_Vir)
            vir_i=fetch_treenode_Vir_ratio(trees,halo);
         else
            vir_i=-1;
         params->n_halos_z[i_M][i_z]++;
         if(tau_form_i >=TAU_FORM_RECOVERED)   params->n_halos_z_recovered_form[i_M][i_z]++;
         if(tau_3to1_i >=TAU_MERGER_RECOVERED) params->n_halos_z_recovered_3to1[i_M][i_z]++;
         if(tau_10to1_i>=TAU_MERGER_RECOVERED) params->n_halos_z_recovered_10to1[i_M][i_z]++;
         if(x_off_i    <=XOFF_RELAXED)         params->n_halos_z_relaxed_x_off[i_M][i_z]++;
         if(f_sub_i    <=FSUB_RELAXED)         params->n_halos_z_relaxed_f_sub[i_M][i_z]++;
         if(vir_i      <=VIR_RELAXED)          params->n_halos_z_relaxed_Vir_ratio[i_M][i_z]++;
         if((x_off_i   <=XOFF_RELAXED) &&
            (f_sub_i   <=FSUB_RELAXED) &&       
            (vir_i     <=VIR_RELAXED))         params->n_halos_z_relaxed_all[i_M][i_z]++;
      }
      // Add halo to tau trends
      if(i_snap>=i_snap_lo&&i_snap<=i_snap_hi){
         add_item_to_trend(((process_trees_params_local *)params)->trends_tau_form[i_M], GBP_ADD_ITEM_TO_TREND_DEFAULT,halo);
         add_item_to_trend(((process_trees_params_local *)params)->trends_tau_3to1[i_M], GBP_ADD_ITEM_TO_TREND_DEFAULT,halo);
         add_item_to_trend(((process_trees_params_local *)params)->trends_tau_10to1[i_M],GBP_ADD_ITEM_TO_TREND_DEFAULT,halo);
      }
   }
}

// ** Write the results here **
void process_trees_fctn_fin_local(tree_info *trees,void *params_in,int mode,int i_type);
void process_trees_fctn_fin_local(tree_info *trees,void *params_in,int mode,int i_type){
   // Create an alias for the passed void pointer
   process_trees_params_local *params=(process_trees_params_local *)params_in;
   // Clean-up markers
   char filename_output_root_root[MAX_FILENAME_LENGTH];
   if(i_type==0){
      sprintf(filename_output_root_root,"%s_groups",params->filename_output_root);
      free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
   }
   else{
      sprintf(filename_output_root_root,"%s_subgroups",params->filename_output_root);
      free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
   }

   // Finalize recovery/relaxation arrays
   int n_M=params->mass_binning->ordinate->hist->n_bins;
   for(int i_M=0;i_M<n_M;i_M++){
      int n_z=params->trends_z[i_M]->ordinate->hist->n_bins;
      SID_Allreduce(SID_IN_PLACE,params->n_halos_z[i_M],                  n_z,SID_INT,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,params->n_halos_z_recovered_form[i_M],   n_z,SID_INT,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,params->n_halos_z_recovered_3to1[i_M],   n_z,SID_INT,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,params->n_halos_z_recovered_10to1[i_M],  n_z,SID_INT,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,params->n_halos_z_relaxed_x_off[i_M],    n_z,SID_INT,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,params->n_halos_z_relaxed_f_sub[i_M],    n_z,SID_INT,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,params->n_halos_z_relaxed_Vir_ratio[i_M],n_z,SID_INT,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,params->n_halos_z_relaxed_all[i_M],      n_z,SID_INT,SID_SUM,SID.COMM_WORLD);
   }

   // Write recovery/relaxation arrays
   for(int i_M=0;i_M<n_M;i_M++){
      char  filename_out[MAX_FILENAME_LENGTH];
      char mass_text_lo[64];
      char mass_text_hi[64];
      float_to_text((float)histogram_bin_x_lo(params->mass_binning->ordinate->hist,i_M),1,mass_text_lo);
      float_to_text((float)histogram_bin_x_hi(params->mass_binning->ordinate->hist,i_M),1,mass_text_hi);
      sprintf(filename_out,"%s_logM_%s_to_%s_relaxed.txt",filename_output_root_root,mass_text_lo,mass_text_hi);
      FILE *fp_out=fopen(filename_out,"w");
      int   n_z   =params->trends_z[i_M]->ordinate->hist->n_bins;
      int   i_column=1;
      fprintf(fp_out,"# Column (%02d): Snapshot\n",                                    i_column++);
      fprintf(fp_out,"#        (%02d): Redshift\n",                                    i_column++);
      fprintf(fp_out,"#        (%02d): No. of halos at z\n",                           i_column++);
      fprintf(fp_out,"#        (%02d): No. of halos w/ tau_form     >=%.1lf\n",        i_column++,TAU_FORM_RECOVERED);
      fprintf(fp_out,"#        (%02d): No. of halos w/ tau_3to1     >=%.1lf\n",        i_column++,TAU_MERGER_RECOVERED);
      fprintf(fp_out,"#        (%02d): No. of halos w/ tau_10to1    >=%.1lf\n",        i_column++,TAU_MERGER_RECOVERED);
      fprintf(fp_out,"#        (%02d): No. of halos w/ x_off        <=%.2lf\n",        i_column++,XOFF_RELAXED);
      fprintf(fp_out,"#        (%02d): No. of halos w/ f_sub        <=%.1lf\n",        i_column++,FSUB_RELAXED);
      if(params->flag_calc_Vir)
         fprintf(fp_out,"#        (%02d): No. of halos w/ Virial Ratio <=%.1lf\n",        i_column++,VIR_RELAXED);
      fprintf(fp_out,"#        (%02d): No. of halos meeting all relaxation criteria\n",i_column++);
      if(params->flag_calc_Vir)
         for(int i_z=0;i_z<n_z;i_z++)
            fprintf(fp_out,"%3d %6.2lf %6d %6d %6d %6d %6d %6d %6d %6d\n",
                           trees->snap_list[i_z],
                           trees->z_list[i_z],
                           params->n_halos_z[i_M][i_z],
                           params->n_halos_z_recovered_form[i_M][i_z],
                           params->n_halos_z_recovered_3to1[i_M][i_z],
                           params->n_halos_z_recovered_10to1[i_M][i_z],
                           params->n_halos_z_relaxed_x_off[i_M][i_z],
                           params->n_halos_z_relaxed_f_sub[i_M][i_z],
                           params->n_halos_z_relaxed_Vir_ratio[i_M][i_z],
                           params->n_halos_z_relaxed_all[i_M][i_z]);
      else
         for(int i_z=0;i_z<n_z;i_z++)
            fprintf(fp_out,"%3d %6.2lf %6d %6d %6d %6d %6d %6d %6d\n",
                           trees->snap_list[i_z],
                           trees->z_list[i_z],
                           params->n_halos_z[i_M][i_z],
                           params->n_halos_z_recovered_form[i_M][i_z],
                           params->n_halos_z_recovered_3to1[i_M][i_z],
                           params->n_halos_z_recovered_10to1[i_M][i_z],
                           params->n_halos_z_relaxed_x_off[i_M][i_z],
                           params->n_halos_z_relaxed_f_sub[i_M][i_z],
                           params->n_halos_z_relaxed_all[i_M][i_z]);
      fclose(fp_out);
   }

   // Free recovery/relaxation arrays
   for(int i_M=0;i_M<n_M;i_M++){
      SID_free(SID_FARG params->n_halos_z[i_M]);
      SID_free(SID_FARG params->n_halos_z_recovered_form[i_M]);
      SID_free(SID_FARG params->n_halos_z_recovered_3to1[i_M]);
      SID_free(SID_FARG params->n_halos_z_recovered_10to1[i_M]);
      SID_free(SID_FARG params->n_halos_z_relaxed_x_off[i_M]);
      SID_free(SID_FARG params->n_halos_z_relaxed_f_sub[i_M]);
      SID_free(SID_FARG params->n_halos_z_relaxed_Vir_ratio[i_M]);
   }
   SID_free(SID_FARG params->n_halos_z);
   SID_free(SID_FARG params->n_halos_z_recovered_form);
   SID_free(SID_FARG params->n_halos_z_recovered_3to1);
   SID_free(SID_FARG params->n_halos_z_recovered_10to1);
   SID_free(SID_FARG params->n_halos_z_relaxed_x_off);
   SID_free(SID_FARG params->n_halos_z_relaxed_f_sub);
   SID_free(SID_FARG params->n_halos_z_relaxed_Vir_ratio);

   // Finalize & write trend results
   finalize_trend(params->mass_binning); // needed for the check on bin_count (below) to work w/ n_proc>1
   for(int i_M=0;i_M<n_M;i_M++){
      if(params->mass_binning->ordinate->hist->bin_count[i_M]>50){
         // Set the root of the output filename
         char filename_output_root[MAX_FILENAME_LENGTH];
         char mass_text_lo[64];
         char mass_text_hi[64];
         float_to_text((float)histogram_bin_x_lo(params->mass_binning->ordinate->hist,i_M),1,mass_text_lo);
         float_to_text((float)histogram_bin_x_hi(params->mass_binning->ordinate->hist,i_M),1,mass_text_hi);
         sprintf(filename_output_root,"%s_logM_%s_to_%s",filename_output_root_root,mass_text_lo,mass_text_hi);
         // Write the trends
         write_trend_ascii(params->trends_z[i_M],        filename_output_root);
         write_trend_ascii(params->trends_tau_form[i_M], filename_output_root);
         write_trend_ascii(params->trends_tau_3to1[i_M], filename_output_root);
         write_trend_ascii(params->trends_tau_10to1[i_M],filename_output_root);
      }
      // Clean-up the trends
      free_trend(&(params->trends_z[i_M]));
      free_trend(&(params->trends_tau_form[i_M]));
      free_trend(&(params->trends_tau_3to1[i_M]));
      free_trend(&(params->trends_tau_10to1[i_M]));
   }
   SID_free(SID_FARG params->trends_z);
   SID_free(SID_FARG params->trends_tau_form);
   SID_free(SID_FARG params->trends_tau_3to1);
   SID_free(SID_FARG params->trends_tau_10to1);

   // Write and clean-up the mass-binning trend
   write_trend_ascii(params->mass_binning,filename_output_root_root);
   free_trend(&(params->mass_binning));
}

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  char filename_SSimPL_root[MAX_FILENAME_LENGTH];
  char filename_halo_version_root[MAX_FILENAME_LENGTH];
  char filename_trees_name[MAX_FILENAME_LENGTH];
  char filename_data_root[MAX_FILENAME_LENGTH];
  char filename_output_root[MAX_FILENAME_LENGTH];
  int i_arg       =1;
  int flag_calc_Vir=TRUE;
  strcpy(filename_SSimPL_root,       argv[i_arg++]);
  strcpy(filename_halo_version_root,argv[i_arg++]);
  strcpy(filename_trees_name,       argv[i_arg++]);
  if(argc==8)
     strcpy(filename_data_root,      argv[i_arg++]);
  else
     flag_calc_Vir=FALSE;
  strcpy(filename_output_root,      argv[i_arg++]);
  double z_lo_tau_trends      =atof(argv[i_arg++]);
  double z_hi_tau_trends      =atof(argv[i_arg++]);

  // Set the halo and tree filename roots
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_root,filename_trees_name);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_root,filename_halo_version_root);

  SID_log("Generating treenode markers & analysis of merger trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Perform analysis
  tree_info *trees;
  read_trees(filename_SSimPL_root,
             filename_halo_version_root,
             filename_trees_name,
             TREE_MODE_DEFAULT,
             &trees);

  // Read catalogs
  read_trees_catalogs(trees,READ_TREES_CATALOGS_BOTH);

  // Read NFW concentrations (if given)
  if(flag_calc_Vir){
     // Set the filename root
     char filename_sim_root[MAX_FILENAME_LENGTH];
     strcpy(filename_sim_root,filename_SSimPL_root);
     strip_path(filename_sim_root);
     char filename_NFW_root[MAX_FILENAME_LENGTH];
     sprintf(filename_NFW_root,"%s/NFW_concentrations_%s",filename_data_root,filename_sim_root);

     read_trees_data(trees,
                     filename_NFW_root,
                     READ_TREES_DATA_SUBGROUPS,
                     SID_FLOAT,
                     "c_NFW");

     // Set group concentrations to central subgroup concentrations
     float **group_concentrations_local;
     float **subgroup_concentrations_local;
     init_trees_data(trees,(void ***)&group_concentrations_local,sizeof(float),INIT_TREE_DATA_GROUPS,"c_NFW_groups");
     subgroup_concentrations_local=(float **)ADaPS_fetch(trees->data,"c_NFW_subgroups");
     for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
        tree_node_info *group_current=trees->first_neighbour_groups[i_snap];
        while(group_current!=NULL){
           tree_node_info *subgroup_current=group_current->substructure_first;
           tree_node_info *subgroup_central=subgroup_current;
           while(subgroup_current!=NULL){
              if(subgroup_current->n_particles>subgroup_central->n_particles)
                 subgroup_central=subgroup_current;
              subgroup_current=subgroup_current->substructure_next;
           }
           if(subgroup_central!=NULL){
              group_concentrations_local[group_current->snap_tree][group_current->neighbour_index]=
                 subgroup_concentrations_local[subgroup_central->snap_tree][subgroup_central->neighbour_index];
           }
           else
              group_concentrations_local[group_current->snap_tree][group_current->neighbour_index]=0.;
           group_current=group_current->next_neighbour;
        }
     }
  }

  // Populate the structure that gets passed to the calculation
  process_trees_params_local params;
  params.flag_calc_Vir=flag_calc_Vir;
  strcpy(params.filename_output_root,filename_output_root);
  params.snap_hi_tau_trends=find_treesnap_z(trees,z_lo_tau_trends); // snap_tree and z_list run in opposite orders
  params.snap_lo_tau_trends=find_treesnap_z(trees,z_hi_tau_trends); // snap_tree and z_list run in opposite orders
  SID_log("Trends with tau will be compiled between (z,snap)=(%lf,%d) and (z,snap)=(%lf,%d).",SID_LOG_COMMENT,
          trees->z_list[params.snap_hi_tau_trends],params.snap_hi_tau_trends,
          trees->z_list[params.snap_lo_tau_trends],params.snap_lo_tau_trends);

  // ** PERFORM the calculation here **
  process_trees_by_snap(trees,&params,PROCESS_TREES_BOTH,0,trees->n_snaps,
                        process_trees_fctn_init_local,
                        process_trees_fctn_init_snap_null,
                        process_trees_fctn_select_null,
                        process_trees_fctn_analyze_local,
                        process_trees_fctn_fin_snap_null,
                        process_trees_fctn_fin_local);

  // Clean-up
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

