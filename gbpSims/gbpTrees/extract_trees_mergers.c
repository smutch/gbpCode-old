#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

typedef struct walk_trees_local_info walk_trees_local_info;
struct walk_trees_local_info{
   float   M_sat;
   int     n_snaps;
   double *a_list;
   int     i_snap_host;
   int     n_halos;
   int     n_halos_total;
   int     n_mergers;
   int     n_mergers_total;
   FILE   *fp_out_mergers;
   FILE   *fp_out_halos;
   double  z_min;
   double  z_max;
   int     n_z;
   double  z_step;
   double  M_min;
   double  M_max;
   int     n_M;
   double  M_step;
   double  theta_min;
   double  theta_max;
   int     n_theta;
   double  theta_step;
   int    *halo_array;
   int    *mergers_array;
};

int walk_trees_recursive_local(halo_properties_SAGE_info *tree,int i_node,walk_trees_local_info *data){
  int        i_current;
  int        i_progenitor;
  int        n_halos_written=0;
  float      M_halo;
  float      theta;

  M_halo=(tree[i_node].M_vir*1e10);

  // This information gets changed every step DOWN the tree
  int i_M;
  data->M_sat      =0.;
  data->i_snap_host=tree[i_node].snap_num;
  i_M   =(int)((take_log10(M_halo)-data->M_min)/data->M_step);
  if(i_M>=0 && i_M<data->n_M)
     data->halo_array[(data->i_snap_host)*data->n_M+i_M]++;

  // Write halo information
  fprintf(data->fp_out_halos,"%le %le\n",data->a_list[data->i_snap_host],M_halo);
  data->n_halos++;
  data->n_halos_total++;

  // Perform the walk
  n_halos_written++;
  i_current   =tree[i_node].progenitor_first;
  i_progenitor=0;
  while(i_current>=0){
    n_halos_written+=walk_trees_recursive_local(tree,i_current,data);

    // Perform actions on mergers
    if(i_progenitor>0){
       int i_theta;
       theta=data->M_sat/M_halo;
       fprintf(data->fp_out_mergers,"%le %le %le\n",data->a_list[data->i_snap_host],M_halo,theta);
       i_M    =(int)((take_log10(data->M_sat)-data->M_min)/data->M_step);
       i_theta=(int)((take_log10(theta)-data->theta_min)/data->theta_step);
       if(i_M>=0 && i_M<data->n_M && i_theta>=0 && i_theta<data->n_theta){
          int index;
          index=(data->i_snap_host*data->n_M+i_M)*data->n_theta+i_theta;
          data->mergers_array[index]++;
//printf("%3d(%3d) %3d(%3d) %3d(%3d) [%7d,%7d] -- %le %le\n",data->i_snap_host,data->n_snaps,i_M,data->n_M,i_theta,data->n_theta,index,data->mergers_array[index],M_halo,data->M_sat);
       }
       data->n_mergers++;
       data->n_mergers_total++;
    }

    i_progenitor++;    
    i_current=tree[i_current].progenitor_next;
  }

  // This information gets changed every step UP the tree
  data->M_sat=MAX(data->M_sat,M_halo);

  return(n_halos_written);
}

int main(int argc, char *argv[]){
  char        filename_root_in[256];
  char        filename_root_in_strip[256];
  char        filename_root_out[256];
  char        filename_out_mergers[256];
  char        filename_out_halos[256];
  char        filename_trees[256];
  char        filename_a_list[256];
  FILE       *fp_in;
  FILE       *fp_out;
  char       *line=NULL;
  size_t      line_length=0;
  int         i_read;
  cosmo_info *cosmo;
  walk_trees_local_info data;
 
  SID_init(&argc,&argv,NULL);

  // Initialize cosmology
  init_cosmo_std(&cosmo);

  // Fetch user inputs
  if(argc!=3)
    SID_trap_error("Incorrect syntax",ERROR_SYNTAX);
  strcpy(filename_root_in, argv[1]);
  strcpy(filename_root_out,argv[2]);

  strcpy(filename_root_in_strip,filename_root_in);
  strip_path(filename_root_in_strip);
  SID_log("Extracting mergers from trees {%s} and writing results to {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_root_in_strip,filename_root_out);

  // Read snapshot expansion factor list
  if(SID.I_am_Master){
    sprintf(filename_a_list,"%s/%s.a_list",filename_root_in,filename_root_in_strip);
    SID_log("Reading snapshot list {%s}...",SID_LOG_OPEN,filename_a_list);

    // Read the original list
    fp_in       =fopen(filename_a_list, "r");
    data.n_snaps=count_lines_data(fp_in);
    data.a_list =(double *)SID_malloc(sizeof(double)*data.n_snaps);
    for(i_read=0;i_read<data.n_snaps;i_read++){
      grab_next_line_data(fp_in,&line,&line_length);
      grab_double(line,1,&(data.a_list[i_read]));
    }
    fclose(fp_in);

    SID_log("Done.",SID_LOG_CLOSE);
  }

  // Set binning parameters
  data.z_min        = 0.;
  data.z_max        = 5.;
  data.n_z          =data.n_snaps;
  data.z_step       =(data.z_max-data.z_min)/(float)data.n_z;
  data.M_min        =take_log10(1e10);
  data.M_max        =take_log10(1e15);
  data.n_M          =50;
  data.M_step       =(data.M_max-data.M_min)/(float)data.n_M;
  data.theta_min    =take_log10(0.001);
  data.theta_max    =take_log10(1.000);
  data.n_theta      =50;
  data.theta_step   =(data.theta_max-data.theta_min)/(float)data.n_theta;
  data.halo_array   =(int *)SID_calloc(sizeof(int)*data.n_M*data.n_z);
  data.mergers_array=(int *)SID_calloc(sizeof(int)*data.n_M*data.n_z*data.n_theta);

  // Determine how many tree files we need to process
  int   n_files;
  fp_in  =NULL;
  n_files=0;
  do{
    if(fp_in!=NULL)
       fclose(fp_in);
    sprintf(filename_trees,"%s/vertical/subgroups/%s.subgroup_trees.%d",filename_root_in,filename_root_in_strip,n_files);
    fp_in=fopen(filename_trees,"r");
    if(fp_in!=NULL)
       n_files++;
  } while(fp_in!=NULL);

  // Open output files
  sprintf(filename_out_halos,  "%s.halos",  filename_root_out);
  sprintf(filename_out_mergers,"%s.mergers",filename_root_out);
  data.fp_out_halos  =fopen(filename_out_halos,  "w");
  data.fp_out_mergers=fopen(filename_out_mergers,"w");

  // Write file headers
  fprintf(data.fp_out_mergers,"# List of mergers for {%s}\n",filename_root_in_strip);
  fprintf(data.fp_out_mergers,"# Column (1): Expansion factor\n");
  fprintf(data.fp_out_mergers,"#        (2): Host mass\n");
  fprintf(data.fp_out_mergers,"#        (3): Merger ratio\n");
  fprintf(data.fp_out_halos,  "# List of halos for {%s}\n",filename_root_in_strip);
  fprintf(data.fp_out_halos,  "# Column (1): Expansion factor\n");
  fprintf(data.fp_out_halos,  "#        (2): Host mass\n");

  // Read each tree file in turn and write the desired output
  int i_file;
  data.n_halos_total  =0;
  data.n_mergers_total=0;
  for(i_file=0;i_file<n_files;i_file++){
  //for(i_file=0;i_file<1;i_file++){
     SID_log("Processing file %d of %d...",SID_LOG_OPEN,i_file+1,n_files);

     // Open file
     sprintf(filename_trees,"%s/vertical/subgroups/%s.subgroup_trees.%d",filename_root_in,filename_root_in_strip,i_file);
     fp_in=fopen(filename_trees,"r");

     // Read tree file header
     int        n_trees;
     int        n_halos;
     int       *n_halos_tree;
     int        n_halos_max;
     halo_properties_SAGE_info *tree;

     fread(&n_trees,sizeof(int),1,fp_in);
     fread(&n_halos,sizeof(int),1,fp_in);
     n_halos_tree=(int *)SID_malloc(sizeof(int)*n_trees);
     fread(n_halos_tree,sizeof(int),n_trees,fp_in);
     calc_max(n_halos_tree,&n_halos_max,n_trees,SID_INT,CALC_MODE_DEFAULT);
     tree=(halo_properties_SAGE_info *)SID_malloc(sizeof(halo_properties_SAGE_info)*n_halos_max);

     // Loop over each tree
     int i_tree;
     data.n_halos  =0;
     data.n_mergers=0;
     for(i_tree=0;i_tree<n_trees;i_tree++){
        // Read tree
        fread(tree,sizeof(halo_properties_SAGE_info),n_halos_tree[i_tree],fp_in);

        // Walk the tree and generate the desired output
        int i_halo=0;
        while(i_halo<n_halos_tree[i_tree])
           i_halo+=walk_trees_recursive_local(tree,i_halo,&data);

     }

     // Close file
     fclose(fp_in);
     SID_free(SID_FARG n_halos_tree);
     SID_log("# of halos  =%d",SID_LOG_COMMENT,data.n_halos);
     SID_log("# of mergers=%d",SID_LOG_COMMENT,data.n_mergers);
     SID_log("Done.",SID_LOG_CLOSE);
  }

  // Close the output file
  fclose(data.fp_out_mergers);
  fclose(data.fp_out_halos);

  SID_log("Total number of halos  =%d",SID_LOG_COMMENT,data.n_halos_total);
  SID_log("Total number of mergers=%d",SID_LOG_COMMENT,data.n_mergers_total);

  // Write binned arrays
  SID_log("Writing binned arrays...",SID_LOG_OPEN);
  sprintf(filename_out_halos,"%s.mergers_binned",filename_root_out);
  data.fp_out_halos=fopen(filename_out_halos,  "w");
  fwrite(&(data.n_M),       sizeof(int),1,                             data.fp_out_halos);
  fwrite(&(data.n_z),       sizeof(int),1,                             data.fp_out_halos);
  fwrite(&(data.n_theta),   sizeof(int),1,                             data.fp_out_halos);
  fwrite(data.halo_array,   sizeof(int),data.n_M*data.n_z,             data.fp_out_halos);

int i_test;
int j_test;
for (i_test=0,j_test=0;i_test<data.n_M*data.n_z*data.n_theta;i_test++) {if((data.mergers_array)[i_test]>20 && j_test<200) {fprintf(stderr,"%3d %d\n",i_test,(data.mergers_array)[i_test]);j_test++;}}
  fwrite(data.mergers_array,sizeof(int),data.n_M*data.n_z*data.n_theta,data.fp_out_halos);
  fclose(data.fp_out_halos);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  SID_free(SID_FARG data.halo_array);
  SID_free(SID_FARG data.mergers_array);
  if(SID.I_am_Master)
     SID_free(SID_FARG data.a_list);
  ADaPS_free(SID_FARG cosmo);
  SID_free(SID_FARG line);

  SID_exit(ERROR_NONE);
}

