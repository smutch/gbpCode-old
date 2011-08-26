#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>
#include <gbpCosmo.h>

// millennium halo input structure 
typedef struct mil_halo_info mil_halo_info;
struct mil_halo_info {
  // merger tree pointers 
  int descendant;
  int progenitor_first;
  int progenitor_next;
  int FirstHaloInFOFgroup;
  int NextHaloInFOFgroup;

  // properties of halo 
  int   Len;
  float M_Mean200, M_vir, M_TopHat;  // Mean 200 values (Mvir=M_Crit200)
  float pos[3];
  float vel[3];
  float VelDisp;
  float Vmax;
  float Spin[3];
  long long MostBoundID;

  // original position in subfind output 
  int snap_num;
  int FileNr;
  int SubhaloIndex;
  float SubHalfMass;
};

int main(int argc, char *argv[]){
  char        filename_tree_in[MAX_FILENAME_LENGTH];
  char        filename_snap_in[MAX_FILENAME_LENGTH];
  char       *line=NULL;
  size_t      line_length=0;
  int         n_trees;
  int         n_halos_total;
  int        *n_halos;
  int         i_tree;
  FILE       *fp;
  void         *halos;
  halo_info     halo;
  size_t        halo_size;
  int        *snap_num;
  size_t     *snap_num_index;
  int         i_snap,i_halo,j_halo,k_halo;
  int         n_halos_snap;
  int        *group_halo_first;
  int         group_halo_last;
  size_t     *group_halo_first_index;
  int        *snap_index;
  int descendant_min,descendant_max;
  int progenitor_first_min,progenitor_first_max;
  int progenitor_next_min,progenitor_next_max;
  int group_halo_first_min,group_halo_first_max;
  int group_halo_next_min,group_halo_next_max;
  int snap_num_min,snap_num_max;
  int halo_index_min,halo_index_max;
  int n_gal=0;
  int max_snap=0;
  int n_halos_max;
  int n_subtrees;
  double statistic;
  int    depth_first_index;
  int    current;
  int     n_interp;
  double *x_int;
  double *y_int;
  interp_info *interp;
  cosmo_info *cosmo;
  FILE       *fp_in;

  int     n_snaps; 
  double *a_list;
  double *t_a_list;
  int     n_int;
  double  a_min;
  double  a_max;
  double  da;
  double  a;
  int     i_int;
  double *M_list;
  double *x_list;
  double *y_list;
  double *z_list;
  double *vx_list;
  double *vy_list;
  double *vz_list;
  double *t_list;
  double  a_descendant;
  double  t_descendant;
  double  M_vir_descendant;
  double  x_descendant;
  double  y_descendant;
  double  z_descendant;
  double  vx_descendant;
  double  vy_descendant;
  double  vz_descendant;
  int     progenitor;
  int     i_orbit;
  int     n_orbit;
  double  ratio_merge;
  double  v_merge;
  int     flag_progenitor_first;

  int          *branch;
  int          *branch_count;
  int           i_branch;
  int           j_branch;
  int           n_branches;
  double      **x_branch;
  double      **y_branch;
  double      **z_branch;
  double      **vx_branch;
  double      **vy_branch;
  double      **vz_branch;
  double      **t_branch;
  interp_info **x_interp;
  interp_info **y_interp;
  interp_info **z_interp;
  interp_info **vx_interp;
  interp_info **vy_interp;
  interp_info **vz_interp;
  int          *branch_root;
  int          *branch_desc;
  double      **a_branch;
  double      **M_branch;
  interp_info **M_interp;
  double        M_desc;
  int           flag_mil;


  SID_init(&argc,&argv,NULL);

  init_cosmo_std(&cosmo);

  // Fetch user inputs
  strcpy(filename_tree_in,argv[1]);
  strcpy(filename_snap_in,argv[2]);
  SID_log("Computing tree stats for {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_tree_in);
  if(argc>3){
    halo_size=sizeof(mil_halo_info);
    flag_mil=TRUE;
    SID_log("Using the Millennium data structure.",SID_LOG_COMMENT);
  }
  else{
    halo_size=sizeof(halo_info);
    flag_mil=FALSE;
    SID_log("Using the GiggleZ data structure.",SID_LOG_COMMENT);
  }

  // Read snapshot expansion factor list
  SID_log("Reading snapshot list {%s}...",SID_LOG_OPEN,filename_snap_in);
  fp_in   =fopen(filename_snap_in, "r");
  n_snaps =count_lines_data(fp_in);
  a_list  =(double *)SID_malloc(sizeof(double)*n_snaps);
  t_a_list=(double *)SID_malloc(sizeof(double)*n_snaps);
  for(i_snap=0;i_snap<n_snaps;i_snap++){
    grab_next_line_data(fp_in,&line,&line_length);
    grab_double(line,1,&(a_list[i_snap]));
  }
  fclose(fp_in);
  SID_log("Done.",SID_LOG_CLOSE);

  // Computing t(a)
  SID_log("Computing t(a) look-up table...",SID_LOG_OPEN);
  n_int=250;
  x_int=(double *)SID_malloc(sizeof(double)*n_int);
  y_int=(double *)SID_malloc(sizeof(double)*n_int);
  a_min=1e-5;
  a_max=1.;
  da   =(a_max-a_min)/((double)(n_int-2));
  x_int[0]=0.;
  y_int[0]=0.;
  for(i_int=1,a=a_min;i_int<n_int;i_int++,a+=da){
    x_int[i_int]=a;
    y_int[i_int]=1./(a*H_convert(H_z(z_of_a(a),cosmo))*S_PER_YEAR);
  }
  init_interpolate(x_int,y_int,(size_t)n_int,gsl_interp_cspline,&interp);
  for(i_snap=0;i_snap<n_snaps;i_snap++){
    t_a_list[i_snap]=interpolate_integral(interp,0.,a_list[i_snap]);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  t_list =(double *)SID_malloc(sizeof(double)*n_snaps);
  M_list =(double *)SID_malloc(sizeof(double)*n_snaps);
  x_list =(double *)SID_malloc(sizeof(double)*n_snaps);
  y_list =(double *)SID_malloc(sizeof(double)*n_snaps);
  z_list =(double *)SID_malloc(sizeof(double)*n_snaps);
  vx_list=(double *)SID_malloc(sizeof(double)*n_snaps);
  vy_list=(double *)SID_malloc(sizeof(double)*n_snaps);
  vz_list=(double *)SID_malloc(sizeof(double)*n_snaps);

  SID_log("Processing trees...",SID_LOG_OPEN|SID_LOG_TIMER);
  fp=fopen(filename_tree_in,"r");
  fread(&n_trees,      sizeof(int),1,fp);
  fread(&n_halos_total,sizeof(int),1,fp);
  SID_log("# of trees: %d",SID_LOG_COMMENT,n_trees);
  SID_log("# of halos: %d",SID_LOG_COMMENT,n_halos_total);
  n_halos=(int *)SID_malloc(sizeof(int)*n_trees);
  fread(n_halos,sizeof(int),n_trees,fp);
  calc_max(n_halos,&n_halos_max,n_trees,SID_INT,CALC_MODE_DEFAULT);
  SID_log("max_size:   %d",SID_LOG_COMMENT,n_halos_max);
  halos               =             SID_malloc(halo_size*n_halos_max);
  snap_num            =(int       *)SID_malloc(sizeof(int)*n_halos_max);
  snap_index          =(int       *)SID_malloc(sizeof(int)*n_halos_max);
  group_halo_first    =(int       *)SID_malloc(sizeof(int)*n_halos_max);
  branch              =(int       *)SID_malloc(sizeof(int)*n_halos_max);
//n_trees=1;
  for(i_tree=0;i_tree<n_trees;i_tree++){

    // Read this tree's halos
    fread(halos,halo_size,n_halos[i_tree],fp);

    // Count the branches
    depth_first_index=0;
    n_branches       =0;
    if(flag_mil){
      while(depth_first_index<n_halos[i_tree]){
        if(((mil_halo_info *)halos)[depth_first_index].progenitor_first<0)
          n_branches++;
        depth_first_index++;
      }
    }
    else{
      while(depth_first_index<n_halos[i_tree]){
        if(((halo_info *)halos)[depth_first_index].progenitor_first<0)
          n_branches++;
        depth_first_index++;
      }
    }

    // Count the halos in each branch
    branch_count=(int *)SID_malloc(sizeof(int)*n_branches);
    for(i_branch=0;i_branch<n_branches;i_branch++)
      branch_count[i_branch]=0;
    depth_first_index=0;
    i_branch         =0;
    if(flag_mil){
      while(depth_first_index<n_halos[i_tree]){
        branch[depth_first_index]=i_branch;
        branch_count[i_branch]++;
        if(((mil_halo_info *)halos)[depth_first_index].progenitor_first<0)
          i_branch++;
        depth_first_index++;
      }
    }
    else{
      while(depth_first_index<n_halos[i_tree]){
        branch[depth_first_index]=i_branch;
        branch_count[i_branch]++;
        if(((halo_info *)halos)[depth_first_index].progenitor_first<0)
          i_branch++;
        depth_first_index++;
      }
    }

    // Initialize branch orbit arrays
    branch_root =(int     *)SID_malloc(sizeof(int)     *n_branches);
    branch_desc =(int     *)SID_malloc(sizeof(int)     *n_branches);
    a_branch    =(double **)SID_malloc(sizeof(double *)*n_branches);
    t_branch    =(double **)SID_malloc(sizeof(double *)*n_branches);
    M_branch    =(double **)SID_malloc(sizeof(double *)*n_branches);
    x_branch    =(double **)SID_malloc(sizeof(double *)*n_branches);
    y_branch    =(double **)SID_malloc(sizeof(double *)*n_branches);
    z_branch    =(double **)SID_malloc(sizeof(double *)*n_branches);
    vx_branch   =(double **)SID_malloc(sizeof(double *)*n_branches);
    vy_branch   =(double **)SID_malloc(sizeof(double *)*n_branches);
    vz_branch   =(double **)SID_malloc(sizeof(double *)*n_branches);
    for(i_branch=0;i_branch<n_branches;i_branch++){
      a_branch[i_branch] =(double *)SID_malloc(sizeof(double)*branch_count[i_branch]);
      t_branch[i_branch] =(double *)SID_malloc(sizeof(double)*branch_count[i_branch]);
      M_branch[i_branch] =(double *)SID_malloc(sizeof(double)*branch_count[i_branch]);
      x_branch[i_branch] =(double *)SID_malloc(sizeof(double)*branch_count[i_branch]);
      y_branch[i_branch] =(double *)SID_malloc(sizeof(double)*branch_count[i_branch]);
      z_branch[i_branch] =(double *)SID_malloc(sizeof(double)*branch_count[i_branch]);
      vx_branch[i_branch]=(double *)SID_malloc(sizeof(double)*branch_count[i_branch]);
      vy_branch[i_branch]=(double *)SID_malloc(sizeof(double)*branch_count[i_branch]);
      vz_branch[i_branch]=(double *)SID_malloc(sizeof(double)*branch_count[i_branch]);
    }

    // Populate branch orbit arrays
    for(i_branch=0,depth_first_index=0;i_branch<n_branches;i_branch++){
      branch_root[i_branch]=depth_first_index;
      if(flag_mil){
        if(i_tree==0)
          branch_desc[i_branch]=0;
        else
          branch_desc[i_branch]=branch[((mil_halo_info *)halos)[depth_first_index].descendant];
        for(j_branch=0;j_branch<branch_count[i_branch];j_branch++,depth_first_index++){
          a_branch[i_branch][j_branch] =a_list[((mil_halo_info *)halos)[depth_first_index].snap_num];
          t_branch[i_branch][j_branch] =t_a_list[((mil_halo_info *)halos)[depth_first_index].snap_num];
          M_branch[i_branch][j_branch] =((mil_halo_info *)halos)[depth_first_index].M_vir;
          x_branch[i_branch][j_branch] =((mil_halo_info *)halos)[depth_first_index].pos[0];
          y_branch[i_branch][j_branch] =((mil_halo_info *)halos)[depth_first_index].pos[1];
          z_branch[i_branch][j_branch] =((mil_halo_info *)halos)[depth_first_index].pos[2];
          vx_branch[i_branch][j_branch]=((mil_halo_info *)halos)[depth_first_index].vel[0];
          vy_branch[i_branch][j_branch]=((mil_halo_info *)halos)[depth_first_index].vel[1];
          vz_branch[i_branch][j_branch]=((mil_halo_info *)halos)[depth_first_index].vel[2];
        }
      }
      else{
        if(i_tree==0)
          branch_desc[i_branch]=0;
        else
          branch_desc[i_branch]=branch[((halo_info *)halos)[depth_first_index].descendant];
        for(j_branch=0;j_branch<branch_count[i_branch];j_branch++,depth_first_index++){
          a_branch[i_branch][j_branch] =a_list[((halo_info *)halos)[depth_first_index].snap_num];
          t_branch[i_branch][j_branch] =t_a_list[((halo_info *)halos)[depth_first_index].snap_num];
          M_branch[i_branch][j_branch] =((halo_info *)halos)[depth_first_index].M_vir;
          x_branch[i_branch][j_branch] =((halo_info *)halos)[depth_first_index].pos[0];
          y_branch[i_branch][j_branch] =((halo_info *)halos)[depth_first_index].pos[1];
          z_branch[i_branch][j_branch] =((halo_info *)halos)[depth_first_index].pos[2];
          vx_branch[i_branch][j_branch]=((halo_info *)halos)[depth_first_index].vel[0];
          vy_branch[i_branch][j_branch]=((halo_info *)halos)[depth_first_index].vel[1];
          vz_branch[i_branch][j_branch]=((halo_info *)halos)[depth_first_index].vel[2];
        }
      }
    }
    if(depth_first_index!=n_halos[i_tree])
      SID_trap_error("Halo counts don't match (i.e. %d=%d)",ERROR_LOGIC,depth_first_index,n_halos[i_tree]);

    // Construct branch interpolations
    M_interp =(interp_info **)SID_malloc(sizeof(interp_info *)*n_branches);
    x_interp =(interp_info **)SID_malloc(sizeof(interp_info *)*n_branches);
    y_interp =(interp_info **)SID_malloc(sizeof(interp_info *)*n_branches);
    z_interp =(interp_info **)SID_malloc(sizeof(interp_info *)*n_branches);
    vx_interp=(interp_info **)SID_malloc(sizeof(interp_info *)*n_branches);
    vy_interp=(interp_info **)SID_malloc(sizeof(interp_info *)*n_branches);
    vz_interp=(interp_info **)SID_malloc(sizeof(interp_info *)*n_branches);
    for(i_branch=0;i_branch<n_branches;i_branch++){
      if(branch_count[i_branch]>4){
        init_interpolate(a_branch[i_branch],M_branch[i_branch], (size_t)branch_count[i_branch],gsl_interp_cspline,&(M_interp[i_branch]));
        init_interpolate(a_branch[i_branch],x_branch[i_branch], (size_t)branch_count[i_branch],gsl_interp_cspline,&(x_interp[i_branch]));
        init_interpolate(a_branch[i_branch],y_branch[i_branch], (size_t)branch_count[i_branch],gsl_interp_cspline,&(y_interp[i_branch]));
        init_interpolate(a_branch[i_branch],z_branch[i_branch], (size_t)branch_count[i_branch],gsl_interp_cspline,&(z_interp[i_branch]));
        init_interpolate(a_branch[i_branch],vx_branch[i_branch],(size_t)branch_count[i_branch],gsl_interp_cspline,&(vx_interp[i_branch]));
        init_interpolate(a_branch[i_branch],vy_branch[i_branch],(size_t)branch_count[i_branch],gsl_interp_cspline,&(vy_interp[i_branch]));
        init_interpolate(a_branch[i_branch],vz_branch[i_branch],(size_t)branch_count[i_branch],gsl_interp_cspline,&(vz_interp[i_branch]));
      }
    }

    // Analyze branches
    if(flag_mil){
      if(n_branches>1){
        for(i_branch=0;i_branch<n_branches;i_branch++){
          if(branch_count[branch_desc[i_branch]]>4){
            M_desc=interpolate(M_interp[branch_desc[i_branch]],a_branch[i_branch][0]);
            printf("%le ",a_list[((mil_halo_info *)halos)[branch_root[i_branch]].snap_num]);
            printf("%le ",t_branch[i_branch][0]);
            printf("%le ",M_desc);
            printf("%le ",M_desc/((mil_halo_info *)halos)[branch_root[i_branch]].M_vir);
            printf("\n");
          }
        }
      }
    }
    else{
      if(n_branches>1){
        for(i_branch=0;i_branch<n_branches;i_branch++){
          if(branch_count[branch_desc[i_branch]]>4){
            M_desc=interpolate(M_interp[branch_desc[i_branch]],a_branch[i_branch][0]);
            printf("%le ",a_list[((halo_info *)halos)[branch_root[i_branch]].snap_num]);
            printf("%le ",t_branch[i_branch][0]);
            printf("%le ",M_desc);
            printf("%le ",M_desc/((halo_info *)halos)[branch_root[i_branch]].M_vir);
            printf("\n");
          }
        }
      }
    }

    // Clean-up
    for(i_branch=0;i_branch<n_branches;i_branch++){
      SID_free(SID_FARG a_branch[i_branch]);
      SID_free(SID_FARG M_branch[i_branch]);
      SID_free(SID_FARG t_branch[i_branch]);
      SID_free(SID_FARG x_branch[i_branch]);
      SID_free(SID_FARG y_branch[i_branch]);
      SID_free(SID_FARG z_branch[i_branch]);
      SID_free(SID_FARG vx_branch[i_branch]);
      SID_free(SID_FARG vy_branch[i_branch]);
      SID_free(SID_FARG vz_branch[i_branch]);
    }
    SID_free(SID_FARG a_branch);
    SID_free(SID_FARG t_branch);
    SID_free(SID_FARG M_branch);
    SID_free(SID_FARG x_branch);
    SID_free(SID_FARG y_branch);
    SID_free(SID_FARG z_branch);
    SID_free(SID_FARG vx_branch);
    SID_free(SID_FARG vy_branch);
    SID_free(SID_FARG vz_branch);
    SID_free(SID_FARG branch_root);
    SID_free(SID_FARG branch_desc);
    for(i_branch=0;i_branch<n_branches;i_branch++){
      if(branch_count[i_branch]>4){
        free_interpolate(SID_FARG x_interp[i_branch]);
        free_interpolate(SID_FARG y_interp[i_branch]);
        free_interpolate(SID_FARG z_interp[i_branch]);
        free_interpolate(SID_FARG vx_interp[i_branch]);
        free_interpolate(SID_FARG vy_interp[i_branch]);
        free_interpolate(SID_FARG vz_interp[i_branch]);
      }
    }
    SID_free(SID_FARG x_interp);
    SID_free(SID_FARG y_interp);
    SID_free(SID_FARG z_interp);
    SID_free(SID_FARG vx_interp);
    SID_free(SID_FARG vy_interp);
    SID_free(SID_FARG vz_interp);
    SID_free(SID_FARG branch_count);
  }
  fclose(fp);
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  SID_log("Cleaning-up...",SID_LOG_OPEN);
  SID_free(SID_FARG t_a_list);
  SID_free(SID_FARG a_list);
  SID_free(SID_FARG t_list);
  SID_free(SID_FARG M_list);
  SID_free(SID_FARG x_list);
  SID_free(SID_FARG y_list);
  SID_free(SID_FARG z_list);
  SID_free(SID_FARG vx_list);
  SID_free(SID_FARG vy_list);
  SID_free(SID_FARG vz_list);

  SID_free((void **)&snap_num);
  SID_free((void **)&snap_index);
  SID_free((void **)&group_halo_first);
  SID_free((void **)&n_halos);
  SID_free((void **)&halos);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
