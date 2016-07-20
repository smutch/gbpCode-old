#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

void write_treenode_list_markers(tree_info *trees,const char *filename_out_root,treenode_list_info *list){
  tree_node_info **list_in         =list->list;
  int              n_list_in       =list->n_list_local;
  int              flag_groups_list=list->flag_groups_list;
  char            *catalog_name    =list->catalog_name;

  // Open file
  char  filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_out,"%s_%s_markers.txt",filename_out_root,list->catalog_name);
  FILE *fp_props_out=fopen(filename_out,"w");

  SID_log("Writing treenode_list markers to {%s}...",SID_LOG_OPEN,filename_out);

  // Write the header
  write_treenode_list_markers_header(trees,list,fp_props_out);

  // Create arrays holding halo_IDs and tree_case flags
  int n_list_max;
  calc_max_global(&n_list_in,&n_list_max,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
  int *halo_ID_list  =(int *)SID_malloc(sizeof(int)*n_list_max);
  int *tree_case_list=(int *)SID_malloc(sizeof(int)*n_list_max);

  int j_list=0;
  for(int i_rank=0;i_rank<SID.n_proc;i_rank++){
     if(i_rank==SID.My_rank){
        for(int i_list=0;i_list<n_list_in;i_list++){
           halo_ID_list[i_list]  =list_in[i_list]->halo_ID;
           tree_case_list[i_list]=list_in[i_list]->tree_case;
        }
     }
     SID_Bcast(halo_ID_list,  sizeof(int)*n_list_max,i_rank,SID.COMM_WORLD);
     SID_Bcast(tree_case_list,sizeof(int)*n_list_max,i_rank,SID.COMM_WORLD);
     // Master Rank does all the writing
     if(SID.My_rank==i_rank || SID.I_am_Master){
        int n_list_i;
        // Generate properties
        if(i_rank==0)
           n_list_i=n_list_in;
        else{
           SID_Sendrecv(&n_list_in,
                        1,
                        SID_INT,
                        MASTER_RANK,
                        1918270,
                        &n_list_i,
                        1,
                        SID_INT,
                        i_rank,
                        1918270,
                        SID.COMM_WORLD);
        }
        for(int i_list=0;i_list<n_list_i;i_list++,j_list++){
           // Point to the halo to be processed 
           tree_markers_info *markers;
           tree_node_info    *current_halo;
           if(i_rank==SID.My_rank){
              current_halo=list_in[i_list];
              // Find some special nodes for this listed halo
              markers=fetch_treenode_precomputed_markers(trees,current_halo);
           }

           // Write properties
           int n_write;
           if(SID.I_am_Master)
              fprintf(fp_props_out,"%4d %4d %4d",j_list,halo_ID_list[i_list],tree_case_list[i_list]);
           if(!flag_groups_list)
              n_write=12;
           else
              n_write=9; // Don't write the halos at the end which pertain only to subgroups
           for(int i_write=0;i_write<n_write;i_write++){
              // Compute properties
              int    i_z_node;
              int    idx_node;
              double t_node;
              double z_node;
              double M_node;
              double M_node_parent;
              int    n_p_node;
              if(SID.My_rank==i_rank){
                 tree_node_info *node_write=NULL;
                 char write_name[32];
                 switch(i_write){
                    case 0:
                       sprintf(write_name,"halo");
                       if(current_halo!=NULL) node_write=current_halo;
                       break;
                    case 1:
                       sprintf(write_name,"descendant");
                       if(current_halo!=NULL) node_write=current_halo->descendant;
                       break;
                    case 2:
                       sprintf(write_name,"main progenitor");
                       if(current_halo!=NULL) node_write=current_halo->progenitor_first;
                       break;
                    case 3:
                       sprintf(write_name,"peak_mass");
                       if(markers!=NULL) node_write=markers->peak_mass;
                       break;
                    case 4:
                       sprintf(write_name,"half_peak_mass");
                       if(markers!=NULL) node_write=markers->half_peak_mass;
                       break;
                    case 5:
                       sprintf(write_name,"root");
                       if(markers!=NULL) node_write=markers->branch_root;
                       break;
                    case 6:
                       sprintf(write_name,"leaf");
                       if(markers!=NULL) node_write=markers->branch_leaf;
                       break;
                    case 7:
                       sprintf(write_name,"merger_remnant_33pc");
                       if(markers!=NULL) node_write=markers->merger_33pc_remnant;
                       break;
                    case 8:
                       sprintf(write_name,"merger_remnant_10pc");
                       if(markers!=NULL) node_write=markers->merger_10pc_remnant;
                       break;
                    case 9:
                       sprintf(write_name,"parent");
                       if(current_halo!=NULL) node_write=current_halo->parent_top;
                       break;
                    case 10:
                       sprintf(write_name,"first_became_satellite");
                       if(markers!=NULL) node_write=markers->first_became_satellite;
                       break;
                    case 11:
                       sprintf(write_name,"joined_current_parent");
                       if(markers!=NULL) node_write=markers->joined_current_parent_top;
                       break;
                 }

                 if(node_write!=NULL){
                    i_z_node     =node_write->snap_tree;
                    idx_node     =node_write->neighbour_index;
                    t_node       =trees->t_list[i_z_node];
                    z_node       =trees->z_list[i_z_node];
                    M_node       =fetch_treenode_M_vir(trees,node_write);
                    M_node_parent=fetch_treenode_M_vir(trees,node_write->parent_top);
                    n_p_node     =fetch_treenode_n_particles(trees,node_write);
                 }
                 else{
                    i_z_node       =-1;
                    idx_node       =-1;
                    t_node         =-1.;
                    z_node         =-1.;
                    M_node         = 0.;
                    M_node_parent  = 0.;
                    n_p_node       = 0;
                 }
              }

              // Write properties
              if(i_rank!=0){
                 SID_Sendrecv(&i_z_node,     1,SID_INT,   MASTER_RANK,1918271,&i_z_node,     1,SID_INT,   i_rank,1918271,SID.COMM_WORLD);
                 SID_Sendrecv(&idx_node,     1,SID_INT,   MASTER_RANK,1918272,&idx_node,     1,SID_INT,   i_rank,1918272,SID.COMM_WORLD);
                 SID_Sendrecv(&t_node,       1,SID_DOUBLE,MASTER_RANK,1918273,&t_node,       1,SID_DOUBLE,i_rank,1918273,SID.COMM_WORLD);
                 SID_Sendrecv(&z_node,       1,SID_DOUBLE,MASTER_RANK,1918274,&z_node,       1,SID_DOUBLE,i_rank,1918274,SID.COMM_WORLD);
                 SID_Sendrecv(&M_node,       1,SID_DOUBLE,MASTER_RANK,1918275,&M_node,       1,SID_DOUBLE,i_rank,1918275,SID.COMM_WORLD);
                 SID_Sendrecv(&M_node_parent,1,SID_DOUBLE,MASTER_RANK,1918276,&M_node_parent,1,SID_DOUBLE,i_rank,1918276,SID.COMM_WORLD);
              }
              if(SID.I_am_Master){
                 int snap_node=-1;
                 if(i_z_node>=0)
                    snap_node=trees->snap_list[i_z_node];
                 fprintf(fp_props_out," %4d %7d %10.3le %5.2lf %6.3lf %7d",
                                                              snap_node,
                                                              idx_node,
                                                              t_node/S_PER_YEAR,
                                                              z_node,
                                                              take_log10(M_node),
                                                              n_p_node);
                 if(!flag_groups_list)
                    fprintf(fp_props_out," %6.3lf",take_log10(M_node_parent));
              }
           } // i_write
           if(SID.I_am_Master) fprintf(fp_props_out,"\n");
        }
     } // if i_rank
     SID_Barrier(SID.COMM_WORLD);
  } // for i_rank
  SID_free(SID_FARG halo_ID_list);
  SID_free(SID_FARG tree_case_list);
  fclose(fp_props_out);
  SID_log("Done.",ERROR_LOGIC);
}

