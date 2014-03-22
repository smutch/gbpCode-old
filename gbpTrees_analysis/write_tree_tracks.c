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

void write_tree_tracks(tree_info *trees,tree_node_info **list_in,int n_list_in,int mode,const char *catalog_name){
  SID_log("Processing %d halos in catalog {%s}...",SID_LOG_OPEN,n_list_in,catalog_name);

  // Fetch properties
  halo_properties_info **group_properties   =(halo_properties_info **)ADaPS_fetch(trees->data,"properties_groups");
  halo_properties_info **subgroup_properties=(halo_properties_info **)ADaPS_fetch(trees->data,"properties_subgroups");

  // Are we processing groups?
  int flag_processing_groups=FALSE;
  halo_properties_info **halo_properties;
  if(mode==1){
     flag_processing_groups=TRUE;
     halo_properties       =group_properties;
  }
  else{
     flag_processing_groups=FALSE;
     halo_properties       =subgroup_properties;
  }

  // Take the halo pointers back to their start
  int              n_list     =n_list_in;
  tree_node_info **list       =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_list_in);
  int             *track_index=(int             *)SID_malloc(sizeof(int)*n_list_in);
  for(int i_list=0;i_list<n_list;i_list++){
     track_index[i_list]=i_list; // default
     find_treenode_leaf(trees,list_in[i_list],&(list[i_list]));
  }

  // Count fragmented halos
  int n_frag=0;
  for(int i_list=0;i_list<n_list_in;i_list++){
    if(check_mode_for_flag(list_in[i_list]->tree_case,TREE_CASE_FRAGMENTED_STRAYED)  ||
       check_mode_for_flag(list_in[i_list]->tree_case,TREE_CASE_FRAGMENTED_EXCHANGED)||
       check_mode_for_flag(list_in[i_list]->tree_case,TREE_CASE_FRAGMENTED_RETURNED))
       n_frag++;
  }
  if(n_frag>0)
     SID_log("%d fragmented...",SID_LOG_CONTINUE,n_frag);

  // Remove duplicates (These are ugly N^2 algorythms.  Fix them sometime.)
  for(int i_list=0;i_list<n_list_in;i_list++)
     track_index[i_list]=i_list;
  for(int i_list=0;i_list<n_list;i_list++){
     for(int j_list=(i_list+1);j_list<n_list;j_list++){
        if(list[i_list]==list[j_list]){
           n_list--;
           for(int k_list=j_list;k_list<n_list;k_list++)
              list[k_list]=list[k_list+1];
           int old_index=track_index[j_list];
           for(int k_list=0;k_list<n_list;k_list++){
              if(track_index[k_list]>old_index)
                 track_index[k_list]--;
           }
           track_index[j_list]=track_index[i_list];
        }
     }
  }
  for(int i_list=0;i_list<n_list_in;i_list++){
     for(int j_list=0;j_list<n_list;j_list++){
        if(list_in[i_list]==list[j_list])
           track_index[i_list]=j_list;
     }
  }
  if(n_list!=n_list_in)
     SID_log("%d removed as duplicates...",SID_LOG_CONTINUE,n_list_in-n_list);

  // Master Rank does all the writing
  FILE *fp_tracks_out=NULL;
  FILE *fp_props_out =NULL;
  if(SID.I_am_Master){
     // Create and open the output files
     char   filename_tracks_out[MAX_FILENAME_LENGTH];
     char   filename_props_out[MAX_FILENAME_LENGTH];
     sprintf(filename_tracks_out,"%s_tracks.dat",catalog_name);
     sprintf(filename_props_out, "%s_props.txt", catalog_name);
     fp_tracks_out=fopen(filename_tracks_out,"w");
     fp_props_out =fopen(filename_props_out, "w");

     // Write header for tracks file
     fwrite(&n_list,          sizeof(int),   1,             fp_tracks_out);
     fwrite(&(trees->n_snaps),sizeof(int),   1,             fp_tracks_out);
     fwrite(trees->snap_list, sizeof(int),   trees->n_snaps,fp_tracks_out);
     fwrite(trees->z_list,    sizeof(double),trees->n_snaps,fp_tracks_out);
     fwrite(trees->t_list,    sizeof(double),trees->n_snaps,fp_tracks_out);

     // Write header for props file
     int n_write;
     int i_write;
     int i_column;
     if(!flag_processing_groups)
        n_write=6;
     else
        n_write=4; // Don't write the halos at the end which pertain only to subgroups
     for(i_write=0,i_column=1;i_write<n_write;i_write++){
        char write_name[32];
        switch(i_write){
           case 0:
              sprintf(write_name,"obs");
              break;
           case 1:
              sprintf(write_name,"main_progenitor");
              break;
           case 2:
              sprintf(write_name,"form");
              break;
           case 3:
              sprintf(write_name,"last");
              break;
           case 4:
              sprintf(write_name,"accrete_last");
              break;
           case 5:
              sprintf(write_name,"accrete_first");
              break;
        }
        if(i_write==0){
           if(flag_processing_groups)
              fprintf(fp_props_out,"# Properties for group catalog {%s}\n",catalog_name);
           else
              fprintf(fp_props_out,"# Properties for subgroup catalog {%s}\n",catalog_name);
           fprintf(fp_props_out,"#\n");
           fprintf(fp_props_out,"# Column (%02d): Catalog item number\n",i_column,write_name);i_column++;
           fprintf(fp_props_out,"#        (%02d): Track index\n",        i_column,write_name);i_column++;
        }
        fprintf(fp_props_out,"#        (%02d): Snapshot No. at t_%s\n",               i_column,write_name);i_column++;
        fprintf(fp_props_out,"#        (%02d): Index No.    at t_%s\n",               i_column,write_name);i_column++;
        fprintf(fp_props_out,"#        (%02d): t_%s\n",                               i_column,write_name);i_column++;
        fprintf(fp_props_out,"#        (%02d): z_%s\n",                               i_column,write_name);i_column++;
        fprintf(fp_props_out,"#        (%02d): log_10(M_%s(z=z_%s) [M_sol])\n",       i_column,write_name,write_name);i_column++;
        if(!flag_processing_groups){
        fprintf(fp_props_out,"#        (%02d): log_10(M_parent_%s(z=z_%s) [M_sol])\n",i_column,write_name,write_name);i_column++;
        }
        fprintf(fp_props_out,"#        (%02d): n_p(z=z_%s)\n",                        i_column,write_name);i_column++;
     }
  }

  // Allocate some temporary arrays for the tracks
  int    *i_z_track=(int    *)SID_malloc(sizeof(int)   *trees->n_snaps);
  int    *idx_track=(int    *)SID_malloc(sizeof(int)   *trees->n_snaps);
  double *M_track  =(double *)SID_malloc(sizeof(double)*trees->n_snaps);
  double *x_track  =(double *)SID_malloc(sizeof(double)*trees->n_snaps);
  double *y_track  =(double *)SID_malloc(sizeof(double)*trees->n_snaps);
  double *z_track  =(double *)SID_malloc(sizeof(double)*trees->n_snaps);
  double *vx_track =(double *)SID_malloc(sizeof(double)*trees->n_snaps);
  double *vy_track =(double *)SID_malloc(sizeof(double)*trees->n_snaps);
  double *vz_track =(double *)SID_malloc(sizeof(double)*trees->n_snaps);

  // Process z_obs halos
  int k_z_obs=0;
  for(int i_rank=0;i_rank<SID.n_proc;i_rank++){
     if(SID.My_rank==i_rank || SID.I_am_Master){
        int n_list_i;
        // Generate properties
        if(i_rank==0)
           n_list_i=n_list_in;
        else
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
        for(int i_list=0;i_list<n_list_i;i_list++){
           // Point to the halo to be processed 
           tree_node_info *current_halo=list_in[i_list];

           // Find some special nodes for this listed halo
           tree_node_info *descendant_last;
           tree_node_info *progenitor_main;
           tree_node_info *progenitor_formation;
           tree_node_info *progenitor_first_accretion;
           tree_node_info *progenitor_last_accretion;
           find_treenode_last_descendant(trees,current_halo,&descendant_last);
           find_treenode_main_progenitor(trees,current_halo,&progenitor_main);
           find_treenode_accretion(      trees,current_halo,&progenitor_first_accretion,&progenitor_last_accretion);
           find_treenode_formation(      trees,current_halo,0.5,&progenitor_formation);

           if(descendant_last->snap_tree==(trees->n_snaps-1))
              descendant_last=NULL;

           // Write properties
           int n_write;
           if(i_rank==0)
              fprintf(fp_props_out,"%4d %4d",i_list,track_index[i_list]);
           if(!flag_processing_groups)
              n_write=6;
           else
              n_write=4; // Don't write the halos at the end which pertain only to subgroups
           for(int i_write=0;i_write<n_write;i_write++){
              // Compute properties
              int    i_z_node;
              int    idx_node;
              int    idx_node_parent;
              double t_node;
              double z_node;
              double M_node;
              double M_node_parent;
              int    n_p_node;
              if(SID.My_rank==i_rank){
                 tree_node_info *node_write;
                 char write_name[32];
                 switch(i_write){
                    case 0:
                       sprintf(write_name,"obs");
                       node_write=current_halo;
                       break;
                    case 1:
                       sprintf(write_name,"main_progenitor");
                       node_write=progenitor_main;
                       break;
                    case 2:
                       sprintf(write_name,"form");
                       node_write=progenitor_formation;
                       break;
                    case 3:
                       sprintf(write_name,"last");
                       node_write=descendant_last;
                       break;
                    case 4:
                       sprintf(write_name,"accrete_last");
                       node_write=progenitor_last_accretion;
                       break;
                    case 5:
                       sprintf(write_name,"accrete_first");
                       node_write=progenitor_first_accretion;
                       break;
                 }

                 if(node_write!=NULL){
                    i_z_node=node_write->snap_tree;
                    idx_node=node_write->neighbour_index;
                    t_node  =trees->t_list[i_z_node];
                    z_node  =trees->z_list[i_z_node];
                    M_node  =halo_properties[i_z_node][idx_node].M_vir;
                    n_p_node=halo_properties[i_z_node][idx_node].n_particles;
                    if(node_write->parent!=NULL && !flag_processing_groups){
                       idx_node_parent=node_write->parent->snap_tree;
                       M_node_parent  =group_properties[i_z_node][idx_node_parent].M_vir;
                    }
                    else{
                       idx_node_parent=-1;
                       M_node_parent  = 0.;
                    }
                 }
                 else{
                    i_z_node=-1;
                    idx_node=-1;
                    t_node  =-1.;
                    z_node  =-1.;
                    M_node  = 0.;
                    n_p_node= 0;
                    idx_node_parent=-1;
                    M_node_parent  = 0.;
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
                 if(!flag_processing_groups)
                    fprintf(fp_props_out," %6.3lf",take_log10(M_node_parent));
              }
           } // i_write
           if(SID.I_am_Master) fprintf(fp_props_out,"\n");
        }

        // Generate tracks
        if(i_rank==0)
           n_list_i=n_list;
        else
           SID_Sendrecv(&n_list,
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
        for(int i_list=0;i_list<n_list_i;i_list++){
           // Point to the halo to be processed 
           tree_node_info *current_halo=list[i_list];
           // Compute track
           int n_track=0;
           if(SID.My_rank==i_rank){
              tree_node_info *current_track=current_halo;
              while(current_track!=NULL && check_treenode_if_main_progenitor(current_track)){
                 i_z_track[n_track]=current_track->snap_tree;
                 idx_track[n_track]=current_track->neighbour_index;
                 x_track[n_track]  =halo_properties[i_z_track[n_track]][idx_track[n_track]].position_MBP[0];
                 y_track[n_track]  =halo_properties[i_z_track[n_track]][idx_track[n_track]].position_MBP[1];
                 z_track[n_track]  =halo_properties[i_z_track[n_track]][idx_track[n_track]].position_MBP[2];
                 vx_track[n_track] =halo_properties[i_z_track[n_track]][idx_track[n_track]].velocity_COM[0];
                 vy_track[n_track] =halo_properties[i_z_track[n_track]][idx_track[n_track]].velocity_COM[1];
                 vz_track[n_track] =halo_properties[i_z_track[n_track]][idx_track[n_track]].velocity_COM[2];
                 M_track[n_track]  =halo_properties[i_z_track[n_track]][idx_track[n_track]].M_vir;
                 n_track++;
                 current_track=current_track->descendant;
              }
           }

           // Write track
           if(i_rank!=0){
              SID_Sendrecv(&n_track,       1,SID_INT,MASTER_RANK,1918370,&n_track,       1,SID_INT,i_rank,1918370,SID.COMM_WORLD);
              SID_Sendrecv(i_z_track,n_track,SID_INT,MASTER_RANK,1918371,i_z_track,n_track,SID_INT,i_rank,1918371,SID.COMM_WORLD);
              SID_Sendrecv(idx_track,n_track,SID_INT,MASTER_RANK,1918372,idx_track,n_track,SID_INT,i_rank,1918372,SID.COMM_WORLD);
              SID_Sendrecv(x_track,  n_track,SID_INT,MASTER_RANK,1918373,x_track,  n_track,SID_INT,i_rank,1918373,SID.COMM_WORLD);
              SID_Sendrecv(y_track,  n_track,SID_INT,MASTER_RANK,1918374,y_track,  n_track,SID_INT,i_rank,1918374,SID.COMM_WORLD);
              SID_Sendrecv(z_track,  n_track,SID_INT,MASTER_RANK,1918375,z_track,  n_track,SID_INT,i_rank,1918375,SID.COMM_WORLD);
              SID_Sendrecv(vx_track, n_track,SID_INT,MASTER_RANK,1918376,vx_track, n_track,SID_INT,i_rank,1918376,SID.COMM_WORLD);
              SID_Sendrecv(vy_track, n_track,SID_INT,MASTER_RANK,1918377,vy_track, n_track,SID_INT,i_rank,1918377,SID.COMM_WORLD);
              SID_Sendrecv(vz_track, n_track,SID_INT,MASTER_RANK,1918378,vz_track, n_track,SID_INT,i_rank,1918378,SID.COMM_WORLD);
              SID_Sendrecv(M_track,  n_track,SID_INT,MASTER_RANK,1918379,M_track,  n_track,SID_INT,i_rank,1918379,SID.COMM_WORLD);
           }
           if(SID.I_am_Master){
              fwrite(&n_track, sizeof(int),   1,      fp_tracks_out);
              fwrite(i_z_track,sizeof(int),   n_track,fp_tracks_out);
              fwrite(idx_track,sizeof(int),   n_track,fp_tracks_out);
              fwrite(x_track,  sizeof(double),n_track,fp_tracks_out);
              fwrite(y_track,  sizeof(double),n_track,fp_tracks_out);
              fwrite(z_track,  sizeof(double),n_track,fp_tracks_out);
              fwrite(vx_track, sizeof(double),n_track,fp_tracks_out);
              fwrite(vy_track, sizeof(double),n_track,fp_tracks_out);
              fwrite(vz_track, sizeof(double),n_track,fp_tracks_out);
              fwrite(M_track,  sizeof(double),n_track,fp_tracks_out);
           }

        } // for i_list
     } // if i_rank
     SID_Barrier(SID.COMM_WORLD);
  } // for i_rank
  if(SID.I_am_Master){
     fclose(fp_tracks_out);
     fclose(fp_props_out);
  }

  // Clean-up
  SID_free(SID_FARG track_index);
  SID_free(SID_FARG list);
  SID_free(SID_FARG i_z_track);
  SID_free(SID_FARG idx_track);
  SID_free(SID_FARG M_track);
  SID_free(SID_FARG x_track);
  SID_free(SID_FARG y_track);
  SID_free(SID_FARG z_track);
  SID_free(SID_FARG vx_track);
  SID_free(SID_FARG vy_track);
  SID_free(SID_FARG vz_track);

  SID_log("Done.",SID_LOG_CLOSE);
}

