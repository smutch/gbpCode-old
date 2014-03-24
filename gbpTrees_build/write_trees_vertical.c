#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int write_forest_vertical_halos_recursive_local(tree_info *trees,tree_node_info *node,SID_fp *fp_out);
int write_forest_vertical_halos_recursive_local(tree_info *trees,tree_node_info *node,SID_fp *fp_out){
  int                      n_halos_written=0;

  // Fetch the halo properties 
  halo_properties_SAGE_info *halo;
  if(node->parent==NULL)
     halo=&(trees->group_properties_SAGE[node->snap_tree][node->neighbour_index]);
  else
     halo=&(trees->subgroup_properties_SAGE[node->snap_tree][node->neighbour_index]);

  // Perform write of halo properties
  SID_fwrite(halo,sizeof(halo_properties_SAGE_info),1,fp_out);

  // Increment the counter
  n_halos_written++;

  // Recurse over all progenitors.  Naturally leads
  //    to depth-first ordering.
  tree_node_info *current=node->progenitor_first;
  while(current!=NULL){
    n_halos_written+=write_forest_vertical_halos_recursive_local(trees,current,fp_out);
    current=current->progenitor_next;
  }
  return(n_halos_written);
}

void write_trees_vertical(tree_info     *trees,
                          double         box_size,
                          int            grid_size,
                          const char    *filename_root_trees_out){

  // Create a directory for the results
  char filename_root_out[MAX_FILENAME_LENGTH];
  sprintf(filename_root_out,"%s/vertical",filename_root_trees_out);
  mkdir(filename_root_out,02755);

  // Calculatet he total number of files to be written
  int n_files_write=grid_size*grid_size*grid_size;

  // Do groups first, then subgroups
  for(int i_type=0;i_type<2;i_type++){
     halo_properties_SAGE_info **properties;
     tree_node_info            **first_node_in_forest;
     char                        group_text_prefix[5];
     int                        *n_halos_forest_local;
     int                         n_forests_local;
     int                         n_forests;
     int                         n_halos_local;
     int                         n_halos;
     switch(i_type){
        case 0:
           sprintf(group_text_prefix,"");
           properties          =trees->group_properties_SAGE;
           first_node_in_forest=trees->first_in_forest_groups;
           n_halos_forest_local=trees->n_groups_forest_local;
           n_forests_local     =trees->n_forests_local;
           n_forests           =trees->n_forests;
           n_halos_local       =trees->n_groups_trees_local;
           n_halos             =trees->n_groups_trees;
           break;
        case 1:
           sprintf(group_text_prefix,"sub");
           properties          =trees->subgroup_properties_SAGE;
           first_node_in_forest=trees->first_in_forest_subgroups;
           n_halos_forest_local=trees->n_subgroups_forest_local;
           n_forests_local     =trees->n_forests_local;
           n_forests           =trees->n_forests;
           n_halos_local       =trees->n_subgroups_trees_local;
           n_halos             =trees->n_subgroups_trees;
           break;
     }
     SID_log("Writing %d %sgroup forests (structure size=%lld bytes)...",SID_LOG_OPEN|SID_LOG_TIMER,
                                                                         n_forests,group_text_prefix,sizeof(halo_properties_SAGE_info));

     // Allocate some need3ed arrays
     int *file_out    =(int *)SID_malloc(sizeof(int)*n_forests_local);
     int *halo_count  =(int *)SID_calloc(sizeof(int)*n_files_write);
     int *forest_count=(int *)SID_calloc(sizeof(int)*n_files_write);

     // Perform forest/halo counts
     SID_log("Performing forest and halo counts...",SID_LOG_OPEN);
     for(int i_forest=0;i_forest<n_forests_local;i_forest++){
        tree_node_info            *current_halo=first_node_in_forest[i_forest];
        halo_properties_SAGE_info *halo_properties=&(properties[current_halo->snap_tree][current_halo->neighbour_index]);
        // Decide which file this forest belongs to
        int i_x=(int)((float)grid_size*(halo_properties->pos[0]/box_size));i_x=MIN(i_x,grid_size-1);
        int i_y=(int)((float)grid_size*(halo_properties->pos[1]/box_size));i_y=MIN(i_y,grid_size-1);
        int i_z=(int)((float)grid_size*(halo_properties->pos[2]/box_size));i_z=MIN(i_z,grid_size-1);
        int i_g=(i_z*grid_size+i_y)*grid_size+i_x;
        forest_count[i_g]++;
        file_out[i_forest]=i_g;
        halo_count[i_g]+=n_halos_forest_local[i_forest];
     }
     SID_Allreduce(SID_IN_PLACE,halo_count,  n_files_write,SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,forest_count,n_files_write,SID_INT,SID_SUM,SID.COMM_WORLD);
     // Sanity checks
     int forest_count_total;
     int halo_count_total;
     calc_sum(forest_count,&forest_count_total,n_files_write,SID_INT,CALC_MODE_DEFAULT);
     calc_sum(halo_count,  &halo_count_total,  n_files_write,SID_INT,CALC_MODE_DEFAULT);
     if(forest_count_total!=n_forests)
        SID_trap_error("Invalid forest count (ie. %d!=%d)",ERROR_LOGIC,forest_count_total,n_forests);
     if(halo_count_total!=n_halos)
        SID_trap_error("Invalid halo count (ie. %d!=%d)",ERROR_LOGIC,halo_count_total,n_halos);
     SID_log("Done.",SID_LOG_CLOSE);

     // Write the file headers
     SID_log("Writing headers...",SID_LOG_OPEN);
     for(int i_rank=0;i_rank<SID.n_proc;i_rank++){
        // Each rank takes it's turn
        if(SID.My_rank==i_rank){
           // Open files
           SID_fp *fp_out;
           fp_out=(SID_fp *)SID_malloc(sizeof(SID_fp)*n_files_write);
           for(int i_file=0;i_file<n_files_write;i_file++){
              char   filename_out[MAX_FILENAME_LENGTH];
              sprintf(filename_out,"%s/%sgroup_forests_%03d.dat",filename_root_out,group_text_prefix,i_file);
              if(i_rank==0){
                 SID_fopen(filename_out,"w",&(fp_out[i_file]));
                 SID_fwrite(&(forest_count[i_file]),sizeof(int),1,&(fp_out[i_file]));
                 SID_fwrite(&(halo_count[i_file]),  sizeof(int),1,&(fp_out[i_file]));
              }
              else
                 SID_fopen(filename_out,"a",&(fp_out[i_file]));
           }
           // Write the number of halos per forest
           for(int i_forest=0;i_forest<n_forests_local;i_forest++)
              SID_fwrite(&(n_halos_forest_local[i_forest]),sizeof(int),1,&(fp_out[file_out[i_forest]]));
           // Close files
           for(int i_file=0;i_file<n_files_write;i_file++)
              SID_fclose(&(fp_out[i_file]));
           SID_free(SID_FARG fp_out);
        }
        SID_Barrier(SID.COMM_WORLD);
     }
     SID_log("Done.",SID_LOG_CLOSE);

     // Write the halos
     SID_log("Writing halos...",SID_LOG_OPEN);
     int n_halos_written  =0;
     int n_forests_written=0;
     for(int i_rank=0;i_rank<SID.n_proc;i_rank++){
        // Each rank takes it's turn
        if(SID.My_rank==i_rank){
           // Open files
           SID_fp *fp_out;
           fp_out=(SID_fp *)SID_malloc(sizeof(SID_fp)*n_files_write);
           for(int i_file=0;i_file<n_files_write;i_file++){
              char   filename_out[MAX_FILENAME_LENGTH];
              sprintf(filename_out,"%s/%sgroup_forests_%03d.dat",filename_root_out,group_text_prefix,i_file);
              SID_fopen(filename_out,"a",&(fp_out[i_file]));
           }
           // Write the halos
           for(int i_forest=0;i_forest<n_forests_local;i_forest++){
              int             n_halos_forest_written=0;
              tree_node_info *current=first_node_in_forest[i_forest];
              while(current!=NULL){
                 if(current->descendant==NULL)
                    n_halos_forest_written+=write_forest_vertical_halos_recursive_local(trees,current,&(fp_out[file_out[i_forest]]));
                 current=current->next_in_forest;
              }
              if(n_halos_forest_written!=n_halos_forest_local[i_forest])
                 SID_trap_error("Number of halos written is not right (i.e. %d!=%d) for forest #%d",ERROR_LOGIC,
                                n_halos_written,
                                n_halos_forest_local[i_forest],i_forest);
              n_halos_written  +=n_halos_forest_written;
              n_forests_written++;
           }
           // Close files
           for(int i_file=0;i_file<n_files_write;i_file++)
              SID_fclose(&(fp_out[i_file]));
           SID_free(SID_FARG fp_out);
        }
        SID_Barrier(SID.COMM_WORLD);
     }
     calc_sum_global(&n_halos_written,  &n_halos_written,  1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
     calc_sum_global(&n_forests_written,&n_forests_written,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
     SID_log("Done.",SID_LOG_CLOSE);

     // Clean-up
     SID_free(SID_FARG file_out);
     SID_free(SID_FARG halo_count);
     SID_free(SID_FARG forest_count);

     // Write some information to the log
     SID_log("Halos   written: %d",SID_LOG_COMMENT,n_halos_written);
     SID_log("Forests written: %d",SID_LOG_COMMENT,n_forests_written);
     SID_log("Done.",SID_LOG_CLOSE);
  } // i_type
}

