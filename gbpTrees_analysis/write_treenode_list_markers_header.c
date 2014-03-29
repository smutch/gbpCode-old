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

void write_treenode_list_markers_header(tree_info *trees,treenode_list_info *list,FILE *fp_props_out){
  tree_node_info **list_in         =list->list;
  int              n_list_in       =list->n_list;
  int              flag_groups_list=list->flag_groups_list;
  char            *catalog_name    =list->catalog_name;

  // Master Rank does all the writing
  if(SID.I_am_Master){
     // Write header for props file
     int n_write;
     int i_write;
     int i_column;
     if(!flag_groups_list)
        n_write=9;
     else
        n_write=6; // Don't write the halos at the end which pertain only to subgroups
     for(i_write=0,i_column=1;i_write<n_write;i_write++){
        char write_name[32];
        switch(i_write){
           case 0:
              sprintf(write_name,"halo");
              break;
           case 1:
              sprintf(write_name,"main_progenitor");
              break;
           case 2:
              sprintf(write_name,"peak_mass");
              break;
           case 3:
              sprintf(write_name,"half_peak_mass");
              break;
           case 4:
              sprintf(write_name,"root");
              break;
           case 5:
              sprintf(write_name,"leaf");
              break;
           case 6:
              sprintf(write_name,"parent");
              break;
           case 7:
              sprintf(write_name,"accrete_last");
              break;
           case 8:
              sprintf(write_name,"accrete_first");
              break;
        }

        if(i_write==0){
           if(flag_groups_list)
              fprintf(fp_props_out,"# Properties for group catalog {%s}\n",catalog_name);
           else
              fprintf(fp_props_out,"# Properties for subgroup catalog {%s}\n",catalog_name);
           fprintf(fp_props_out,"#\n");
           fprintf(fp_props_out,"# Column (%02d): Catalog item number\n",      i_column,write_name);i_column++;
           fprintf(fp_props_out,"#        (%02d): Halo ID\n",                  i_column,write_name);i_column++;
           fprintf(fp_props_out,"#        (%02d): Tree case bit-wise-switch\n",i_column,write_name);i_column++;
        }
        fprintf(fp_props_out,"#        (%02d): Snapshot No. at t_%s\n",               i_column,write_name);i_column++;
        fprintf(fp_props_out,"#        (%02d): Index No.    at t_%s\n",               i_column,write_name);i_column++;
        fprintf(fp_props_out,"#        (%02d): t_%s\n",                               i_column,write_name);i_column++;
        fprintf(fp_props_out,"#        (%02d): z_%s\n",                               i_column,write_name);i_column++;
        fprintf(fp_props_out,"#        (%02d): log_10(M_%s(z=z_%s) [M_sol])\n",       i_column,write_name,write_name);i_column++;
        if(!flag_groups_list){
        fprintf(fp_props_out,"#        (%02d): log_10(M_parent_%s(z=z_%s) [M_sol])\n",i_column,write_name,write_name);i_column++;
        }
        fprintf(fp_props_out,"#        (%02d): n_p(z=z_%s)\n",                        i_column,write_name);i_column++;
     }
  }
}

