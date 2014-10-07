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

void write_treenode_list_data_header(tree_info *trees,treenode_list_info *list,FILE *fp_props_out){
  tree_node_info **list_in         =list->list;
  int              flag_groups_list=list->flag_groups_list;
  char            *catalog_name    =list->catalog_name;

  // Master Rank does all the writing
  if(SID.I_am_Master){
     // Write header for props file
     int i_write =0;
     int i_column=1;
     ADaPS *current=list->data;
     while(current!=NULL){
        if(i_write==0){
           if(flag_groups_list)
              fprintf(fp_props_out,"# Properties for group catalog {%s}\n",catalog_name);
           else
              fprintf(fp_props_out,"# Properties for subgroup catalog {%s}\n",catalog_name);
           fprintf(fp_props_out,"#\n");
           fprintf(fp_props_out,"# Column (%02d): Catalog item number\n",i_column);i_column++;
        }
        fprintf(fp_props_out,"#        (%02d): %s\n",i_column,current->name);i_column++;
        i_write++;
        current=current->next;
     }
  }
}

