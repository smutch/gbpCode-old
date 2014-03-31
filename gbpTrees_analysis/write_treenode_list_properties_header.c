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

void write_treenode_list_properties_header(tree_info *trees,treenode_list_info *list,FILE *fp_props_out){
  tree_node_info **list_in         =list->list;
  int              n_list_in       =list->n_list;
  int              flag_groups_list=list->flag_groups_list;
  char            *catalog_name    =list->catalog_name;

  // Master Rank does all the writing
  if(SID.I_am_Master){
     // Write header for props file
     int i_column=0;
     for(int i_write=0;i_write<WRITE_TREENODE_LIST_PROPERTIES_N;i_write++){
        char write_name[128];
        write_treenode_list_properties_set_ith(trees,
                                               NULL,
                                               i_write,
                                               write_name,
                                               NULL,
                                               NULL,
                                               NULL);
        if(i_write==0){
           if(flag_groups_list)
              fprintf(fp_props_out,"# Properties for group catalog {%s}\n",catalog_name);
           else
              fprintf(fp_props_out,"# Properties for subgroup catalog {%s}\n",catalog_name);
           fprintf(fp_props_out,"#\n");
           fprintf(fp_props_out,"# Column (%02d): Catalog item number\n",i_column);i_column++;
        }
        fprintf(fp_props_out,"#        (%02d): %s\n",i_column,write_name);i_column++;
     }
  }
}

