#include <stdio.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

void init_treenode_trend_coordinate(tree_info *trees,trend_info *trend,const char *name){
  if(!strcmp(name,"z"))
     init_trend_coordinate(trend,name,trees,init_tree_property_z,free_tree_property_z,calc_tree_property_index_z);
  else if(!strcmp(name,"logM"))
     init_trend_coordinate(trend,name,trees,init_tree_property_logM,free_tree_property_logM,calc_tree_property_index_logM);
  else if(!strcmp(name,"xoff"))
     init_trend_coordinate(trend,name,trees,init_tree_property_xoff,free_tree_property_xoff,calc_tree_property_index_xoff);
  else if(!strcmp(name,"SSFctn"))
     init_trend_coordinate(trend,name,trees,init_tree_property_SSFctn,free_tree_property_SSFctn,calc_tree_property_index_SSFctn);
  else if(!strcmp(name,"tau_form"))
     init_trend_coordinate(trend,name,trees,init_tree_property_tau,free_tree_property_tau,calc_tree_property_index_tau_form);
  else if(!strcmp(name,"tau_3to1"))
     init_trend_coordinate(trend,name,trees,init_tree_property_tau,free_tree_property_tau,calc_tree_property_index_tau_3to1);
  else if(!strcmp(name,"tau_10to1"))
     init_trend_coordinate(trend,name,trees,init_tree_property_tau,free_tree_property_tau,calc_tree_property_index_tau_10to1);
  else
     SID_trap_error("Invalid property {%s} specified in init_treenode_trend().",ERROR_LOGIC,name);
}
