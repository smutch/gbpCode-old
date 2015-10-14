#include <stdio.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>

void init_halo_trend_coordinate(halo_trend_info *halo_trend_data,trend_info *trend,const char *name){
  if(!strcmp(name,"z"))
     init_trend_coordinate(trend,name,halo_trend_data,init_halo_trend_property_z,free_halo_trend_property_z,calc_halo_trend_property_index_z);
  else if(!strcmp(name,"logM_FoF"))
     init_trend_coordinate(trend,name,halo_trend_data,init_halo_trend_property_logM_FoF,free_halo_trend_property_logM_FoF,calc_halo_trend_property_index_logM_FoF);
  else if(!strcmp(name,"SSFctn"))
     init_trend_coordinate(trend,name,halo_trend_data,init_halo_trend_property_SSFctn,free_halo_trend_property_SSFctn,calc_halo_trend_property_index_SSFctn);
  else
     SID_trap_error("Invalid property {%s} specified in init_halo_trend().",ERROR_LOGIC,name);
}
