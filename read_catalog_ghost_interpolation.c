#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void read_catalog_ghost_interpolation(tree_horizontal_ghost_group_info     *groups,
                                      halo_properties_info                **group_properties,
                                      int                                   n_groups,
                                      tree_horizontal_ghost_subgroup_info  *subgroups,
                                      halo_properties_info                **subgroup_properties,
                                      int                                   n_subgroups,
                                      char                                 *filename_cat_root_in,
                                      int                                   i_file,
                                      int                                   i_read,
                                      int                                   n_wrap){
   // Read subgroups first
   fp_catalog_info fp_properties;
   int             i_subgroup;
   fopen_catalog(filename_cat_root_in,
                 i_read,
                 READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES,
                 &fp_properties);
   if(fp_properties.n_halos_total!=n_subgroups)
      SID_trap_error("Mismatched subgroup counts (ie. %d!=%d) in read_catalog_ghost_interpolation().",ERROR_LOGIC,fp_properties.n_halos_total,n_subgroups);
   for(i_subgroup=0;i_subgroup<n_subgroups;i_subgroup++){
      if(subgroups[i_subgroup].file_offset>1 || check_mode_for_flag(subgroups[i_subgroup].type,TREE_CASE_FOUND)){
         if(subgroup_properties[i_subgroup]==NULL)
            subgroup_properties[i_subgroup]=(halo_properties_info *)SID_malloc(sizeof(halo_properties_info));
         fread_catalog_file(&fp_properties,NULL,subgroup_properties[i_subgroup],NULL,i_subgroup);
      }
      else
         SID_free(SID_FARG subgroup_properties[i_subgroup]);
   }
   fclose_catalog(&fp_properties);

   // Then read the groups
   int i_group;
   fopen_catalog(filename_cat_root_in,
                 i_read,
                 READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES,
                 &fp_properties);
   if(fp_properties.n_halos_total!=n_groups)
      SID_trap_error("Mismatched group counts (ie. %d!=%d) in read_catalog_ghost_interpolation().",ERROR_LOGIC,fp_properties.n_halos_total,n_groups);
   for(i_group=0;i_group<n_groups;i_group++){
      if(groups[i_group].file_offset>1 || check_mode_for_flag(groups[i_group].type,TREE_CASE_FOUND)){
         if(group_properties[i_group]==NULL)
            group_properties[i_group]=(halo_properties_info *)SID_malloc(sizeof(halo_properties_info));
         fread_catalog_file(&fp_properties,NULL,group_properties[i_group],NULL,i_group);
      }
      else
         SID_free(SID_FARG group_properties[i_group]);
   }
   fclose_catalog(&fp_properties);

}

