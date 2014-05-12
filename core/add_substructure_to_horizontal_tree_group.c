#include <gbpLib.h>
#include <gbpTrees_build.h>

void add_substructure_to_horizontal_tree_group(tree_horizontal_ghost_group_info    *group,
                                               tree_horizontal_ghost_subgroup_info *subgroup_descendant,
                                               tree_horizontal_ghost_subgroup_info *subgroup){
   tree_horizontal_extended_info *current;

   // Create a linked list of substructures for each group
   if(group->last_substructure!=NULL)
      group->last_substructure->next_substructure=subgroup;
   else
      group->first_substructure=subgroup;
   group->last_substructure   =subgroup;
   subgroup->next_substructure=NULL;

   // Increment the number of substructures the group has (needed when we write the tree file)
   group->n_subgroups++;

   // This pointer is needed to fetch the file_index value just before we perform the write
   subgroup->descendant=subgroup_descendant;
}

