#include <gbpLib.h>
#include <gbpTrees.h>

void process_ghosts(tree_horizontal_ghost_group_info    **groups_in,
                    tree_horizontal_ghost_subgroup_info **subgroups_in,
                    int         *n_groups,    
                    int         *n_subgroups, 
                    int        **n_subgroups_group, 
                    int         *n_group_ghosts,
                    int         *n_subgroup_ghosts,
                    int         *n_group_ghosts_used,
                    int         *n_subgroup_ghosts_used,
                    int          i_file,
                    int          i_read,
                    int          j_file,
                    int          i_file_start,
                    int          n_search,
                    int          n_wrap,
                    int          n_files,
                    double      *a_list,
                    cosmo_info **cosmo){
   int file;
   int i_offset;
   int j_offset;
   int ghost_index;
   int i_group;
   int i_subgroup;
   int j_subgroup;

   // If n_group_ghosts & n_subgroup_ghosts are specified, then
   //    we are counting ghosts, n_group_ghosts_used & n_subgroup_ghosts_used
   //    are specified, then we are creating ghosts.
   int flag_create_ghosts;
   if(n_group_ghosts==NULL && n_subgroup_ghosts==NULL)
      flag_create_ghosts=TRUE;
   else if(n_group_ghosts_used==NULL && n_subgroup_ghosts_used==NULL)
      flag_create_ghosts=FALSE;
   else
      SID_trap_error("process_ghosts() has been called with ambiguous parameters.",ERROR_LOGIC);

   // Create an array that will hold the descendant line of each 
   //   group setting where subgroup ghosts will be placed.
   tree_horizontal_ghost_group_info **group_descendants;
   group_descendants=(tree_horizontal_ghost_group_info **)SID_malloc(sizeof(tree_horizontal_ghost_group_info *)*n_search);

   // Process one group at a time
   tree_horizontal_ghost_group_info group_dummy; // used when this function is called for ghost counting
   for(i_group=0,i_subgroup=0;i_group<n_groups[j_file];i_group++){
      tree_horizontal_ghost_group_info *group_descendant=NULL;
      int group_final_file_index;
      int group_file_offset;
      int group_ghost_base_type;
      int flag_descenadnt_is_new_ghost=FALSE;
      int flag_ghost_merger           =FALSE;

      // Loop over the search range, building group ghost-chains (if needed) and the descendant history of this group
      group_descendant =&(groups_in[i_file%n_wrap][i_group]);
      group_file_offset=groups_in[i_file%n_wrap][i_group].file_offset;
      flag_ghost_merger=check_mode_for_flag(group_descendant->type,TREE_CASE_MERGER);
      for(i_offset=1,file=i_file+1;i_offset<n_search && i_offset<group_file_offset && file<=i_file_start;i_offset++,file++){

         if(i_offset==1){
            group_final_file_index=group_descendant->file_index;
            group_ghost_base_type =TREE_CASE_GHOST;
            // Set the type for the ghost halos.  A merger flag may be added later for the last one.
            if(check_mode_for_flag(group_descendant->type,TREE_CASE_FRAGMENTED_LOST))
               group_ghost_base_type|=TREE_CASE_FRAGMENTED_LOST;
            if(check_mode_for_flag(group_descendant->type,TREE_CASE_FRAGMENTED_RETURNED))
               group_ghost_base_type|=TREE_CASE_FRAGMENTED_RETURNED;
            if(check_mode_for_flag(group_descendant->type,TREE_CASE_FRAGMENTED_EXCHANGED))
               group_ghost_base_type|=TREE_CASE_FRAGMENTED_EXCHANGED;
            if(check_mode_for_flag(group_descendant->type,TREE_CASE_SPUTTERED))
               group_ghost_base_type|=TREE_CASE_SPUTTERED;
            if(check_mode_for_flag(group_descendant->type,TREE_CASE_MAIN_PROGENITOR))
               group_ghost_base_type|=TREE_CASE_MAIN_PROGENITOR;
            if(check_mode_for_flag(group_descendant->type,TREE_CASE_MERGER)){
               flag_ghost_merger=TRUE;
               group_descendant->type&=(~TREE_CASE_MERGER); // We want the last ghost in this chain to be tagged a merger, not it's base halo
            }
         }

         if(flag_create_ghosts){
            ghost_index=n_groups[j_file-i_offset]+n_group_ghosts_used[file];

            // If this is the first ghost for this ghost-chain, correct the information in the
            //    base halo so that it points properly to the ghost(s) we are creating ...
            if(i_offset==1){
               group_descendant->file_offset  =1;
               group_descendant->file_index   =ghost_index;
               group_descendant->descendant_id=group_descendant->id; // The descendant will change from the halo's id only for the last ghost (set later)
            }
            group_descendant=&(groups_in[file%n_wrap][ghost_index]); 
         }
         else{
            ghost_index     =n_groups[j_file-i_offset]+n_group_ghosts[file];
            group_descendant=&group_dummy;
         }

         // Create the new halo (ghost *groups* are trivially added to the end of the list)
         group_descendant->halo_index =ghost_index; 
         group_descendant->id         =groups_in[i_file%n_wrap][i_group].id;
         group_descendant->type       =group_ghost_base_type; 
         group_descendant->tree_id    =groups_in[i_file%n_wrap][i_group].tree_id;
         group_descendant->file_offset=1;

         // Store the information we will need for interpolating the ghost halo properties
         if(flag_create_ghosts){
            group_descendant->interp.file_start =i_file;
            group_descendant->interp.index_start=i_group;
            group_descendant->interp.file_stop  =i_file+group_file_offset;
            group_descendant->interp.index_stop =group_final_file_index;
            if(groups_in[i_file%n_wrap][i_group].file_index<0)
               SID_trap_error("Invalid stopping state for group interpolation: file=%d i_file=%d index=%d type=%d",ERROR_LOGIC,
                              file,i_file,group_descendant->interp.index_stop,group_descendant->type);
            group_descendant->interp.time_start =deltat_a(cosmo,a_list[j_file],0.)/S_PER_YEAR;
            group_descendant->interp.time_stop  =deltat_a(cosmo,a_list[j_file-group_file_offset],0.)/S_PER_YEAR;
         }

         // If this is the last in this chain of ghost groups, force the descendant to the known solution
         if(i_offset==(group_file_offset-1)){ 
            group_descendant->file_index   =group_final_file_index; 
            group_descendant->descendant_id=groups_in[i_file%n_wrap][i_group].descendant_id; 
         }
         // Keep the base ID for the ghosts. If this is a merger, the last ghost's type will be changed later.
         else{
            if(flag_create_ghosts)
               group_descendant->file_index=n_groups[j_file-i_offset-1]+n_group_ghosts_used[file+1]; 
            else
               group_descendant->file_index=n_groups[j_file-i_offset-1]+n_group_ghosts[file+1]; 
            group_descendant->descendant_id=groups_in[i_file%n_wrap][i_group].id; 
         }

         // Initially, ghost groups have no substructure
         group_descendant->n_subgroups       =0;
         group_descendant->first_substructure=NULL;
         group_descendant->last_substructure =NULL;

         if(flag_create_ghosts)
            n_group_ghosts_used[file]++;
         else
            n_group_ghosts[file]++;

         // Compile the descendant history of this group so we know
         //    where to put the substructures after this
         flag_descenadnt_is_new_ghost=TRUE;
         group_descendants[i_offset] =group_descendant;
      } // i_offset

      // Finalize a group ghost-chain
      if(flag_descenadnt_is_new_ghost && flag_create_ghosts){
         // Set the last ghost halo in its chain to be a merger if the 
         //    halo at the start of the chain is labeled as a merger
         if(flag_ghost_merger)
            group_descendant->type|=TREE_CASE_MERGER;
      }

      // Finish off the group's descendant history.  No ghost construction needed here.
      for(;i_offset<n_search && file<=i_file_start;i_offset++,file++){
         if(group_descendant!=NULL){
            if(group_descendant->file_index>=0)
               group_descendant=&(groups_in[file%n_wrap][group_descendant->file_index]); 
            else
               group_descendant=NULL;
         }
         group_descendants[i_offset]=group_descendant;
      }

      // Process the subgroups now
      for(j_subgroup=0;j_subgroup<n_subgroups_group[i_file%n_wrap][i_group];j_subgroup++,i_subgroup++){
         tree_horizontal_ghost_subgroup_info *subgroup_progenitor=NULL;
         tree_horizontal_ghost_subgroup_info *subgroup_descendant;
         int subgroup_file_offset;
         int subgroup_ghost_base_type;
         int subgroup_final_file_index;
         int flag_progenitor_is_new_ghost;
         int flag_ghost_merger;
         int i_offset;
         flag_progenitor_is_new_ghost=FALSE;
         flag_ghost_merger           =FALSE;
         subgroup_progenitor         =&(subgroups_in[i_file%n_wrap][i_subgroup]);
         subgroup_file_offset        =subgroups_in[i_file%n_wrap][i_subgroup].file_offset;
         flag_ghost_merger           =check_mode_for_flag(subgroups_in[i_file%n_wrap][i_subgroup].type,TREE_CASE_MERGER);
         for(i_offset=1,file=i_file+1;i_offset<subgroup_file_offset && file<=i_file_start;i_offset++,file++){
            if(flag_create_ghosts){
               ghost_index=n_subgroups[j_file-i_offset]+n_subgroup_ghosts_used[i_file+i_offset];

               // If this is the first ghost for this group, do some things ...
               if(i_offset==1){
                  // ... set the descendant pointer of the base halo.  We weren't able to set this before when we
                  //     read the halo in and called add_substructure_to_horizontal_tree_group(). 
                  //     It was set to NULL at that time.
                  if(subgroup_progenitor->descendant!=NULL)
                     SID_trap_error("A substructure (i_group=%d,i_subgroup=%d) with file_offset>1 did not have it's descendant pointer initialized to NULL properly.",ERROR_LOGIC,i_group,i_subgroup);
                  subgroup_progenitor->descendant=&(subgroups_in[(i_file+1)%n_wrap][ghost_index]);

                  // ... change the descendant info of the base halo
                  subgroup_final_file_index         = subgroup_progenitor->file_index;
                  subgroup_progenitor->file_offset  = 1;
                  subgroup_progenitor->halo_index   =-1;                      // This will be set later when we perform the write
                  subgroup_progenitor->file_index   =-1;                      // This will be set later when we perform the write
                  subgroup_progenitor->descendant_id=subgroup_progenitor->id; // The descendant will change from the halo's id only for the last ghost 

                  // ... set the type for the ghost halos.  A merger flag may be added later for the last one.
                  subgroup_ghost_base_type=TREE_CASE_GHOST;
                  if(check_mode_for_flag(subgroup_progenitor->type,TREE_CASE_FRAGMENTED_LOST))
                     subgroup_ghost_base_type|=TREE_CASE_FRAGMENTED_LOST;
                  if(check_mode_for_flag(subgroup_progenitor->type,TREE_CASE_FRAGMENTED_RETURNED))
                     subgroup_ghost_base_type|=TREE_CASE_FRAGMENTED_RETURNED;
                  if(check_mode_for_flag(subgroup_progenitor->type,TREE_CASE_FRAGMENTED_EXCHANGED))
                     subgroup_ghost_base_type|=TREE_CASE_FRAGMENTED_EXCHANGED;
                  if(check_mode_for_flag(subgroup_progenitor->type,TREE_CASE_SPUTTERED))
                     subgroup_ghost_base_type|=TREE_CASE_SPUTTERED;
                  if(check_mode_for_flag(subgroup_progenitor->type,TREE_CASE_MAIN_PROGENITOR))
                     subgroup_ghost_base_type|=TREE_CASE_MAIN_PROGENITOR;
                  if(check_mode_for_flag(subgroup_progenitor->type,TREE_CASE_MERGER)){
                     flag_ghost_merger=TRUE;
                     subgroup_progenitor->type&=(~TREE_CASE_MERGER); // We want the last ghost in this chain to be tagged a merger, not it's base halo
                     subgroup_progenitor->type|=TREE_CASE_SIMPLE;
                  }
               }

               // Store the information we will need for interpolating the ghost halo properties
               subgroups_in[file%n_wrap][ghost_index].interp.file_start =i_file;
               subgroups_in[file%n_wrap][ghost_index].interp.index_start=i_subgroup;
               subgroups_in[file%n_wrap][ghost_index].interp.file_stop  =i_file+subgroup_file_offset;
               subgroups_in[file%n_wrap][ghost_index].interp.index_stop =subgroup_final_file_index;
               if(subgroups_in[file%n_wrap][ghost_index].interp.index_stop<0)
                  SID_trap_error("Invalid stopping state for subgroup interpolation: file=%d i_file=%d index=%d type=%d offset=%d",ERROR_LOGIC,
                                 file,i_file,subgroups_in[file%n_wrap][ghost_index].interp.index_stop,
                                 subgroups_in[file%n_wrap][ghost_index].type,subgroups_in[i_file%n_wrap][i_subgroup].file_offset);
               subgroups_in[file%n_wrap][ghost_index].interp.time_start =deltat_a(cosmo,a_list[j_file],0.)/S_PER_YEAR;
               subgroups_in[file%n_wrap][ghost_index].interp.time_stop  =deltat_a(cosmo,a_list[j_file-group_file_offset],0.)/S_PER_YEAR;

               // Create the new halo
               subgroups_in[file%n_wrap][ghost_index].halo_index =-1; // This will be set later when we perform the write
               subgroups_in[file%n_wrap][ghost_index].id         =subgroups_in[i_file%n_wrap][i_subgroup].id;
               subgroups_in[file%n_wrap][ghost_index].type       =subgroup_ghost_base_type;
               subgroups_in[file%n_wrap][ghost_index].tree_id    =subgroups_in[i_file%n_wrap][i_subgroup].tree_id;
               subgroups_in[file%n_wrap][ghost_index].file_offset=1;

               // If this is the last ghost for this group's ghost chain, force the descendant to the known solution
               if(i_offset==(subgroup_file_offset-1)){
                  subgroups_in[file%n_wrap][ghost_index].file_index   =subgroup_final_file_index;
                  subgroups_in[file%n_wrap][ghost_index].descendant_id=subgroups_in[i_file%n_wrap][i_subgroup].descendant_id;
               }
               // Keep the base ID for the ghosts. If this is a merger, the last ghost will be changed later to the descendant_id.
               else{
                  subgroups_in[file%n_wrap][ghost_index].file_index   =n_subgroups[j_file-i_offset-1]+n_subgroup_ghosts_used[file+1];
                  subgroups_in[file%n_wrap][ghost_index].descendant_id=subgroups_in[i_file%n_wrap][i_subgroup].id;
               }
               subgroup_descendant=&(subgroups_in[(file+1)%n_wrap][subgroups_in[file%n_wrap][ghost_index].file_index]);

               // Add the ghost subgroup to it's group.  Create null group ghosts if needed.
               // Create a new null ghost group that just holds the ghost subgroup...
               if(group_descendants[i_offset]==NULL){
                  int null_ghost_index;
                  null_ghost_index=n_groups[j_file-i_offset]+n_group_ghosts_used[file];
                  groups_in[file%n_wrap][null_ghost_index].halo_index        = null_ghost_index; 
                  groups_in[file%n_wrap][null_ghost_index].id                =-1;
                  groups_in[file%n_wrap][null_ghost_index].type              = TREE_CASE_GHOST|TREE_CASE_GHOST_NULL;
                  groups_in[file%n_wrap][null_ghost_index].tree_id           =-1;
                  groups_in[file%n_wrap][null_ghost_index].file_offset       =-1;
                  groups_in[file%n_wrap][null_ghost_index].file_index        =-1; 
                  groups_in[file%n_wrap][null_ghost_index].descendant_id     =-1; 
                  groups_in[file%n_wrap][null_ghost_index].n_subgroups       = 0;
                  groups_in[file%n_wrap][null_ghost_index].first_substructure= NULL;
                  groups_in[file%n_wrap][null_ghost_index].last_substructure = NULL;
                  add_substructure_to_horizontal_tree_group(&(groups_in[file%n_wrap][null_ghost_index]),
                                                            subgroup_descendant,
                                                            &(subgroups_in[file%n_wrap][ghost_index]));
                  n_group_ghosts_used[file]++;
               }
               // ... or add the substructure to an existing group.
               else{
                  add_substructure_to_horizontal_tree_group(group_descendants[i_offset],
                                                            subgroup_descendant,
                                                            &(subgroups_in[file%n_wrap][ghost_index]));
               }
               subgroup_progenitor         =&(subgroups_in[file%n_wrap][ghost_index]);
               flag_progenitor_is_new_ghost=TRUE;
               n_subgroup_ghosts_used[file]++;
            }
            else{
               if(group_descendants[i_offset]==NULL){
                  n_group_ghosts[file]++;
               }
               n_subgroup_ghosts[file]++;
            }
         } // i_offset
         if(flag_progenitor_is_new_ghost && flag_ghost_merger && flag_create_ghosts)
            subgroup_progenitor->type|=TREE_CASE_MERGER;
      } // j_subgroup
   }
   SID_free(SID_FARG group_descendants);
//int i_test;
//fprintf(stderr,"\n");
//if(flag_create_ghosts) for(i_test=0;i_test<n_files;i_test++) fprintf(stderr,"%3d %6d %6d\n",i_test,n_group_ghosts_used[i_test],n_subgroup_ghosts_used[i_test]);
//else                   for(i_test=0;i_test<n_files;i_test++) fprintf(stderr,"%3d %6d %6d\n",i_test,n_group_ghosts[i_test],     n_subgroup_ghosts[i_test]);
}

