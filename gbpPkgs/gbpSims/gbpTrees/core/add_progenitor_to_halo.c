#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void add_progenitor_to_halo(tree_horizontal_info **halos,
                            int                    i_file,
                            int                    i_halo,
                            int                    j_file,
                            int                    j_halo,
                            float                  score,
                            int                    flag_two_way,
                            int                    flag_back_match,
                            int                   *max_id,
                            int                    n_wrap){
   tree_horizontal_info *halos_i=halos[i_file%n_wrap];
   tree_horizontal_info *halos_j=halos[j_file%n_wrap];
   int                   file_offset=j_file-i_file;

   // Set progenitor IDs and pointers ...
   match_info old_progenitor;
   match_info new_progenitor;
   // ... create new progenitor ...
   new_progenitor.halo           =&(halos_i[i_halo]);
   new_progenitor.score          =score;
   new_progenitor.flag_two_way   =flag_two_way;
   new_progenitor.flag_back_match=flag_back_match;
   // ... increment counter ...
   halos_j[j_halo].n_progenitors++;
   // ... set main progenitor id ...
   new_progenitor.halo->main_progenitor_id=halos_j[j_halo].id;
   // ... set tree id ...
   new_progenitor.halo->tree_id=halos_j[j_halo].tree_id;
   // ... create initial progenitor ...
   if(halos_j[j_halo].n_progenitors==1){
      // ... set initial progenitor ...
      memcpy(&(halos_j[j_halo].first_progenitor),&new_progenitor,sizeof(match_info));
      // ... set initial progenitor id ...
      new_progenitor.halo->id=halos_j[j_halo].id;
      // ... set initial progenitor type ...
      new_progenitor.halo->type|=TREE_CASE_MAIN_PROGENITOR;
   }
   // ... else add a next progenitor ...
   else{
      // Set remnant flag
      halos_j[j_halo].type|=TREE_CASE_REMNANT;
      // If we have a better main-progenitor, insert it at the 
      //   beginning of the list and swap IDs with the main progenitor so that
      //   the correct halo gets the main progenitor ID and all others get a 
      //   new one.  Progenitor order can be decided many ways.  Here we choose
      //   to rank by immediacy and n_particles.  We choose this to minimize the cases of
      //   a small halo becoming the main progenitor due only to a good centrally
      //   weighted match.  One advanage of this is that we minimize the cases
      //   of fragmented halos being identified as 'exchanged' or of being
      //   propagated past the point where they merge.
      memcpy(&old_progenitor,&(halos_j[j_halo].first_progenitor),sizeof(match_info));
      int file_offset_new=(halos_j[j_halo].file)-(new_progenitor.halo->file);
      int file_offset_old=(halos_j[j_halo].file)-(old_progenitor.halo->file);
      if(file_offset_new<=0 || file_offset_old<=0)
         SID_trap_error("Progenitor file offset is invalid (ie either %d<=0 or %d<=0).",ERROR_LOGIC,file_offset_new,file_offset_old);
      if(check_if_match_is_bigger(&(halos_j[j_halo]),&old_progenitor,&new_progenitor)){
      //if(file_offset_new<file_offset_old || (file_offset_new==file_offset_old && new_progenitor.halo->n_particles>old_progenitor.halo->n_particles)){
         // ... let the new main progenitor inherit the descendant's id ...
         halos_i[i_halo].id=halos_j[j_halo].id;
         // ... and give the old main progenitor the new id ...
         change_horizontal_ID_recursive(old_progenitor.halo,old_progenitor.halo->id,(*max_id)++);
         // ... set new main progenitor type ...
         old_progenitor.halo->type&=(~TREE_CASE_MAIN_PROGENITOR);
         new_progenitor.halo->type|=  TREE_CASE_MAIN_PROGENITOR;
         // ... set new main progenitor ...
         memcpy(&(halos_j[j_halo].first_progenitor),                      &new_progenitor,sizeof(match_info));
         memcpy(&(halos_j[j_halo].first_progenitor.halo->next_progenitor),&old_progenitor,sizeof(match_info));
      }
      // ... else just add the new halo to the end of the list and create a new ID for it (if this is not a strayed halo).
      else{
         // ... set new non-main progenitor ...
         memcpy(&(halos_j[j_halo].last_progenitor.halo->next_progenitor),&new_progenitor,sizeof(match_info));
         // ... set new non-main progenitor id ...
         if(halos_j[j_halo].id>=0)
            halos_i[i_halo].id=(*max_id)++;
         else
            halos_i[i_halo].id=halos_j[j_halo].id;
         // ... set new non-main progenitor type ...
         halos_i[i_halo].type&=(~TREE_CASE_MAIN_PROGENITOR);
      }
   }

   // ... set last progenitor
   memcpy(&(halos_j[j_halo].last_progenitor),&new_progenitor,sizeof(match_info));

   // Set descendant info
   halos_i[i_halo].descendant.halo =&(halos_j[j_halo]);
   halos_i[i_halo].descendant.score=score;

   // Set match-type flags ...                   
   //   Propagate strayed halo flags (needed for identifying 
   //   strayed fragmented halos later, for example)...
   if(check_mode_for_flag(halos_j[j_halo].type,TREE_CASE_STRAYED))
      halos_i[i_halo].type|=TREE_CASE_STRAYED;

   // ... set the flags for emerged halo matches
   if(check_mode_for_flag(halos_j[j_halo].type,TREE_CASE_EMERGED_CANDIDATE)){
      halos_i[i_halo].type|=TREE_CASE_MATCHED_TO_EMERGED;
      halos_j[j_halo].type|=TREE_CASE_EMERGED;
   }

   // ... set the flags for dropped halo matches
   else if(file_offset!=1)
      halos_i[i_halo].type|=TREE_CASE_DROPPED;

   // ... turn off the TREE_CASE_NO_PROGENITORS flag for the descendant ...
   halos_j[j_halo].type&=(~TREE_CASE_NO_PROGENITORS);

   // Mark the halo as processed
   halos_i[i_halo].type&=(~(TREE_CASE_UNPROCESSED));
}
