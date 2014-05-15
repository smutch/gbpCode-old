#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

// Decide what to do when halo_i[i_halo] is matched to halo_j[j_halo]
void set_halo_and_descendant(tree_horizontal_info **halos,
                             int                    i_file,
                             int                    i_halo,
                             int                    j_file,
                             int                    j_halo,
                             float                  score,
                             int                   *max_id,
                             int                    n_wrap){
   tree_horizontal_info *halos_i;
   tree_horizontal_info *halos_j;
   int                   file_offset;
   int                   k_file;
   int                   l_file;
   int                   k_file_temp;
   int                   k_file_main;
   int                   k_size_main;
   int                   flag_process;
   int                   flag_emerged;
   int                   n_p_diff_old;
   int                   n_p_diff_new;

   if(j_file<=i_file)
     SID_trap_error("j_file<=i_file in set_halo_and_descendant().",ERROR_NONE);

   // Process the inputs a bit
   halos_i    =halos[i_file%n_wrap];
   halos_j    =halos[j_file%n_wrap];
   file_offset=j_file-i_file;
   if(file_offset==0)
      SID_trap_error("A zero file offset has been requested.  It should be -ve for roots and +ve otherwise.",ERROR_LOGIC);

   // If the score exceeds the allowed match score, this must be a forced match.  Make it so.
   int flag_forced=(score>MAX_TREE_MATCH_SCORE);

if(halos_i[i_halo].id==0 || halos_j[j_halo].id==0){
char *halo_type_string=NULL;
tree_case_flags_text(halos_i[i_halo].type,"+",&halo_type_string);
fprintf(stderr,"A:%05d/%05d %05d/%05d %le %s\n",halos_i[i_halo].snap,halos_i[i_halo].index,halos_j[j_halo].snap,halos_j[j_halo].index,score,halo_type_string);
SID_free(SID_FARG halo_type_string);
}

   // Set non-bridged halos or finalize bridge matches (ie. apply defaults for bridge progenitors not matched to emerged halos)
   if(!check_mode_for_flag(halos_j[j_halo].type,TREE_CASE_BRIDGED)                       ||
       check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED) ||
       check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_FINALIZE)               ||
       flag_forced){

      // If we are processing a halo already identified as matched to a bridge and we are not
      //    finalizing or forcing this match, only accept this new match if it meets these criteria ...
      flag_process=TRUE;
      flag_emerged=FALSE;
      if(!check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_FINALIZE)               &&
         !flag_forced                                                                       &&
          check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED)){

        // ... we are not matching to a descendant of the initial bridged match ...
        tree_horizontal_info *current;
        current=halos_i[i_halo].bridge_forematch.halo->descendant.halo;
        if(current!=NULL)
           k_file =current->file;
        l_file=k_file;
        while(current!=NULL &&
              k_file>=l_file && k_file<=j_file && // not true when we reach past the rolling array bounds
              flag_process){
           if(current==&(halos_j[j_halo]))
             flag_process=FALSE;
if(halos_i[i_halo].id==0 || halos_j[j_halo].id==0 && !flag_process) fprintf(stderr,"   Denied - Descendant of bridge\n");
           current=current->descendant.halo;
           l_file=k_file;
           if(current!=NULL)
              k_file =current->file;
        }
        flag_emerged=flag_process;
      }

      // If the code above hasn't decided that we shouldn't process this match ...
      if(flag_process){

         // Set progenitor IDs and pointers ...
         match_info old_progenitor;
         match_info new_progenitor;
         // ... create new progenitor ...
         new_progenitor.halo =&(halos_i[i_halo]);
         new_progenitor.score=score;
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
            new_progenitor.halo->type|=  TREE_CASE_MAIN_PROGENITOR;
            new_progenitor.halo->type&=(~TREE_CASE_MERGER);
         }
         // ... else add a next progenitor ...
         else{
            // If we have a better main-progenitor, insert it at the 
            //   beginning of the list and swap IDs with the main progenitor so that
            //   the correct halo gets the main progenitor ID and all others get a new one ...
            memcpy(&old_progenitor,&(halos_j[j_halo].first_progenitor),sizeof(match_info));
            // Progenitor order needs to be set in a way that ensures correct emerged and fragmented halo processing
            //if(score>old_progenitor.score){ 
            if((new_progenitor.halo->n_particles)>(old_progenitor.halo->n_particles)){
               // ... let the new main progenitor inherit the descendant's id ...
               halos_i[i_halo].id=halos_j[j_halo].id;
               // ... and give the old main progenitor the new id ...
               change_horizontal_ID_recursive(old_progenitor.halo,old_progenitor.halo->id,(*max_id)++);
               // ... set new main progenitor type ...
               old_progenitor.halo->type|=  TREE_CASE_MERGER;
               old_progenitor.halo->type&=(~TREE_CASE_MAIN_PROGENITOR);
               old_progenitor.halo->type&=(~TREE_CASE_SIMPLE);
               new_progenitor.halo->type&=(~TREE_CASE_MERGER);
               new_progenitor.halo->type|=  TREE_CASE_MAIN_PROGENITOR;
               new_progenitor.halo->type&=(~TREE_CASE_SIMPLE);
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
               halos_i[i_halo].type|=  TREE_CASE_MERGER;
               halos_i[i_halo].type&=(~TREE_CASE_MAIN_PROGENITOR);
               halos_i[i_halo].type&=(~TREE_CASE_SIMPLE);
            }
         }

         // ... set last progenitor
         memcpy(&(halos_j[j_halo].last_progenitor),&new_progenitor,sizeof(match_info));

         // Set descendant info
         halos_i[i_halo].descendant.halo =&(halos_j[j_halo]);
         halos_i[i_halo].descendant.score=score;

         // Set match-type flags ...                   
         //   Propagate strayed halo flags ...
         if(check_mode_for_flag(halos_j[j_halo].type,TREE_CASE_STRAYED))
            halos_i[i_halo].type|=TREE_CASE_STRAYED;

         // ... set the flags for simple and dropped halo matches
         if(file_offset==1)
            halos_i[i_halo].type|=TREE_CASE_SIMPLE;
         else if(!check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_MATCHED_TO_BRIDGE))
            halos_i[i_halo].type|=TREE_CASE_DROPPED;

         // ... set the flag for emerged halos
         if(flag_emerged)
            halos_j[j_halo].type|=TREE_CASE_EMERGED;

         // ... turn off the TREE_CASE_NO_PROGENITORS flag for the descendant ...
         halos_j[j_halo].type&=(~TREE_CASE_NO_PROGENITORS);

         // Mark the halo as processed
         halos_i[i_halo].type&=(~(TREE_CASE_UNPROCESSED|TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED|TREE_CASE_BRIDGE_FINALIZE));
      }
   }
   // ... else we've matched to a bridge.  Here we just set the info needed
   //     to connect this halo with the bridge it's been matched to.  We'll
   //     use that info to search the bridge's list of possible emergent halos and if we fail
   //     to identify to one of those, we'll finalize this match as the default.  The check
   //     is constructed in such a way that only the first match to a bridge (and not say, to
   //     it's descendants) is used.
   else if(halos_i[i_halo].bridge_forematch.halo==NULL){
      halos_i[i_halo].bridge_forematch.halo =&(halos_j[j_halo]);
      halos_i[i_halo].bridge_forematch.score =score;
      halos_i[i_halo].type                  |=(TREE_CASE_MATCHED_TO_BRIDGE|TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED);
   }
if(halos_i[i_halo].id==0 || halos_j[j_halo].id==0){
char *halo_type_string=NULL;
tree_case_flags_text(halos_i[i_halo].type,"+",&halo_type_string);
fprintf(stderr,"B:%05d/%05d %05d/%05d %le %s\n",halos_i[i_halo].snap,halos_i[i_halo].index,halos_j[j_halo].snap,halos_j[j_halo].index,score,halo_type_string);
SID_free(SID_FARG halo_type_string);
}
}

