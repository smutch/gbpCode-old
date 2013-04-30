#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees.h>

void check_for_fragmented_halos(int k_match,tree_horizontal_info **groups,int n_groups,
                               int i_write,int j_write,int n_wrap){
   if(k_match==0)
      SID_log("Checking for fragmented groups in snapshot #%03d...",SID_LOG_OPEN,j_write);
   else
      SID_log("Checking for fragmented subgroups in snapshot #%03d...",SID_LOG_OPEN,j_write);

   // Finalize halos
   int i_group;
   int n_lost     =0;
   int n_returned =0;
   int n_exchanged=0;
   for(i_group=0;i_group<n_groups;i_group++){
      // Perform some sanity checks on the match_type flag
      if(check_mode_for_flag(groups[i_write%n_wrap][i_group].type,TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED) ||
         check_mode_for_flag(groups[i_write%n_wrap][i_group].type,TREE_CASE_BRIDGE_FINALIZE)               ||
         check_mode_for_flag(groups[i_write%n_wrap][i_group].type,TREE_CASE_UNPROCESSED)                   ||
         check_mode_for_flag(groups[i_write%n_wrap][i_group].type,TREE_CASE_INVALID))
         SID_trap_error("Invalid group match_type flag (%d) for i_group=%d",ERROR_LOGIC,groups[i_write%n_wrap][i_group].type,i_group);

      // Perform checks for fragmented halos here
      if(check_mode_for_flag(groups[i_write%n_wrap][i_group].type,TREE_CASE_FRAGMENTED_NEW)){
         int bridge_id;
         bridge_id=set_match_id(&(groups[i_write%n_wrap][i_group].bridge_backmatch));
         if(groups[i_write%n_wrap][i_group].main_progenitor_id<0){
            groups[i_write%n_wrap][i_group].type|=TREE_CASE_FRAGMENTED_LOST;
            n_lost++;
         }
         else if(groups[i_write%n_wrap][i_group].main_progenitor_id!=bridge_id){
            groups[i_write%n_wrap][i_group].type|=TREE_CASE_FRAGMENTED_EXCHANGED;
            n_exchanged++;
         }
         else if(groups[i_write%n_wrap][i_group].main_progenitor_id==bridge_id){
            groups[i_write%n_wrap][i_group].type|=TREE_CASE_FRAGMENTED_RETURNED;
            n_returned++;
         }
      }
   }
   SID_log("# of new lost      fragmented halos = %-8d",SID_LOG_COMMENT,n_lost);
   SID_log("# of new returned  fragmented halos = %-8d",SID_LOG_COMMENT,n_returned);
   SID_log("# of new exchanged fragmented halos = %-8d",SID_LOG_COMMENT,n_exchanged);
   SID_log("Done.",SID_LOG_CLOSE);
}

