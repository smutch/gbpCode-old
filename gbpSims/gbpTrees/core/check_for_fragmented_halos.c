#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void check_for_fragmented_halos(int k_match,
                                tree_horizontal_info **halos,
                                int n_halos,
                                int i_write,
                                int j_write,
                                int l_write,
                                int n_wrap){
   if(k_match==0)
      SID_log("Checking for fragmented subgroups in snapshot #%03d...",SID_LOG_OPEN,j_write);
   else
      SID_log("Checking for fragmented groups in snapshot #%03d...",SID_LOG_OPEN,j_write);

   // Finalize halos
   int i_halo;
   int n_strayed  =0;
   int n_returned =0;
   int n_exchanged=0;
   for(i_halo=0;i_halo<n_halos;i_halo++){
      // Perform some sanity checks on the match_type flag
      if(check_mode_for_flag(halos[i_write%n_wrap][i_halo].type,TREE_CASE_UNPROCESSED) ||
         check_mode_for_flag(halos[i_write%n_wrap][i_halo].type,TREE_CASE_INVALID))
         SID_trap_error("Invalid halo match_type flag (%d) for i_halo=%d",ERROR_LOGIC,halos[i_write%n_wrap][i_halo].type,i_halo);

      // Perform checks for fragmented halos here
      match_info *back_match =&(halos[i_write%n_wrap][i_halo].bridge_backmatch);
      if(check_mode_for_flag(halos[i_write%n_wrap][i_halo].type,TREE_CASE_NO_PROGENITORS) &&
         (back_match->halo)!=NULL &&
         l_write!=0){ 

         // Decide if this is a fragmented halo or a fragmented halo source
         int                   n_bridge_back_matches=back_match->halo->n_back_matches;
         tree_horizontal_info *frag_source=back_match->halo->back_matches[0].halo; // Sorted by largest descendant size (DESCENDAING ORDER)
         //tree_horizontal_info *frag_source=back_match->halo->descendant.halo; 
         int flag_fragment_source=((&(halos[i_write%n_wrap][i_halo]))==frag_source);

         if(!flag_fragment_source){ 
            // We've identified a fragmented halo.  First, set/unset a couple flags ...
            halos[i_write%n_wrap][i_halo].type|=TREE_CASE_FRAGMENTED_NEW;

            // ... then decide what type of fragmented halo it is.
            int bridge_id;
            int halo_id;
            int descendant_id;
            int main_progenitor_id;
            bridge_id         =set_match_id(back_match);
            halo_id           =halos[i_write%n_wrap][i_halo].id;
            descendant_id     =set_match_id(&(halos[i_write%n_wrap][i_halo].descendant));
            main_progenitor_id=halos[i_write%n_wrap][i_halo].main_progenitor_id;
            if(halo_id<0 || check_mode_for_flag(halos[i_write%n_wrap][i_halo].type,TREE_CASE_STRAYED)){
               halos[i_write%n_wrap][i_halo].type|=TREE_CASE_FRAGMENTED_STRAYED;
               n_strayed++;
            }
            else if(bridge_id==main_progenitor_id){
               halos[i_write%n_wrap][i_halo].type|=TREE_CASE_FRAGMENTED_RETURNED;
               n_returned++;
            }
            else{
               halos[i_write%n_wrap][i_halo].type|=TREE_CASE_FRAGMENTED_EXCHANGED;
               n_exchanged++;
            }
         }
      }
   }
   if(n_strayed!=0 || n_returned!=0 || n_exchanged!=0){
      SID_log("# of new strayed   fragmented halos = %-8d",SID_LOG_COMMENT,n_strayed);
      SID_log("# of new returned  fragmented halos = %-8d",SID_LOG_COMMENT,n_returned);
      SID_log("# of new exchanged fragmented halos = %-8d",SID_LOG_COMMENT,n_exchanged);
   }
   else
      SID_log("none found...",SID_LOG_CONTINUE);
   SID_log("Done.",SID_LOG_CLOSE);
}

