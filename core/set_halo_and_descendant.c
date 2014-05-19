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
                             int                    n_search,
                             int                    n_wrap,
                             back_match_info       *back_match){
   int                   k_file;
   int                   l_file;
   int                   k_file_temp;
   int                   k_file_main;
   int                   k_size_main;

   // Sanity check
   if(j_file<=i_file)
     SID_trap_error("j_file<=i_file in set_halo_and_descendant().",ERROR_NONE);

   // If we've been passed a back match, it means we're checking emerged halo candidates
   int flag_check_emerged=(back_match!=NULL);

   // Process the inputs a bit
   tree_horizontal_info *halos_i=halos[i_file%n_wrap];
   tree_horizontal_info *halos_j=halos[j_file%n_wrap];;
   int                   file_offset=j_file-i_file;
   if(file_offset==0)
      SID_trap_error("A zero file offset has been requested.  It should be -ve for roots and +ve otherwise.",ERROR_LOGIC);
//int i_read_report=730;
//int i_report     =249;
//if(halos_i[0].snap==i_read_report){
//if(i_halo==i_report){
//char *halo_type_string=NULL;
//tree_case_flags_text(halos_i[i_report].type,"+",&halo_type_string);
//fprintf(stderr,"I-A[%d][%d]: %d %s\n",i_read_report,i_report,halos_i[i_report].n_bridges,halo_type_string);
//SID_free(SID_FARG halo_type_string);
//}
//}
//if(halos_j[0].snap==i_read_report){
//if(j_halo==i_report){
//char *halo_type_string=NULL;
//tree_case_flags_text(halos_j[i_report].type,"+",&halo_type_string);
//fprintf(stderr,"J-A[%d][%d]: %d %s\n",i_read_report,i_report,halos_j[i_report].n_bridges,halo_type_string);
//SID_free(SID_FARG halo_type_string);
//}
//}

   // If we are matching to a bridge ...
   if(check_mode_for_flag(halos_j[j_halo].type,TREE_CASE_BRIDGED)){
      // If we are matching to a bridge for the first time, just set the info needed
      //    to connect this halo with the bridge it's been matched to.  We'll
      //    use that info to search the bridge's list of possible emergent halos and if we fail
      //    to identify to one of those, we'll finalize this match as the default.  The check
      //    is constructed in such a way that only the first match to a bridge (and not say, to
      //    it's descendants) is used.
      if(halos_i[i_halo].bridge_forematch.halo==NULL || flag_check_emerged){
         halos_i[i_halo].bridge_forematch.halo =&(halos_j[j_halo]);
         halos_i[i_halo].bridge_forematch.score =score;
         halos_i[i_halo].type                  |=(TREE_CASE_MATCHED_TO_BRIDGE|TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED);
      }
      // If we are forcing a match to a halo that was matched to a bridge ...
      else if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_FINALIZE)){
         if(halos_i[i_halo].bridge_forematch.halo==NULL)
            SID_trap_error("We have been asked to finalize a bridge match that has not been set.",ERROR_LOGIC);
         add_progenitor_to_halo(halos,
                                i_file,
                                i_halo,
                                j_file,
                                j_halo,
                                score,
                                max_id,
                                n_wrap,
                                FALSE);
      }
   }
   // ... else if we Set non-bridged halos or finalize bridge matches (ie. apply defaults for bridge progenitors not matched to emerged halos)
   else if(flag_check_emerged){

      // Perform some checks to see if we want to make a match to this emerged halo
      int flag_process=TRUE;

      // 1) Because we may be recursively finding emerged halos, make sure we haven't
      //    exceeded the search interval.  If so, ignore this attempted match.
      int total_offset=back_match->file-halos_i[i_halo].file;
      if(total_offset>n_search)
         flag_process=FALSE;

      // 2) We are not matching to a descendant of the initial bridged match ...
      tree_horizontal_info *current;
      current=halos_i[i_halo].bridge_forematch.halo->descendant.halo;
      if(current!=NULL)
         k_file =current->file;
      l_file=k_file;
      // ... loop over the main progenitor line of the original bridge match
      while(current!=NULL &&
            k_file>=l_file && k_file<=j_file && // not true when we reach past the rolling array bounds
            flag_process){
         // ... if we come to the halo we have been asked to match to, cancel the match.
         if(current==&(halos_j[j_halo]))
           flag_process=FALSE;
         // ... if not, keep moving ...
         current=current->descendant.halo;
         l_file=k_file;
         if(current!=NULL)
            k_file =current->file;
      }

      // If the code above hasn't decided that we shouldn't process this match ...
      if(flag_process)
         add_progenitor_to_halo(halos,
                                i_file,
                                i_halo,
                                j_file,
                                j_halo,
                                score,
                                max_id,
                                n_wrap,
                                TRUE);
   }
   // ... else we are just forming a simple match.
   else{
      add_progenitor_to_halo(halos,
                             i_file,
                             i_halo,
                             j_file,
                             j_halo,
                             score,
                             max_id,
                             n_wrap,
                             FALSE);
   }

//if(halos_i[0].snap==i_read_report){
//if(i_halo==i_report){
//char *halo_type_string=NULL;
//tree_case_flags_text(halos_i[i_report].type,"+",&halo_type_string);
//fprintf(stderr,"I-B[%d][%d]: %d %s\n",i_read_report,i_report,halos_i[i_report].n_bridges,halo_type_string);
//SID_free(SID_FARG halo_type_string);
//}
//}
//if(halos_j[0].snap==i_read_report){
//if(j_halo==i_report){
//char *halo_type_string=NULL;
//tree_case_flags_text(halos_j[i_report].type,"+",&halo_type_string);
//fprintf(stderr,"J-B[%d][%d]: %d %s\n",i_read_report,i_report,halos_j[i_report].n_bridges,halo_type_string);
//SID_free(SID_FARG halo_type_string);
//}
//}
}

