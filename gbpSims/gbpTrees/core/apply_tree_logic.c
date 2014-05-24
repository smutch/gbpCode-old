#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void apply_tree_logic(tree_horizontal_info **halos,
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
   tree_horizontal_info *halos_j=halos[j_file%n_wrap];
   int                   file_offset=j_file-i_file;
   if(file_offset==0)
      SID_trap_error("A zero file offset has been requested.  It should be -ve for roots and +ve otherwise.",ERROR_LOGIC);
//int i_read_report=842;
//int i_report     =72;
int i_read_report=562;
int i_report     =20;
int flag_report_i=FALSE;
int flag_report_j=FALSE;
if(halos_i[0].snap==i_read_report && i_halo==i_report) flag_report_i=TRUE;
if(halos_j[0].snap==i_read_report && j_halo==i_report) flag_report_j=TRUE;
if(flag_report_i){
char *halo_type_string_1=NULL;
char *halo_type_string_2=NULL;
tree_case_flags_text(halos_i[i_halo].type,"+",&halo_type_string_1);
tree_case_flags_text(halos_j[j_halo].type,"+",&halo_type_string_2);
fprintf(stderr,"\nREPORT: I-A[%d][%d]->[%d][%d]: %d %d {%s} {%s} -- %d\n",halos_i[i_halo].snap,i_halo,halos_j[j_halo].snap,j_halo,halos_i[i_halo].n_back_matches,halos_j[j_halo].n_back_matches,halo_type_string_1,halo_type_string_2,flag_check_emerged);
SID_free(SID_FARG halo_type_string_1);
SID_free(SID_FARG halo_type_string_2);
}
if(flag_report_j){
char *halo_type_string_1=NULL;
char *halo_type_string_2=NULL;
tree_case_flags_text(halos_i[i_halo].type,"+",&halo_type_string_1);
tree_case_flags_text(halos_j[j_halo].type,"+",&halo_type_string_2);
fprintf(stderr,"\nREPORT: J-A[%d][%d]->[%d][%d]: %d %d {%s} {%s} -- %d\n",halos_i[i_halo].snap,i_halo,halos_j[j_halo].snap,j_halo,halos_i[i_halo].n_back_matches,halos_j[j_halo].n_back_matches,halo_type_string_1,halo_type_string_2,flag_check_emerged);
SID_free(SID_FARG halo_type_string_1);
SID_free(SID_FARG halo_type_string_2);
}

   // If we are finalizing a match to a bridge ...
   if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_BRIDGE_FINALIZE)){
      if(halos_i[i_halo].bridge_forematch_first.halo==NULL)
         SID_trap_error("We have been asked to finalize a bridge match that has not been set.",ERROR_LOGIC);
if(flag_report_i || flag_report_j) fprintf(stderr,"REPORT:  -> testA\n");
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
   // If we are matching to a bridge (either for the first time, or while checking an emerged
   //    halo candidate) then just set the info needed to connect this halo with the bridge it's 
   //    been matched to.  We'll use that info to search the bridge's list of possible emergent 
   //    halos and if we fail to identify to one of those, we'll finalize this match as the 
   //    default.  The check is constructed in such a way that only the first match to a bridge 
   //    (and not say, to it's descendants) is used.
   else if(check_mode_for_flag(halos_j[j_halo].type,TREE_CASE_BRIDGED)){
      if(halos_i[i_halo].bridge_forematch_first.halo==NULL){
if(flag_report_i || flag_report_j) fprintf(stderr,"REPORT:  -> testB1\n");
         halos_i[i_halo].bridge_forematch_first.halo   =&(halos_j[j_halo]);
         halos_i[i_halo].bridge_forematch_first.score  =score;
         halos_i[i_halo].bridge_forematch_default.halo =&(halos_j[j_halo]);
         halos_i[i_halo].bridge_forematch_default.score=score;
         halos_i[i_halo].bridge_forematch.halo         =&(halos_j[j_halo]);
         halos_i[i_halo].bridge_forematch.score        =score;
         halos_i[i_halo].type|=(TREE_CASE_MATCHED_TO_BRIDGE|TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED);
      }
      else if(!check_if_halo_is_descendant(halos_i[i_halo].bridge_forematch_default.halo,&(halos_j[j_halo]),n_search)){
if(flag_report_i || flag_report_j) fprintf(stderr,"REPORT:  -> testB2\n");
         halos_i[i_halo].bridge_forematch_default.halo =&(halos_j[j_halo]);
         halos_i[i_halo].bridge_forematch_default.score=score;
      }
      halos_i[i_halo].bridge_forematch.halo =&(halos_j[j_halo]);
      halos_i[i_halo].bridge_forematch.score=score;
if(flag_report_i || flag_report_j) fprintf(stderr,"REPORT:  -> testB3\n");
   }
   else if(flag_check_emerged){
if(flag_report_i || flag_report_j) fprintf(stderr,"REPORT:  -> testB4\n");
      if(!check_if_halo_is_descendant(halos_i[i_halo].bridge_forematch_default.halo,&(halos_j[j_halo]),n_search)){
if(flag_report_i || flag_report_j) fprintf(stderr,"REPORT:  -> testB5\n");
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
   }
   else{
if(flag_report_i || flag_report_j) fprintf(stderr,"REPORT:  -> testC\n");
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

if(flag_report_i){
char *halo_type_string_1=NULL;
char *halo_type_string_2=NULL;
tree_case_flags_text(halos_i[i_halo].type,"+",&halo_type_string_1);
tree_case_flags_text(halos_j[j_halo].type,"+",&halo_type_string_2);
fprintf(stderr,"REPORT: I-B[%d][%d]->[%d][%d]: %d %d {%s} {%s} -- %d\n",halos_i[i_halo].snap,i_halo,halos_j[j_halo].snap,j_halo,halos_i[i_halo].n_back_matches,halos_j[j_halo].n_back_matches,halo_type_string_1,halo_type_string_2,flag_check_emerged);
SID_free(SID_FARG halo_type_string_1);
SID_free(SID_FARG halo_type_string_2);
}
if(flag_report_j){
char *halo_type_string_1=NULL;
char *halo_type_string_2=NULL;
tree_case_flags_text(halos_i[i_halo].type,"+",&halo_type_string_1);
tree_case_flags_text(halos_j[j_halo].type,"+",&halo_type_string_2);
fprintf(stderr,"REPORT: J-B[%d][%d]->[%d][%d]: %d %d {%s} {%s} -- %d\n",halos_i[i_halo].snap,i_halo,halos_j[j_halo].snap,j_halo,halos_i[i_halo].n_back_matches,halos_j[j_halo].n_back_matches,halo_type_string_1,halo_type_string_2,flag_check_emerged);
SID_free(SID_FARG halo_type_string_1);
SID_free(SID_FARG halo_type_string_2);
}
}

