#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

// Create some default null function pointers
void process_trees_fctn_init_null     (tree_info *trees,void *params,int mode,int i_type){}
void process_trees_fctn_init_snap_null(tree_info *trees,void *params,int mode,int i_type,int flag_init,int i_snap){}
int  process_trees_fctn_select_null   (tree_info *trees,void *params,int mode,int i_type,int flag_init,tree_node_info *halo){return(TRUE);}
void process_trees_fctn_analyze_null  (tree_info *trees,void *params,int mode,int i_type,int flag_init,tree_node_info *halo){}
void process_trees_fctn_fin_snap_null (tree_info *trees,void *params,int mode,int i_type,int flag_init,int i_snap){}
void process_trees_fctn_fin_null      (tree_info *trees,void *params,int mode,int i_type){}

void process_trees_by_snap(tree_info  *trees,
                           void       *params,
                           int         mode,
                           int         i_snap_lo,
                           int         n_snap_process,
                           void      (*init_fctn)     (tree_info *trees,void *params,int mode,int i_type),
                           void      (*init_snap_fctn)(tree_info *trees,void *params,int mode,int i_type,int flag_init,int i_snap),
                           int       (*select_fctn)   (tree_info *trees,void *params,int mode,int i_type,int flag_init,tree_node_info *halo),
                           void      (*analyze_fctn)  (tree_info *trees,void *params,int mode,int i_type,int flag_init,tree_node_info *halo),
                           void      (*fin_snap_fctn) (tree_info *trees,void *params,int mode,int i_type,int flag_init,int i_snap),
                           void      (*fin_fctn)      (tree_info *trees,void *params,int mode,int i_type)){

   // Perform sanity check on the passed snapshot range
   int i_snap_hi=i_snap_lo+n_snap_process-1;
   if(i_snap_hi<i_snap_lo || i_snap_hi>trees->n_snaps || i_snap_lo<0)
      SID_trap_error("Invalid snapshot range (%d to %d; %d in trees) passed to process_trees_by_snap().",ERROR_LOGIC,
                     i_snap_lo,i_snap_hi,trees->n_snaps);

   // Loop over each halo type in turn
   for(int i_type=0;i_type<2;i_type++){
      int flag_do=TRUE;
      if(i_type==0){
         if(check_mode_for_flag(mode,PROCESS_TREES_GROUPS))
            SID_log("Performing GROUP analysis...",SID_LOG_OPEN|SID_LOG_TIMER);
         else
            flag_do=FALSE;
      }
      else if(i_type==1){
         if(check_mode_for_flag(mode,PROCESS_TREES_SUBGROUPS))
            SID_log("Performing SUBGROUP analysis...",SID_LOG_OPEN|SID_LOG_TIMER);
         else
            flag_do=FALSE;
      }

      if(flag_do){
         // Initialize fctn initializeation flags
         int flag_init_init_snap=TRUE;
         int flag_init_analyze  =TRUE;
         int flag_init_select   =TRUE;
         int flag_init_fin_snap =TRUE;

         // Initialize
         if(init_fctn!=process_trees_fctn_init_null){
            SID_log("Initializing...",SID_LOG_OPEN|SID_LOG_TIMER);
            init_fctn(trees,params,mode,i_type);
            SID_log("Done.",SID_LOG_CLOSE);
         }

         // Loop over each snapshot in turn
         for(int i_snap=i_snap_lo;i_snap<=i_snap_hi;i_snap++){
            SID_log("Processing snapshot #%03d of %03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_snap+1,trees->n_snaps);

            // Initialize redshift
            if(init_snap_fctn!=process_trees_fctn_init_snap_null)
               init_snap_fctn(trees,params,mode,i_type,flag_init_init_snap,i_snap);
            flag_init_init_snap=FALSE;

            // Process redshift
            tree_node_info *current_halo;
            if(i_type==0)
               current_halo=trees->first_neighbour_groups[i_snap];
            else if(i_type==1)
               current_halo=trees->first_neighbour_subgroups[i_snap];
            while(current_halo!=NULL){
               // Perfrom select-and-analyze on each halo here
               if(select_fctn==process_trees_fctn_select_null){
                  analyze_fctn(trees,params,mode,i_type,flag_init_analyze,current_halo);
                  flag_init_analyze=FALSE;
               }
               else if(select_fctn(trees,params,mode,i_type,flag_init_select,current_halo)){
                  analyze_fctn(trees,params,mode,i_type,flag_init_analyze,current_halo);
                  flag_init_analyze=FALSE;
               }
               flag_init_select=FALSE;
               current_halo=current_halo->next_neighbour;
            }

            // Finalize redshift
            if(fin_snap_fctn!=process_trees_fctn_fin_snap_null)
               fin_snap_fctn(trees,params,mode,i_type,flag_init_fin_snap,i_snap);
            flag_init_fin_snap=FALSE;

            SID_log("Done.",SID_LOG_CLOSE);
         }

         // Finalize
         if(fin_fctn!=process_trees_fctn_fin_null){
            SID_log("Finalizing...",SID_LOG_OPEN|SID_LOG_TIMER);
            fin_fctn(trees,params,mode,i_type);
            SID_log("Done.",SID_LOG_CLOSE);
         }
      }
      SID_log("Done.",SID_LOG_CLOSE);
   }
}

