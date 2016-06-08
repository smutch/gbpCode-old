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

// We want to measure merger ratios (zeta) using the particle counts at the time when
//    the secondary reaches its peak size.  Find the halos involved at that time.  Because
//    of skips in the trees, etc. the halos may not be *exactly* contemporaneous, but we can
//    usually tolerate that since the change of *peak* mass accross skips should be small.
int  fetch_treenode_merger_info(tree_info *trees,tree_node_info **halo_secondary,tree_node_info **halo_primary,
                                double *zeta,double *v_rel,double *sig_v_p,double *sig_v_s){

   // Sanity checks
   if((*halo_primary)==NULL)
      SID_trap_error("Primary halo passed to fetch_treenode_merger_info() is not defined.",ERROR_LOGIC);
   if((*halo_secondary)==NULL)
      SID_trap_error("Secondary halo passed to fetch_treenode_merger_info() is not defined.",ERROR_LOGIC);
   //if(!check_mode_for_flag((*halo_primary)->tree_case,TREE_CASE_MERGER_PRIMARY)) 
   //   SID_trap_error("Primary halo passed to fetch_treenode_merger_info() is not marked with TREE_CASE_MERGER_PRIMARY.",ERROR_LOGIC);
   //if(!check_mode_for_flag((*halo_secondary)->tree_case,TREE_CASE_MERGER)) 
   //   SID_trap_error("Secondary halo passed to fetch_treenode_merger_info() is not marked with TREE_CASE_MERGER.",ERROR_LOGIC);

   // Find peak-mass secondary and primary
   int flag_success=FALSE;
   tree_markers_info *markers_secondary=fetch_treenode_precomputed_markers(trees,(*halo_secondary));
   (*halo_secondary)                   =markers_secondary->peak_mass;
   flag_success=find_treenode_snap_equals_given(trees,
                                                (*halo_primary),
                                                (*halo_secondary)->snap_tree,
                                                halo_primary,
                                                TREE_PROGENITOR_ORDER_N_PARTICLES_INCLUSIVE_PEAK);

   // Compute merger ratio (zeta)
   if((*halo_secondary)!=NULL && (*halo_primary)!=NULL){
      if(zeta!=NULL){
         // Fetch peak particle counts for both halos
         int n_p_peak_secondary=(*halo_secondary)->n_particles_inclusive_peak;
         int n_p_peak_primary  =(*halo_primary)->n_particles_inclusive_peak;

         // Make sure zeta<=1
         int n_p_lo=MIN(n_p_peak_primary,n_p_peak_secondary);
         int n_p_hi=MAX(n_p_peak_primary,n_p_peak_secondary);
         (*zeta)=(double)n_p_lo/(double)n_p_hi;
      }

      // Compute relative velocity
      if(v_rel!=NULL){
         double vx_primary  =fetch_treenode_vx(trees,(*halo_primary));
         double vy_primary  =fetch_treenode_vy(trees,(*halo_primary));
         double vz_primary  =fetch_treenode_vz(trees,(*halo_primary));
         double vx_secondary=fetch_treenode_vx(trees,(*halo_secondary));
         double vy_secondary=fetch_treenode_vy(trees,(*halo_secondary));
         double vz_secondary=fetch_treenode_vz(trees,(*halo_secondary));
         (*v_rel)    =sqrt(pow(vx_secondary-vx_primary,2.)+pow(vy_secondary-vy_primary,2.)+pow(vz_secondary-vz_primary,2.));
      }

      // Fetch velocity dispersions
      if(sig_v_p!=NULL) (*sig_v_p)=fetch_treenode_sigmav(trees,(*halo_primary));
      if(sig_v_s!=NULL) (*sig_v_s)=fetch_treenode_sigmav(trees,(*halo_secondary));
   }
   else{
      flag_success=FALSE;
      if(zeta!=NULL)    (*zeta)   =-1;
      if(v_rel!=NULL)   (*v_rel)  =0.;
      if(sig_v_p!=NULL) (*sig_v_p)=0.;
      if(sig_v_s!=NULL) (*sig_v_s)=0.;
   }

   return(flag_success);
}

