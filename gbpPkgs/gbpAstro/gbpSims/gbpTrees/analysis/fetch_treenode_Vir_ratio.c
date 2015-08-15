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

double fetch_treenode_Vir_ratio(tree_info *trees,tree_node_info *halo){
   double virial_ratio=-1.;
   if(halo!=NULL){
      double R_vir=fetch_treenode_Rvir (trees,halo);
      double r_s  =fetch_treenode_c_NFW(trees,halo);
      double c    =R_vir/r_s;
      if(c>0. && r_s>0.){
         R_vir*=M_PER_MPC;
         r_s  *=M_PER_MPC;
         double M_vir  =fetch_treenode_Mvir(trees,halo)*M_SOL;
         double sigma_v=1e3*fetch_treenode_sigmav(trees,halo);
         double g      =1./(take_ln(1.+c)-c/(1.+c));
         double pot_fac=c*g*g*G_NEWTON*M_vir*M_vir/2.0/R_vir;
         double pot    =pot_fac*(1.0-1.0/(1.0+c)/(1.0+c)-2*take_ln(1.0+c)/(1.0+c));
         double kin    =0.5*M_vir*sigma_v*sigma_v;
         virial_ratio  =2*kin/pot;
      }
   }
   return(virial_ratio);
}

