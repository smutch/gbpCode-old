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
      double M_vir   = fetch_treenode_Mvir(trees,halo);
      double R_vir   = fetch_treenode_Rvir(trees,halo);
      double sigma_v = fetch_treenode_sigmav(trees,halo);
      double r_s     = fetch_treenode_c_NFW(trees,halo);
      double c       = R_vir/r_s;
      double g       = 1./(take_log10(1.+c)-c/(1.+c));
      double pot_fac = c*g*g*G_NEWTON*M_vir*M_vir/2.0/R_vir;
      double pot     = pot_fac*(1.0-1.0/(1.0+c)/(1.0+c)-2*take_log10(1.0+c)/(1.0+c));
      double kin     = 0.5*M_vir*sigma_v*sigma_v;
      virial_ratio   = 2*kin/pot;
//if(M_vir>1e8 && M_vir<3e8) printf("%10.5le %10.5lf %10.5lf -- %10.5le -- %10.3le %10.3le %10.5lf\n",M_vir,R_vir,sigma_v,c,pot,kin,virial_ratio);
   }
   return(virial_ratio);
}

