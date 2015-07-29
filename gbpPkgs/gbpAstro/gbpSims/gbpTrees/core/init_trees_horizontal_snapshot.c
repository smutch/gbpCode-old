#include <gbpLib.h>
#include <gbpTrees_build.h>

void init_trees_horizontal_snapshot(tree_horizontal_info *halos_i,int n_halos_i,int i_read,int i_file,int n_halos_max){
   int i_halo;
   for(i_halo=0;i_halo<n_halos_max;i_halo++){
      halos_i[i_halo].file                          =  i_file;
      halos_i[i_halo].snap                          =  i_read;
      halos_i[i_halo].index                         =  (size_t)i_halo;
      halos_i[i_halo].n_back_matches                =    0;
      halos_i[i_halo].descendant.halo               = NULL;
      halos_i[i_halo].descendant.score              =   0.;
      halos_i[i_halo].first_progenitor.halo         = NULL;
      halos_i[i_halo].first_progenitor.score        =   0.;
      halos_i[i_halo].first_progenitor.flag_two_way =FALSE;
      halos_i[i_halo].last_progenitor.halo          = NULL;
      halos_i[i_halo].last_progenitor.score         =   0.;
      halos_i[i_halo].last_progenitor.flag_two_way  =FALSE;
      halos_i[i_halo].next_progenitor.halo          = NULL;
      halos_i[i_halo].next_progenitor.score         =   0.;
      halos_i[i_halo].next_progenitor.flag_two_way  =FALSE;
      halos_i[i_halo].forematch_first.halo          = NULL;
      halos_i[i_halo].forematch_first.score         =   0.;
      halos_i[i_halo].forematch_first.flag_two_way  =FALSE;
      halos_i[i_halo].forematch_default.halo        = NULL;
      halos_i[i_halo].forematch_default.score       =   0.;
      halos_i[i_halo].forematch_default.flag_two_way=FALSE;
      halos_i[i_halo].forematch_best.halo           = NULL;
      halos_i[i_halo].forematch_best.score          =   0.;
      halos_i[i_halo].forematch_best.flag_two_way   =FALSE;
      halos_i[i_halo].bridge_backmatch.halo         = NULL;
      halos_i[i_halo].bridge_backmatch.score        =   0.;
      halos_i[i_halo].bridge_backmatch.flag_two_way =FALSE;
      SID_free(SID_FARG halos_i[i_halo].back_matches);
      if(i_halo<n_halos_i)
         halos_i[i_halo].type=TREE_CASE_UNPROCESSED|TREE_CASE_NO_PROGENITORS;
      else
         halos_i[i_halo].type=TREE_CASE_INVALID;
      halos_i[i_halo].id                =-1;
      halos_i[i_halo].main_progenitor_id=-1;
      halos_i[i_halo].tree_id           =-1;
      halos_i[i_halo].n_particles       = 0;
      halos_i[i_halo].n_particles_parent= 0;
      halos_i[i_halo].n_particles_largest_descendant= 0;
      halos_i[i_halo].n_progenitors                 = 0;
   }
}

