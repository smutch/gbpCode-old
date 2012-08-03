#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpClustering.h>

#define CFUNC_ADD_PAIR_DD   2
#define CFUNC_ADD_PAIR_DR   4
#define CFUNC_ADD_PAIR_RR   8
#define CFUNC_ADD_PAIR_D1D2 16
#define CFUNC_SELF_MATCH    32

double calc_CFUNC_local(double DD,double DR,double RR,cfunc_info *cfunc);
double calc_CFUNC_local(double DD,double DR,double RR,cfunc_info *cfunc){
  if(RR>0)
    return((DD*(double)(cfunc->F_random*cfunc->F_random)-2*DR*(double)(cfunc->F_random)+RR)/RR);
  else
    return(0.);
}

void add_pair_CFUNC_local(double x_i,double y_i,double z_i,
                          double x_j,double y_j,double z_j,
                          int zone_i,int zone_j,
                          int mode, cfunc_info *cfunc);
void add_pair_CFUNC_local(double x_i,double y_i,double z_i,
                          double x_j,double y_j,double z_j,
                          int zone_i,int zone_j,
                          int mode, cfunc_info *cfunc){
   // Compute the separations in the periodic box
   double dx,dy,dz;
   dx =d_periodic(x_i-x_j,cfunc->box_size);
   dy =d_periodic(y_i-y_j,cfunc->box_size);
   dz =d_periodic(z_i-z_j,cfunc->box_size);

   // ... 1D case ...
   double sep_1D;
   sep_1D=sqrt(dx*dx+dy*dy+dz*dz);
   if(sep_1D<cfunc->r_max_1D){

     // Set the arrays we are adding to
     long long **array;
     long long **larray;
     switch(mode){
        case CFUNC_ADD_PAIR_DD:
           array =cfunc->DD_1D;
           larray=cfunc->DD_l1D;
           break;
        case CFUNC_ADD_PAIR_DR:
           array =cfunc->DR_1D;
           larray=cfunc->DR_l1D;
           break;
        case CFUNC_ADD_PAIR_RR:
           array =cfunc->RR_1D;
           larray=cfunc->RR_l1D;
           break;
     }

     // Add to the linear array if we are within bounds
     int bin_1D;
     bin_1D=(int)((sep_1D)/cfunc->dr_1D); //r_min=0 in this case
     if(bin_1D>=0 && bin_1D<cfunc->n_1D){
       int i_jack=0;
       array[i_jack++][bin_1D]++;
       for(;i_jack<=cfunc->n_jack_total;i_jack++){
         if(zone_i!=i_jack && zone_j!=i_jack)
           array[i_jack][bin_1D]++;
       }
     }

     // Add to the logarythmic array if we are within bounds
     int bin_l1D;
     bin_l1D=(int)((take_log10(sep_1D)-cfunc->r_min_l1D)/cfunc->dr_l1D);
     if(bin_l1D>=0 && bin_l1D<cfunc->n_1D){
       int i_jack=0;
       larray[i_jack++][bin_1D]++;
       for(;i_jack<=cfunc->n_jack_total;i_jack++){
         if(zone_i!=i_jack && zone_j!=i_jack)
           larray[i_jack][bin_1D]++;
       }
     }
   }

   // ... 2D case ...
   double sep_2D_x;
   double sep_2D_y;
   int    bin_2D_x;
   int    bin_2D_y;
   int    bin_2D;
   sep_2D_x=sqrt(dx*dx+dy*dy);
   bin_2D_x=(int)((sep_2D_x-cfunc->r_min_2D)/cfunc->dr_2D);
   if(bin_2D_x>=0 && bin_2D_x<cfunc->n_2D){
      sep_2D_y=sqrt(dz*dz);
      bin_2D_y=(int)((sep_2D_x-cfunc->r_min_2D)/cfunc->dr_2D);
      if(bin_2D_y>=0 && bin_2D_y<cfunc->n_2D){

         // Compute bin
         bin_2D_x=(int)((sep_2D_x-cfunc->r_min_2D)/cfunc->dr_2D);
         bin_2D_y=(int)((sep_2D_y-cfunc->r_min_2D)/cfunc->dr_2D);
         bin_2D=bin_2D_y*cfunc->n_2D+bin_2D_x;

         // Set the arrays we are adding to
         long long **array;
         switch(mode){
            case CFUNC_ADD_PAIR_DD:
               array =cfunc->DD_2D;
               break;
            case CFUNC_ADD_PAIR_DR:
               array =cfunc->DR_2D;
               break;
            case CFUNC_ADD_PAIR_RR:
               array =cfunc->RR_2D;
               break;
         } 

         // Add to the array if we are within bounds
         if(bin_2D>=0 && bin_2D<cfunc->n_2D){
            int i_jack=0;
            array[i_jack++][bin_2D]++;
            for(;i_jack<=cfunc->n_jack_total;i_jack++){
               if(zone_i!=i_jack && zone_j!=i_jack)
                  array[i_jack][bin_2D]++;
            }
         }
      }
   }
}
void calc_pairs_local(char       *species_name1,
                      char       *species_name2,
                      int         i_rank,
                      int         i_run,
                      int         mode,
                      plist_info *plist,
                      cfunc_info *cfunc);
void calc_pairs_local(char       *species_name1,
                      char       *species_name2,
                      int         i_rank,
                      int         i_run,
                      int         mode,
                      plist_info *plist,
                      cfunc_info *cfunc){
   SID_log("Processing %s-%s pairs...",SID_LOG_OPEN|SID_LOG_TIMER|SID_LOG_CHECKPOINT,species_name1,species_name2);
   size_t i_data1;
   size_t j_data2;
   size_t index_i;
   size_t index_j;
   size_t j_data1_start=0;
   size_t j_data2_start=0;
   int    i_key_rank;
   size_t PHK_last_local;
   size_t PHK_last_rank;

   // Interpret the mode
   int flag_self_match;
   int flag_pair_type=-1;
   if(check_mode_for_flag(mode,CFUNC_SELF_MATCH))
      flag_self_match=TRUE;
   else
      flag_self_match=FALSE;
   if(check_mode_for_flag(mode,CFUNC_ADD_PAIR_DD)){
      flag_pair_type=CFUNC_ADD_PAIR_DD;
   }
   else if(check_mode_for_flag(mode,CFUNC_ADD_PAIR_DR)){
      if(flag_pair_type>=0) SID_trap_error("Multiple pair types specified.",ERROR_LOGIC);
      flag_pair_type=CFUNC_ADD_PAIR_DR;
   }
   else if(check_mode_for_flag(mode,CFUNC_ADD_PAIR_RR)){
      if(flag_pair_type>=0) SID_trap_error("Multiple pair types specified.",ERROR_LOGIC);
      flag_pair_type=CFUNC_ADD_PAIR_RR;
   }
   else
      SID_trap_error("Unsupported pair type specified.",ERROR_LOGIC);

   // Fetch some stuff we need
   GBPREAL *x_data1_local;
   GBPREAL *y_data1_local;
   GBPREAL *z_data1_local;
   size_t  *PHK_data1_local;
   size_t  *PHK_idx_data1_local;
   size_t   n_data1_local;
   size_t   n_data2_rank;
   int     *zone_data1_local;
   GBPREAL *x_data2_rank;
   GBPREAL *y_data2_rank;
   GBPREAL *z_data2_rank;
   size_t  *PHK_data2_rank;
   size_t  *PHK_idx_data2_rank;
   int     *zone_data2_rank;
   x_data1_local      = (GBPREAL *)ADaPS_fetch(plist->data,"x_%s",        species_name1);
   y_data1_local      = (GBPREAL *)ADaPS_fetch(plist->data,"y_%s",        species_name1);
   z_data1_local      = (GBPREAL *)ADaPS_fetch(plist->data,"z_%s",        species_name1);
   PHK_data1_local    = (size_t  *)ADaPS_fetch(plist->data,"PHK_%s",      species_name1);
   PHK_idx_data1_local= (size_t  *)ADaPS_fetch(plist->data,"PHK_index_%s",species_name1);
   zone_data1_local   = (int     *)ADaPS_fetch(plist->data,"zone_%s",     species_name1);
   if(i_rank==0){
      n_data1_local     =((size_t  *)ADaPS_fetch(plist->data,"n_%s",species_name1))[0];
      n_data2_rank      =((size_t  *)ADaPS_fetch(plist->data,"n_%s",species_name2))[0];
      x_data2_rank      = (GBPREAL *)ADaPS_fetch(plist->data,"x_%s",          species_name2);
      y_data2_rank      = (GBPREAL *)ADaPS_fetch(plist->data,"y_%s",          species_name2);
      z_data2_rank      = (GBPREAL *)ADaPS_fetch(plist->data,"z_%s",          species_name2);
      PHK_data2_rank    = (size_t  *)ADaPS_fetch(plist->data,"PHK_%s",        species_name2);
      PHK_idx_data2_rank= (size_t  *)ADaPS_fetch(plist->data,"PHK_index_%s",  species_name2);
      zone_data2_rank   = (int     *)ADaPS_fetch(plist->data,"zone_%s",       species_name2);
   }
   else{
      n_data1_local     =((size_t  *)ADaPS_fetch(plist->data,"n_boundary_%s",species_name1))[0];
      n_data2_rank      =((size_t  *)ADaPS_fetch(plist->data,"n_boundary_%s",species_name2))[0];
      x_data2_rank      = (GBPREAL *)ADaPS_fetch(plist->data,"x_xchg_%s",              species_name2);
      y_data2_rank      = (GBPREAL *)ADaPS_fetch(plist->data,"y_xchg_%s",              species_name2);
      z_data2_rank      = (GBPREAL *)ADaPS_fetch(plist->data,"z_xchg_%s",              species_name2);
      PHK_data2_rank    = (size_t  *)ADaPS_fetch(plist->data,"PHK_xchg_%s",            species_name2);
      PHK_idx_data2_rank= (size_t  *)ADaPS_fetch(plist->data,"PHK_index_xchg_%s",      species_name2);
      zone_data2_rank   = (int     *)ADaPS_fetch(plist->data,"zone_xchg_%s",           species_name2);
   }

   // Make the z-coordinate the redshift-space coordinate
   GBPREAL *x_data1_local_swap;
   GBPREAL *y_data1_local_swap;
   GBPREAL *z_data1_local_swap;
   GBPREAL *x_data2_rank_swap;
   GBPREAL *y_data2_rank_swap;
   GBPREAL *z_data2_rank_swap;
   x_data1_local_swap=x_data1_local;
   y_data1_local_swap=y_data1_local;
   z_data1_local_swap=z_data1_local;
   x_data2_rank_swap =x_data2_rank;
   y_data2_rank_swap =y_data2_rank;
   z_data2_rank_swap =z_data2_rank;
   switch(i_run){
      case 1:
        x_data1_local=y_data1_local_swap;
        y_data1_local=z_data1_local_swap;
        z_data1_local=x_data1_local_swap;
        x_data2_rank =y_data2_rank_swap;
        y_data2_rank =z_data2_rank_swap;
        z_data2_rank =x_data2_rank_swap;
        break;
      case 2:
        x_data1_local=x_data1_local_swap;
        y_data1_local=z_data1_local_swap;
        z_data1_local=y_data1_local_swap;
        x_data2_rank =x_data2_rank_swap;
        y_data2_rank =z_data2_rank_swap;
        z_data2_rank =y_data2_rank_swap;
        break;
      case 3:
        x_data1_local=x_data1_local_swap;
        y_data1_local=y_data1_local_swap;
        z_data1_local=z_data1_local_swap;
        x_data2_rank =x_data2_rank_swap;
        y_data2_rank =y_data2_rank_swap;
        z_data2_rank =z_data2_rank_swap;
        break;
   }

   PHK_t  *PHK_volume=NULL;
   size_t *index_PHK_volume=NULL;
   int     n_PHK_volume;
   PHK_last_local=PHK_N_KEYS_3D(cfunc->n_bits_PHK)+1;
   for(i_data1=0;i_data1<n_data1_local;i_data1++){
      index_i=PHK_idx_data1_local[i_data1];

      // Determine the PHKs within the local volume of i_data1
      if(PHK_last_local!=PHK_data1_local[index_i]){

         // Change the last key to the new current key
         PHK_last_local=PHK_data1_local[index_i];

         // Determine the local-volume keys
         if(PHK_volume!=NULL)
            SID_free(SID_FARG PHK_volume);
         compute_PHK_volume_keys(cfunc->n_bits_PHK,
                                 (PHK_t)PHK_last_local,
                                 0,cfunc->PHK_width,
                                 &n_PHK_volume,&PHK_volume);
         
         // Compute the indices to the first rank item with each neighbouring key.
         //    Only keep keys that have something in them.
         size_t i_key_rank;
         int    n_key_rank_use;
         if(index_PHK_volume!=NULL)
            SID_free(SID_FARG index_PHK_volume);
         index_PHK_volume=(size_t *)SID_malloc(sizeof(size_t)*n_PHK_volume);
         for(i_key_rank=0,n_key_rank_use=0;i_key_rank<(size_t)n_PHK_volume;i_key_rank++){
            index_PHK_volume[n_key_rank_use]=find_index(PHK_data2_rank,(size_t)PHK_volume[i_key_rank],n_data2_rank,PHK_idx_data2_rank);
            PHK_volume[n_key_rank_use]      =PHK_volume[i_key_rank];
            index_j                         =PHK_idx_data2_rank[index_PHK_volume[n_key_rank_use]];
            if(PHK_data2_rank[index_j]==PHK_volume[i_key_rank]) n_key_rank_use++;
         }
         n_PHK_volume=n_key_rank_use;
      }

      // Check pairs
      for(i_key_rank=0;i_key_rank<n_PHK_volume;i_key_rank++){
         size_t PHK_i;
         PHK_i   =(size_t)PHK_volume[i_key_rank];
         j_data2=index_PHK_volume[i_key_rank];
         index_j =PHK_idx_data2_rank[j_data2];
         if(flag_self_match){
            while(PHK_data2_rank[index_j]==PHK_i){
               if((i_rank==0 && index_i<index_j) || SID.My_rank<((SID.My_rank+i_rank)%SID.n_proc)) // Avoids double counting
                  add_pair_CFUNC_local(x_data1_local[index_i],y_data1_local[index_i],z_data1_local[index_i],
                                       x_data2_rank[index_j], y_data2_rank[index_j], z_data2_rank[index_j],
                                       zone_data1_local[index_i],zone_data2_rank[index_j],
                                       flag_pair_type,
                                       cfunc);
               j_data2++;
               if(j_data2>=n_data2_rank) break;
               index_j=PHK_idx_data2_rank[j_data2];
            }
         }
         else{
            while(PHK_data2_rank[index_j]==PHK_i){
               add_pair_CFUNC_local(x_data1_local[index_i],y_data1_local[index_i],z_data1_local[index_i],
                                    x_data2_rank[index_j], y_data2_rank[index_j], z_data2_rank[index_j],
                                    zone_data1_local[index_i],zone_data2_rank[index_j],
                                    flag_pair_type,
                                    cfunc);
               j_data2++;
               if(j_data2>=n_data2_rank) break;
               index_j=PHK_idx_data2_rank[j_data2];
            }
         }
      }
   }

   // Clean-up
   if(PHK_volume!=NULL)
      SID_free(SID_FARG PHK_volume);
   if(index_PHK_volume!=NULL)
      SID_free(SID_FARG index_PHK_volume);

   SID_log("Done.",SID_LOG_CLOSE);
}

void compute_cfunc(plist_info  *plist,
                   char        *species_name,
                   char        *random_name,
                   cfunc_info  *cfunc,
                   int          i_run){
  int         i_bin,j_bin;
  int         i_jack;
  int         i_x_jack,i_y_jack,i_z_jack;
  int         i_rank;
  int         j_rank;
  size_t      i_data;
  size_t      j_data;
  size_t      n_data_local;
  size_t      n_data;
  size_t      i_random;
  size_t      j_random;
  size_t      n_random_local;
  size_t      n_random;
  size_t      n_data_allocate;
  size_t      n_random_allocate;
  size_t      n_data_rank;
  size_t      n_temp;
  GBPREAL    *x_data_rank;
  GBPREAL    *y_data_rank;
  GBPREAL    *z_data_rank;
  size_t      n_random_rank;
  GBPREAL    *x_random_rank;
  GBPREAL    *y_random_rank;
  GBPREAL    *z_random_rank;
  size_t     *x_data_rank_index;
  size_t     *x_random_rank_index;
  size_t     *PHK_data_rank;
  size_t     *PHK_idx_data_rank;
  size_t     *PHK_bidx_data_rank;
  size_t     *PHK_idx_data_local;
  size_t     *PHK_bidx_data_local;
  size_t     *PHK_random_rank;
  size_t     *PHK_idx_random_rank;
  size_t     *PHK_bidx_random_rank;
  size_t     *PHK_idx_random_local;
  size_t     *PHK_bidx_random_local;
  GBPREAL    *x_random_local_swap;
  GBPREAL    *y_random_local_swap;
  GBPREAL    *z_random_local_swap;
  GBPREAL    *x_data_local_swap;
  GBPREAL    *y_data_local_swap;
  GBPREAL    *z_data_local_swap;
  GBPREAL    *x_data_local_temp;
  GBPREAL    *y_data_local_temp;
  GBPREAL    *z_data_local_temp;
  GBPREAL    *x_data_local;
  GBPREAL    *y_data_local;
  GBPREAL    *z_data_local;
  GBPREAL    *vx_data_local;
  GBPREAL    *vy_data_local;
  GBPREAL    *vz_data_local;
  GBPREAL    *x_random_local;
  GBPREAL    *y_random_local;
  GBPREAL    *z_random_local;
  GBPREAL     x_i,y_i,z_i;
  GBPREAL     x_j,y_j,z_j;
  double      mean,std_dev;
  double      dx,dy,dz;
  double      sep_1D;
  double      sep_2D_x;
  double      sep_2D_y;
  int         bin_1D;
  int         bin_l1D;
  int         bin_2D_x;
  int         bin_2D_y;
  double      d_jack;
  int        *zone_data_local;
  int        *zone_random_local;
  int        *zone_data_rank;
  int        *zone_random_rank;
  size_t     *x_data_local_index;
  size_t     *x_random_local_index;
  size_t      j_data_lo_1;
  size_t      j_data_hi_1;
  size_t      j_data_lo_2;
  size_t      j_data_hi_2;
  size_t      j_random_lo_1;
  size_t      j_random_hi_1;
  size_t      j_random_lo_2;
  size_t      j_random_hi_2;
  int         bin_2D;
  double      r_max;
  double      r_min_1D;
  double      r_max_l1D;
  double      dr_l1D;
  double      dr_1D;
  double      dr_2D;
  GBPREAL     x_temp;
  GBPREAL     y_temp;
  GBPREAL     z_temp;
  long long  *temp_array;
  double      bar_i;
  double      bar_j;
  double      bar_ij;
  PHK_t      *PHK_volume=NULL;
  size_t     *index_PHK_volume=NULL;
  int         n_PHK_volume;

  SID_log("Computing correlation function (%d 1D bins and %d 2D bins)...",SID_LOG_OPEN|SID_LOG_TIMER,cfunc->n_1D,cfunc->n_2D);

  // Parse some stuff from the cfunc structure
  int          n_2D_total;
  int          n_jack_total;
  double       box_size;
  double       r_min_l1D;
  double       r_max_1D;
  double       r_min_2D;
  double       r_max_2D;
  int          n_1D;
  int          n_2D;
  int          n_jack;
  long long  **DD_l1D;
  long long  **DR_l1D;
  long long  **RR_l1D;
  long long  **DD_1D;
  long long  **DR_1D;
  long long  **RR_1D;
  long long  **DD_2D;
  long long  **DR_2D;
  long long  **RR_2D;
  double      *CFUNC_l1D;
  double      *dCFUNC_l1D;
  double      *COVMTX_l1D;
  double      *CFUNC_1D;
  double      *dCFUNC_1D;
  double      *COVMTX_1D;
  double      *CFUNC_2D;
  double      *dCFUNC_2D;
  double      *COVMTX_2D;
  n_2D_total  =cfunc->n_2D_total;
  n_jack_total=cfunc->n_jack_total;
  box_size    =cfunc->box_size;
  r_min_l1D   =cfunc->r_min_l1D;
  r_max_1D    =cfunc->r_max_1D;
  r_min_2D    =cfunc->r_min_2D;
  r_max_2D    =cfunc->r_max_2D;
  n_1D        =cfunc->n_1D;
  n_2D        =cfunc->n_2D;
  n_jack      =cfunc->n_jack;
  DD_l1D      =cfunc->DD_l1D;
  DR_l1D      =cfunc->DR_l1D;
  RR_l1D      =cfunc->RR_l1D;
  DD_1D       =cfunc->DD_1D;
  DR_1D       =cfunc->DR_1D;
  RR_1D       =cfunc->RR_1D;
  DD_2D       =cfunc->DD_2D;
  DR_2D       =cfunc->DR_2D;
  RR_2D       =cfunc->RR_2D;
  CFUNC_l1D   =cfunc->CFUNC_l1D[i_run];
  dCFUNC_l1D  =cfunc->dCFUNC_l1D[i_run];
  COVMTX_l1D  =cfunc->COVMTX_l1D[i_run];
  CFUNC_1D    =cfunc->CFUNC_1D[i_run];
  dCFUNC_1D   =cfunc->dCFUNC_1D[i_run];
  COVMTX_1D   =cfunc->COVMTX_1D[i_run];
  CFUNC_2D    =cfunc->CFUNC_2D[i_run];
  dCFUNC_2D   =cfunc->dCFUNC_2D[i_run];
  COVMTX_2D   =cfunc->COVMTX_2D[i_run];

  // Fetch the positions and object counts
  SID_log("Initializing...",SID_LOG_OPEN);
  n_data        =((size_t  *)ADaPS_fetch(plist->data,"n_all_%s",species_name))[0]; 
  n_data_local  =((size_t  *)ADaPS_fetch(plist->data,"n_%s",    species_name))[0]; 
  x_data_local  = (GBPREAL *)ADaPS_fetch(plist->data,"x_%s",    species_name);
  y_data_local  = (GBPREAL *)ADaPS_fetch(plist->data,"y_%s",    species_name);
  z_data_local  = (GBPREAL *)ADaPS_fetch(plist->data,"z_%s",    species_name);

  n_random      =((size_t  *)ADaPS_fetch(plist->data,"n_all_%s",random_name))[0]; 
  n_random_local=((size_t  *)ADaPS_fetch(plist->data,"n_%s",    random_name))[0]; 
  x_random_local= (GBPREAL *)ADaPS_fetch(plist->data,"x_%s",    random_name);
  y_random_local= (GBPREAL *)ADaPS_fetch(plist->data,"y_%s",    random_name);
  z_random_local= (GBPREAL *)ADaPS_fetch(plist->data,"z_%s",    random_name);

  // Fetch PHK information
  size_t  n_data_boundary;
  size_t  n_random_boundary;
  size_t *PHK_data_local;
  size_t *PHK_random_local;
  size_t *PHK_index_data_local;
  size_t *PHK_index_random_local;
  n_data_boundary      =((size_t *)ADaPS_fetch(plist->data,"n_boundary_%s",        species_name))[0];
  n_random_boundary    =((size_t *)ADaPS_fetch(plist->data,"n_boundary_%s",        random_name))[0];
  PHK_data_local       = (size_t *)ADaPS_fetch(plist->data,"PHK_%s",               species_name);
  PHK_random_local     = (size_t *)ADaPS_fetch(plist->data,"PHK_%s",               random_name);
  PHK_idx_data_local   = (size_t *)ADaPS_fetch(plist->data,"PHK_index_%s",         species_name);
  PHK_bidx_data_local  = (size_t *)ADaPS_fetch(plist->data,"PHK_index_boundary_%s",species_name);
  PHK_idx_random_local = (size_t *)ADaPS_fetch(plist->data,"PHK_index_%s",         random_name);
  PHK_bidx_random_local= (size_t *)ADaPS_fetch(plist->data,"PHK_index_boundary_%s",random_name);

  // Make the z-coordinate the redshift-space coordinate
  x_data_local_swap  =x_data_local;
  y_data_local_swap  =y_data_local;
  z_data_local_swap  =z_data_local;
  x_random_local_swap=x_random_local;
  y_random_local_swap=y_random_local;
  z_random_local_swap=z_random_local;
  switch(i_run){
     case 1:
       x_data_local  =y_data_local_swap;
       y_data_local  =z_data_local_swap;
       z_data_local  =x_data_local_swap;
       x_random_local=y_random_local_swap;
       y_random_local=z_random_local_swap;
       z_random_local=x_random_local_swap;
       break;
     case 2:
       x_data_local  =x_data_local_swap;
       y_data_local  =z_data_local_swap;
       z_data_local  =y_data_local_swap;
       x_random_local=x_random_local_swap;
       y_random_local=z_random_local_swap;
       z_random_local=y_random_local_swap;
       break;
     case 3:
       x_data_local  =x_data_local_swap;
       y_data_local  =y_data_local_swap;
       z_data_local  =z_data_local_swap;
       x_random_local=x_random_local_swap;
       y_random_local=y_random_local_swap;
       z_random_local=z_random_local_swap;
       break;
  }

  // Determine which jacknife zone each object belongs to
  d_jack=box_size/(double)n_jack;
  zone_data_local=(int *)SID_malloc(sizeof(int)*n_data_local);
  for(i_data=0;i_data<n_data_local;i_data++){
    i_x_jack=(int)(x_data_local[i_data]/d_jack);
    i_y_jack=(int)(y_data_local[i_data]/d_jack);
    i_z_jack=(int)(z_data_local[i_data]/d_jack);
    if(i_x_jack<0)       i_x_jack=0;
    if(i_x_jack>=n_jack) i_x_jack=n_jack-1;
    if(i_y_jack<0)       i_y_jack=0;
    if(i_y_jack>=n_jack) i_y_jack=n_jack-1;
    if(i_z_jack<0)       i_z_jack=0;
    if(i_z_jack>=n_jack) i_z_jack=n_jack-1;
    zone_data_local[i_data]=(i_x_jack*n_jack*n_jack)+i_y_jack*n_jack+i_z_jack;
  }
  zone_random_local=(int *)SID_malloc(sizeof(int)*n_random_local);
  for(i_random=0;i_random<n_random_local;i_random++){
    i_x_jack=(int)(x_random_local[i_random]/d_jack);
    i_y_jack=(int)(y_random_local[i_random]/d_jack);
    i_z_jack=(int)(z_random_local[i_random]/d_jack);
    if(i_x_jack<0)       i_x_jack=0;
    if(i_x_jack>=n_jack) i_x_jack=n_jack-1;
    if(i_y_jack<0)       i_y_jack=0;
    if(i_y_jack>=n_jack) i_y_jack=n_jack-1;
    if(i_z_jack<0)       i_z_jack=0;
    if(i_z_jack>=n_jack) i_z_jack=n_jack-1;
    zone_random_local[i_random]=(i_x_jack*n_jack*n_jack)+i_y_jack*n_jack+i_z_jack;
  }

  // Store the zones
  ADaPS_store(&(plist->data),(void *)zone_data_local,  "zone_%s",ADaPS_DEFAULT,species_name);
  ADaPS_store(&(plist->data),(void *)zone_random_local,"zone_%s",ADaPS_DEFAULT,random_name);

  // Initialize the pair count arrays
  if(cfunc->flag_compute_RR){
    for(i_jack=0;i_jack<=n_jack_total;i_jack++){
      for(i_bin=0;i_bin<n_1D;i_bin++){
        RR_1D[i_jack][i_bin] =0;
        RR_l1D[i_jack][i_bin]=0;
      }
      for(i_bin=0;i_bin<n_2D*n_2D;i_bin++)
        RR_2D[i_jack][i_bin]=0;
    }
  }
  for(i_jack=0;i_jack<=n_jack_total;i_jack++){
    for(i_bin=0;i_bin<n_1D;i_bin++){
      DD_1D[i_jack][i_bin] =0;
      DD_l1D[i_jack][i_bin]=0;
      DR_1D[i_jack][i_bin] =0;
      DR_l1D[i_jack][i_bin]=0;
    }
    for(i_bin=0;i_bin<n_2D*n_2D;i_bin++){
      DD_2D[i_jack][i_bin]=0;
      DR_2D[i_jack][i_bin]=0;
    }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Loop over all the ranks
  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
    if(SID.n_proc>1)
      SID_log("Processing rank %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_rank+1,SID.n_proc);

    // For the first iteration, process local pairs ...
    if(i_rank==0){
       n_data_rank        =n_data_local;
       x_data_rank        =x_data_local;
       y_data_rank        =y_data_local;
       z_data_rank        =z_data_local;
       PHK_data_rank      =PHK_data_local;
       PHK_idx_data_rank  =PHK_idx_data_local;
       n_random_rank      =n_random_local;
       x_random_rank      =x_random_local;
       y_random_rank      =y_random_local;
       z_random_rank      =z_random_local;
       PHK_random_rank    =PHK_random_local;
       PHK_idx_random_rank=PHK_idx_random_local;
       zone_data_rank     =zone_data_local;
       zone_random_rank   =zone_random_local;
       x_data_rank_index  =x_data_local_index;
       x_random_rank_index=x_random_local_index;
    }
    // ... subsequently, process pairs between boundaries.
    else{
      if(i_rank==1){
        // We only need to work with items on the boundaries now.  Since they are all at
        //   the beginning of the arrays, we can just pretend that the arrays are shorter.
        n_data_local  =n_data_boundary;
        n_random_local=n_random_boundary;

        // Determine the max size of exchanges and report results
        SID_log("Communicating array sizes...",SID_LOG_OPEN|SID_LOG_CHECKPOINT);
        SID_Allreduce(&n_data_local,  &n_data_allocate,  1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
        SID_log("n_data_max  =%lld",SID_LOG_COMMENT,n_data_allocate);
        SID_Allreduce(&n_random_local,&n_random_allocate,1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
        SID_log("n_random_max=%lld",SID_LOG_COMMENT,n_random_allocate);
        SID_log("Done.",SID_LOG_CLOSE);
        SID_log("Allocating arrays...",SID_LOG_OPEN|SID_LOG_CHECKPOINT);

        // Allocate exchange buffers
        x_data_rank        =(GBPREAL *)SID_calloc(sizeof(GBPREAL)*n_data_allocate);
        y_data_rank        =(GBPREAL *)SID_calloc(sizeof(GBPREAL)*n_data_allocate);
        z_data_rank        =(GBPREAL *)SID_calloc(sizeof(GBPREAL)*n_data_allocate);
        PHK_data_rank      =(size_t  *)SID_calloc(sizeof(size_t) *n_data_allocate);
        PHK_idx_data_rank  =(size_t  *)SID_calloc(sizeof(size_t) *n_data_allocate);
        zone_data_rank     =(int     *)SID_calloc(sizeof(int)    *n_data_allocate);
        x_random_rank      =(GBPREAL *)SID_calloc(sizeof(GBPREAL)*n_random_allocate);
        y_random_rank      =(GBPREAL *)SID_calloc(sizeof(GBPREAL)*n_random_allocate);
        z_random_rank      =(GBPREAL *)SID_calloc(sizeof(GBPREAL)*n_random_allocate);
        PHK_random_rank    =(size_t  *)SID_calloc(sizeof(size_t) *n_random_allocate);
        PHK_idx_random_rank=(size_t  *)SID_calloc(sizeof(size_t) *n_random_allocate);
        zone_random_rank   =(int     *)SID_calloc(sizeof(int)    *n_random_allocate);

        // Store the buffers
        ADaPS_store(&(plist->data),(void *)x_data_rank,        "x_xchg_%s",        ADaPS_DEFAULT,species_name);
        ADaPS_store(&(plist->data),(void *)y_data_rank,        "y_xchg_%s",        ADaPS_DEFAULT,species_name);
        ADaPS_store(&(plist->data),(void *)z_data_rank,        "z_xchg_%s",        ADaPS_DEFAULT,species_name);
        ADaPS_store(&(plist->data),(void *)PHK_data_rank,      "PHK_xchg_%s",      ADaPS_DEFAULT,species_name);
        ADaPS_store(&(plist->data),(void *)PHK_idx_data_rank,  "PHK_index_xchg_%s",ADaPS_DEFAULT,species_name);
        ADaPS_store(&(plist->data),(void *)zone_data_rank,     "zone_xchg_%s",     ADaPS_DEFAULT,species_name);
        ADaPS_store(&(plist->data),(void *)x_random_rank,      "x_xchg_%s",        ADaPS_DEFAULT,random_name);
        ADaPS_store(&(plist->data),(void *)y_random_rank,      "y_xchg_%s",        ADaPS_DEFAULT,random_name);
        ADaPS_store(&(plist->data),(void *)z_random_rank,      "z_xchg_%s",        ADaPS_DEFAULT,random_name);
        ADaPS_store(&(plist->data),(void *)PHK_random_rank,    "PHK_xchg_%s",      ADaPS_DEFAULT,random_name);
        ADaPS_store(&(plist->data),(void *)PHK_idx_random_rank,"PHK_index_xchg_%s",ADaPS_DEFAULT,random_name);
        ADaPS_store(&(plist->data),(void *)zone_random_rank,   "zone_xchg_%s",     ADaPS_DEFAULT,random_name);

        SID_log("Done.",SID_LOG_CLOSE);
      }
      // Perform exchange
      SID_log("Performing exchange...",SID_LOG_OPEN|SID_LOG_TIMER);
      exchange_ring_buffer(x_data_local,         sizeof(GBPREAL),n_data_local,  x_data_rank,        &n_data_rank,  i_rank);
      exchange_ring_buffer(y_data_local,         sizeof(GBPREAL),n_data_local,  y_data_rank,        &n_data_rank,  i_rank);
      exchange_ring_buffer(z_data_local,         sizeof(GBPREAL),n_data_local,  z_data_rank,        &n_data_rank,  i_rank);
      exchange_ring_buffer(PHK_data_local,       sizeof(size_t), n_data_local,  PHK_data_rank,      &n_data_rank,  i_rank);
      exchange_ring_buffer(PHK_bidx_data_local,  sizeof(size_t), n_data_local,  PHK_idx_data_rank,  &n_data_rank,  i_rank);
      exchange_ring_buffer(zone_data_local,      sizeof(int),    n_data_local,  zone_data_rank,     &n_data_rank,  i_rank);
      exchange_ring_buffer(x_random_local,       sizeof(GBPREAL),n_random_local,x_random_rank,      &n_random_rank,i_rank);
      exchange_ring_buffer(y_random_local,       sizeof(GBPREAL),n_random_local,y_random_rank,      &n_random_rank,i_rank);
      exchange_ring_buffer(z_random_local,       sizeof(GBPREAL),n_random_local,z_random_rank,      &n_random_rank,i_rank);
      exchange_ring_buffer(PHK_random_local,     sizeof(size_t), n_random_local,PHK_random_rank,    &n_random_rank,i_rank);
      exchange_ring_buffer(PHK_bidx_random_local,sizeof(size_t), n_random_local,PHK_idx_random_rank,&n_random_rank,i_rank);
      exchange_ring_buffer(zone_random_local,    sizeof(int),    n_random_local,zone_random_rank,   &n_random_rank,i_rank);
      SID_log("Done.",SID_LOG_CLOSE);
    }

    // Compute random-random pairs (only done for the first call)
    if(cfunc->flag_compute_RR)
       calc_pairs_local(random_name,random_name,i_rank,i_run,CFUNC_SELF_MATCH|CFUNC_ADD_PAIR_RR,plist,cfunc);

    // Compute data-data pairs
    calc_pairs_local(species_name,species_name,i_rank,i_run,CFUNC_SELF_MATCH|CFUNC_ADD_PAIR_DD,plist,cfunc);

    // Compute data-random pairs
    calc_pairs_local(species_name,random_name,i_rank,i_run,CFUNC_ADD_PAIR_DR,plist,cfunc);

    if(SID.n_proc>1)
      SID_log("Done.",SID_LOG_CLOSE);
  } // i_rank

  if(SID.n_proc>1){
    SID_log("Combining results from separate ranks...",SID_LOG_OPEN);
    for(i_jack=0;i_jack<=n_jack_total;i_jack++){
      SID_Allreduce(SID_IN_PLACE,DD_l1D[i_jack],n_1D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,DR_l1D[i_jack],n_1D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD); 
      SID_Allreduce(SID_IN_PLACE,DD_1D[i_jack], n_1D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,DR_1D[i_jack], n_1D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
      if(cfunc->flag_compute_RR){
        SID_Allreduce(SID_IN_PLACE,RR_l1D[i_jack],n_1D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
        SID_Allreduce(SID_IN_PLACE,RR_1D[i_jack], n_1D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
      }
      SID_Allreduce(SID_IN_PLACE,DD_2D[i_jack],n_2D*n_2D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,DR_2D[i_jack],n_2D*n_2D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
      if(cfunc->flag_compute_RR)
        SID_Allreduce(SID_IN_PLACE,RR_2D[i_jack],n_2D*n_2D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }

  // Compile correlation functions from data and random pair counts
  SID_log("Processing correlation functions...",SID_LOG_OPEN);
  for(i_bin=0;i_bin<n_1D;i_bin++){
    CFUNC_1D[i_bin] =calc_CFUNC_local((double)DD_1D[0][i_bin]/(double)((n_data)*(n_data-1)),
                                      (double)DR_1D[0][i_bin]/(double)((n_data)*(n_random)),
                                      (double)RR_1D[0][i_bin]/(double)((n_random)*(n_random-1)),cfunc);
    CFUNC_l1D[i_bin]=calc_CFUNC_local((double)DD_l1D[0][i_bin]/(double)((n_data)*(n_data-1)),
                                      (double)DR_l1D[0][i_bin]/(double)((n_data)*(n_random)),
                                      (double)RR_l1D[0][i_bin]/(double)((n_random)*(n_random-1)),cfunc);
if(SID.I_am_Master) fprintf(stderr,"%10.3le %7lld %7lld %7lld %10.3le\n",(double)(i_bin*cfunc->dr_1D),DD_1D[0][i_bin],DR_1D[0][i_bin],RR_1D[0][i_bin],CFUNC_1D[i_bin]);
  }
  for(i_bin=0;i_bin<n_2D*n_2D;i_bin++)
    CFUNC_2D[i_bin]=calc_CFUNC_local((double)DD_2D[0][i_bin]/(double)((n_data)*(n_data-1)),
                                     (double)DR_2D[0][i_bin]/(double)((n_data)*(n_random)),
                                     (double)RR_2D[0][i_bin]/(double)((n_random)*(n_random-1)),cfunc);
  SID_log("Done.",SID_LOG_CLOSE);

  // Generate 1D covariance matrices
  SID_log("Generating covariance matrices...",SID_LOG_OPEN);
  // ... process log-space profile first ...
  for(i_bin=0;i_bin<n_1D;i_bin++){
    for(i_jack=1,bar_i=0.;i_jack<=n_jack_total;i_jack++){
      bar_i+=calc_CFUNC_local((double)DD_l1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                              (double)DR_l1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                              (double)RR_l1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)),cfunc);
    }
    bar_i/=(double)n_jack_total;
    for(j_bin=0;j_bin<n_1D;j_bin++){
      for(i_jack=1,bar_j=0.,bar_ij=0.;i_jack<=n_jack_total;i_jack++){
        bar_j +=calc_CFUNC_local((double)DD_l1D[i_jack][j_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_l1D[i_jack][j_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_l1D[i_jack][j_bin]/(double)((n_random)*(n_random-1)),cfunc);
        bar_ij+=calc_CFUNC_local((double)DD_l1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_l1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_l1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)),cfunc)*
                calc_CFUNC_local((double)DD_l1D[i_jack][j_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_l1D[i_jack][j_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_l1D[i_jack][j_bin]/(double)((n_random)*(n_random-1)),cfunc);
      }
      bar_j /=(double)n_jack_total;
      bar_ij/=(double)n_jack_total;
      COVMTX_l1D[i_bin*n_1D+j_bin]=(double)(n_jack_total-1)*(bar_ij-bar_i*bar_j);
    }
  }
  // ... then process the linear-space profile ...
  for(i_bin=0;i_bin<n_1D;i_bin++){
    for(i_jack=1,bar_i=0.;i_jack<=n_jack_total;i_jack++){
      bar_i+=calc_CFUNC_local((double)DD_1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                              (double)DR_1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                              (double)RR_1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)),cfunc);
    }
    bar_i/=(double)n_jack_total;
    for(j_bin=0;j_bin<n_1D;j_bin++){
      for(i_jack=1,bar_j=0.,bar_ij=0.;i_jack<=n_jack_total;i_jack++){
        bar_j +=calc_CFUNC_local((double)DD_1D[i_jack][j_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_1D[i_jack][j_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_1D[i_jack][j_bin]/(double)((n_random)*(n_random-1)),cfunc);
        bar_ij+=calc_CFUNC_local((double)DD_1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)),cfunc)*
                calc_CFUNC_local((double)DD_1D[i_jack][j_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_1D[i_jack][j_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_1D[i_jack][j_bin]/(double)((n_random)*(n_random-1)),cfunc);
      }
      bar_j /=(double)n_jack_total;
      bar_ij/=(double)n_jack_total;
      COVMTX_1D[i_bin*n_1D+j_bin]=(double)(n_jack_total-1)*(bar_ij-bar_i*bar_j);
    }
  }
  // ... and finally, the 2D case ...
  for(i_bin=0;i_bin<n_2D_total;i_bin++){
    for(i_jack=1,bar_i=0.;i_jack<=n_jack_total;i_jack++){
      bar_i+=calc_CFUNC_local((double)DD_2D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                              (double)DR_2D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                              (double)RR_2D[i_jack][i_bin]/(double)((n_random)*(n_random-1)),cfunc);
    }
    bar_i/=(double)n_jack_total;
    for(j_bin=0;j_bin<n_2D_total;j_bin++){
      for(i_jack=1,bar_j=0.,bar_ij=0.;i_jack<=n_jack_total;i_jack++){
        bar_j +=calc_CFUNC_local((double)DD_2D[i_jack][j_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_2D[i_jack][j_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_2D[i_jack][j_bin]/(double)((n_random)*(n_random-1)),cfunc);
        bar_ij+=calc_CFUNC_local((double)DD_2D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_2D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_2D[i_jack][i_bin]/(double)((n_random)*(n_random-1)),cfunc)*
                calc_CFUNC_local((double)DD_2D[i_jack][j_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_2D[i_jack][j_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_2D[i_jack][j_bin]/(double)((n_random)*(n_random-1)),cfunc);
      }
      bar_j /=(double)n_jack_total;
      bar_ij/=(double)n_jack_total;
      COVMTX_2D[i_bin*n_2D_total+j_bin]=(double)(n_jack_total-1)*(bar_ij-bar_i*bar_j);
    }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Compile uncertainties from jacknife-resampling data and random pair counts
  SID_log("Processing jack knifes...",SID_LOG_OPEN);
  // ... 1D case ...
  for(i_bin=0;i_bin<n_1D;i_bin++){
    for(i_jack=1,mean=0.;i_jack<=n_jack_total;i_jack++)
      mean+=calc_CFUNC_local((double)DD_l1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                             (double)DR_l1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                             (double)RR_l1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)),cfunc);    
    mean/=(double)n_jack_total;
    for(i_jack=1,std_dev=0.;i_jack<=n_jack_total;i_jack++)
      std_dev+=pow(mean-calc_CFUNC_local((double)DD_l1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                                         (double)DR_l1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                                         (double)RR_l1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)),cfunc),2.);
    std_dev=sqrt(std_dev/(double)n_jack_total);
    dCFUNC_l1D[i_bin]=std_dev*sqrt((double)(n_jack_total-1));
    for(i_jack=1,mean=0.;i_jack<=n_jack_total;i_jack++)
      mean+=calc_CFUNC_local((double)DD_1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                             (double)DR_1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                             (double)RR_1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)),cfunc);    
    mean/=(double)n_jack_total;
    for(i_jack=1,std_dev=0.;i_jack<=n_jack_total;i_jack++)
      std_dev+=pow(mean-calc_CFUNC_local((double)DD_1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                                         (double)DR_1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                                         (double)RR_1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)),cfunc),2.);
    std_dev=sqrt(std_dev/(double)n_jack_total);
    dCFUNC_1D[i_bin]=std_dev*sqrt((double)(n_jack_total-1));
  }

  // ... 2D case ...
  for(i_bin=0;i_bin<n_2D*n_2D;i_bin++){
    for(i_jack=1,mean=0.;i_jack<=n_jack_total;i_jack++)
      mean+=calc_CFUNC_local((double)DD_2D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                             (double)DR_2D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                             (double)RR_2D[i_jack][i_bin]/(double)((n_random)*(n_random-1)),cfunc);    
    mean/=(double)n_jack_total;
    for(i_jack=1,std_dev=0.;i_jack<=n_jack_total;i_jack++)
      std_dev+=pow(mean-calc_CFUNC_local((double)DD_2D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                                         (double)DR_2D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                                         (double)RR_2D[i_jack][i_bin]/(double)((n_random)*(n_random-1)),cfunc),2.);
    std_dev=sqrt(std_dev/(double)n_jack_total);
    dCFUNC_2D[i_bin]=std_dev*sqrt((double)(n_jack_total-1));
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Don't recompute the RR's unless otherwise instructed
  cfunc->flag_compute_RR=FALSE;  

  SID_log("Done.",SID_LOG_CLOSE);
}


