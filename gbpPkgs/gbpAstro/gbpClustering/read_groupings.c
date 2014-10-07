#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>
#include <gbpPHKs.h>

#define N_BITS_MIN 1

void read_groupings(const char *filename_root,int grouping_number,plist_info *plist,int mode,...){
   // Interpret variable arguments
   double      box_size;
   double      redshift;
   cosmo_info *cosmo;
   double      h_Hubble;
   va_list     vargs;
   va_start(vargs,mode);

   // Interpret mode ...

   // ... domain decomposition ...
   double     r_max;
   int        n_bits_PHK;
   int        PHK_width;
   slab_info *slab;
   if(check_mode_for_flag(mode,READ_GROUPING_SLAB) && check_mode_for_flag(mode,READ_GROUPING_PHK))
      SID_trap_error("Multiple domain decompositions have been set in read_groupings().",ERROR_LOGIC);
   if(check_mode_for_flag(mode,READ_GROUPING_SLAB))
      slab =(slab_info *)va_arg(vargs,slab_info *);
   else if(check_mode_for_flag(mode,READ_GROUPING_PHK)){
      r_max     =(double)va_arg(vargs,double);
      box_size  =(double)va_arg(vargs,double);
      n_bits_PHK=(int)   va_arg(vargs,int);
      PHK_width =(int)   va_arg(vargs,int);
   }
   else
      SID_trap_error("No domain decomposition has been set in read_groupings().",ERROR_LOGIC);

   // ... redshift space? ...
   int    flag_add_zspace_x;
   int    flag_add_zspace_y;
   int    flag_add_zspace_z;
   int    flag_add_zspace;
   int    n_z_dims;
   flag_add_zspace_x=check_mode_for_flag(mode,READ_GROUPING_ADD_VX);
   flag_add_zspace_y=check_mode_for_flag(mode,READ_GROUPING_ADD_VY);
   flag_add_zspace_z=check_mode_for_flag(mode,READ_GROUPING_ADD_VZ);
   n_z_dims=flag_add_zspace_x+flag_add_zspace_y+flag_add_zspace_z;
   if(n_z_dims>1)
      SID_trap_error("More then one redshift-space dimension (%d) has been specified.",ERROR_LOGIC,n_z_dims);
   else if(n_z_dims==1){
      flag_add_zspace=TRUE;
      if(check_mode_for_flag(mode,READ_GROUPING_SLAB))
         box_size=(double)va_arg(vargs,double);
      redshift   =va_arg(vargs,double);
      cosmo      =va_arg(vargs,cosmo_info *);
      h_Hubble   =((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
   }

   // Set filename, open file and count the number of items
   char   filename_in[MAX_FILENAME_LENGTH];
   char   filename_name[MAX_FILENAME_LENGTH];
   FILE  *fp_in;
   size_t n_halos;
   sprintf(filename_in,"%s_grouping_%03d.dat",filename_root,grouping_number);
   strcpy(filename_name,filename_in);
   strip_path(filename_name);
   SID_log("Reading grouping catalog {%s}...",SID_LOG_OPEN,filename_name);
   fp_in  =fopen(filename_in,"r");
   n_halos=(size_t)count_lines_data(fp_in);

   // Set the columns of the file we're reading
   int x_column     =10;
   int y_column     =11;
   int z_column     =12;
   int vx_sub_column=14;
   int vy_sub_column=15;
   int vz_sub_column=16;
   int vx_FoF_column=17;
   int vy_FoF_column=18;
   int vz_FoF_column=19;
   int vx_sys_column=20;
   int vy_sys_column=21;
   int vz_sys_column=22;
   int vx_column    =vx_sub_column; // Used for z-space distortions
   int vy_column    =vy_sub_column; // Used for z-space distortions
   int vz_column    =vz_sub_column; // Used for z-space distortions

   // Set initial values for position bounds of the objects we will read
   GBPREAL x_min;
   GBPREAL x_max;
   GBPREAL y_min;
   GBPREAL y_max;
   GBPREAL z_min;
   GBPREAL z_max;
   GBPREAL vx_min;
   GBPREAL vx_max;
   GBPREAL vy_min;
   GBPREAL vy_max;
   GBPREAL vz_min;
   GBPREAL vz_max;
   x_min = 1e10;
   x_max =-1e10;
   y_min = 1e10;
   y_max =-1e10;
   z_min = 1e10;
   z_max =-1e10;
   vx_min= 1e10;
   vx_max=-1e10;
   vy_min= 1e10;
   vy_max=-1e10;
   vz_min= 1e10;
   vz_max=-1e10;

   // Perform a slab-decomposed read
   size_t   n_halos_allocate;
   size_t   n_halos_local;
   int      i_halo; 
   char    *line=NULL; 
   size_t   line_length=0;
   GBPREAL  x_in;
   GBPREAL  y_in;
   GBPREAL  z_in;
   GBPREAL  vx_in;
   GBPREAL  vy_in;
   GBPREAL  vz_in;
   GBPREAL *x_halos;
   GBPREAL *y_halos;
   GBPREAL *z_halos;
   GBPREAL *vx_halos;
   GBPREAL *vy_halos;
   GBPREAL *vz_halos;
   int      PHK_min_local;
   int      PHK_max_local;
   size_t  *PHK_halo;
   size_t  *read_index;
   if(check_mode_for_flag(mode,READ_GROUPING_SLAB)){
      // Count the number of halos that will be read (slab decomposed)
      for(i_halo=0,n_halos_allocate=0;i_halo<n_halos;i_halo++){
         grab_next_line_data(fp_in,&line,&line_length);
         grab_real(line,x_column,&x_in);
         if(flag_add_zspace_x){
            grab_real(line,vx_column,&vx_in);
            x_in+=(GBPREAL)(1e3*h_Hubble*((double)vx_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
         }
         force_periodic(&x_in,0.,(GBPREAL)box_size);
         if(x_in>=slab->x_min_local && x_in<slab->x_max_local)
            n_halos_allocate++;
      }
      rewind(fp_in);

      // Allocate storage arrays
      x_halos     =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
      y_halos     =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
      z_halos     =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
      vx_halos=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
      vy_halos=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
      vz_halos=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);

      // Read the halos one at a time, making decisions about where they will go
      for(i_halo=0,n_halos_local=0;i_halo<n_halos;i_halo++){
        grab_next_line_data(fp_in,&line,&line_length);
        grab_real(line,x_column, &x_in);
        if(flag_add_zspace_x){
           grab_real(line,vx_column,&vx_in);
           x_in+=(GBPREAL)(1e3*h_Hubble*((double)vx_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
        }
        force_periodic(&x_in,0.,(GBPREAL)box_size);
        if(x_in>=slab->x_min_local && x_in<slab->x_max_local){
           grab_real(line,y_column,&y_in);
           grab_real(line,z_column,&z_in);
           if(flag_add_zspace_y){
              grab_real(line,vy_column,&vy_in);
              y_in+=(GBPREAL)(1e3*h_Hubble*((double)vy_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
           }
           force_periodic(&y_in,0.,(GBPREAL)box_size);
           if(flag_add_zspace_z){
              grab_real(line,vz_column,&vz_in);
              z_in+=(GBPREAL)(1e3*h_Hubble*((double)vz_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
           }
           force_periodic(&z_in,0.,(GBPREAL)box_size);
           x_halos[n_halos_local] =(GBPREAL)x_in;
           y_halos[n_halos_local] =(GBPREAL)y_in;
           z_halos[n_halos_local] =(GBPREAL)z_in;
           grab_real(line,vx_sub_column,&vx_in);
           grab_real(line,vy_sys_column,&vy_in);
           grab_real(line,vz_sys_column,&vz_in);
           vx_halos[n_halos_local]=(GBPREAL)vx_in;
           vy_halos[n_halos_local]=(GBPREAL)vy_in;
           vz_halos[n_halos_local]=(GBPREAL)vz_in;

           // Keep track of ranges for a sanity check
           x_min =MIN(x_min, x_halos[n_halos_local]);
           x_max =MAX(x_max, x_halos[n_halos_local]);
           y_min =MIN(y_min, y_halos[n_halos_local]);
           y_max =MAX(y_max, y_halos[n_halos_local]);
           z_min =MIN(z_min, z_halos[n_halos_local]);
           z_max =MAX(z_max, z_halos[n_halos_local]);
           vx_min=MIN(vx_min,vx_halos[n_halos_local]);
           vx_max=MAX(vx_max,vx_halos[n_halos_local]);
           vy_min=MIN(vy_min,vy_halos[n_halos_local]);
           vy_max=MAX(vy_max,vy_halos[n_halos_local]);
           vz_min=MIN(vz_min,vz_halos[n_halos_local]);
           vz_max=MAX(vz_max,vz_halos[n_halos_local]);

           // Count the number of halos read locally
           n_halos_local++;
        }
      }
   }
   else if(check_mode_for_flag(mode,READ_GROUPING_PHK)){
      // Determine what key size to use.  We want the size of each key's region to exceed
      //   the largest scale we're interested in.
      SID_log("Using %d bit-per-dimension keys (%d keys; size=%le [Mpc/h]).",SID_LOG_COMMENT,
              n_bits_PHK,(int)PHK_N_KEYS_3D(n_bits_PHK),box_size/pow(2.,(double)(n_bits_PHK)));

      // If this is being done in parallel, we need to determine the key range of each core.
      if(SID.n_proc>1){
         // Count the number of halos that will be 
         //    read for key generation & sorting
         int     i_halo;
         for(i_halo=0,n_halos_allocate=0;i_halo<n_halos;i_halo++){
            if(i_halo%SID.n_proc==SID.My_rank) n_halos_allocate++; // Just for setting PHKs
         }

         // Allocate RAM for PHK keys
         PHK_halo=(size_t *)SID_malloc(sizeof(size_t)*n_halos_allocate);

         // Compute keys for the objects we're reading
         int     j_halo;
         GBPREAL z_space_i      =0;
         GBPREAL z_space        =0;
         GBPREAL z_space_local  =0;
         size_t  n_z_space      =0;
         size_t  n_z_space_local=0;
         SID_log("Generating PHKs for domain decomposition...",SID_LOG_OPEN|SID_LOG_TIMER);
         for(i_halo=0,j_halo=0;i_halo<n_halos;i_halo++){
            grab_next_line_data(fp_in,&line,&line_length);
            if(i_halo%SID.n_proc==SID.My_rank){
               // Read x,y,z
               grab_real(line,x_column,&x_in);
               grab_real(line,y_column,&y_in);
               grab_real(line,z_column,&z_in);
               // Apply z-space distortions if needed
               if(flag_add_zspace_x){
                  grab_real(line,vx_column,&vx_in);
                  z_space_i=(GBPREAL)(1e3*h_Hubble*((double)vx_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
                  x_in+=z_space_i;
                  z_space_local+=z_space_i;
                  n_z_space_local++;
               }
               else if(flag_add_zspace_y){
                  grab_real(line,vy_column,&vy_in);
                  z_space_i=(GBPREAL)(1e3*h_Hubble*((double)vy_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
                  y_in+=z_space_i;
                  z_space_local+=z_space_i;
                  n_z_space_local++;
               }
               else if(flag_add_zspace_z){
                  grab_real(line,vz_column,&vz_in);
                  z_space_i=(GBPREAL)(1e3*h_Hubble*((double)vz_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
                  z_in+=z_space_i;
                  z_space_local+=z_space_i;
                  n_z_space_local++;
               }
               force_periodic(&x_in,0.,box_size);
               force_periodic(&y_in,0.,box_size);
               force_periodic(&z_in,0.,box_size);
               // Compute key
               PHK_halo[j_halo++]=(size_t)compute_PHK_from_Cartesian(n_bits_PHK,3,(double)x_in/box_size,(double)y_in/box_size,(double)z_in/box_size);
            }
         }
         SID_Allreduce(&n_z_space_local,&n_z_space,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
         if(n_z_space>0){
            SID_Allreduce(&z_space_local,&z_space,1,SID_DOUBLE,SID_SUM,SID.COMM_WORLD);
            z_space/=(double)n_z_space;
            SID_log("Average z-space dispacement=%le [Mpc/h]",SID_LOG_COMMENT,z_space);
         }
         rewind(fp_in);
         SID_log("Done.",SID_LOG_CLOSE);

         // Sanity check
         if(j_halo!=n_halos_allocate)
            SID_trap_error("The full allocation of halos was not read (ie. %zd!=%zd).",ERROR_LOGIC,j_halo,n_halos_allocate);

         // Perform a global sort of the keys.  PHK_halo_rank will be the
         //   index *globally* of a given key.  We will use this index to
         //   decide on the range of keys for each rank.
         size_t *PHK_halo_index=NULL;
         size_t *PHK_halo_rank =NULL;
         SID_log("Generating global ranking of PHKs...",SID_LOG_OPEN|SID_LOG_TIMER);
         SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
         sort(PHK_halo,      (size_t)n_halos_allocate,&PHK_halo_index,SID_SIZE_T,SORT_GLOBAL,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
         sort(PHK_halo_index,(size_t)n_halos_allocate,&PHK_halo_rank, SID_SIZE_T,SORT_GLOBAL,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
         SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
         SID_log("Done.",SID_LOG_CLOSE);

         // Perfrom sort of local PHKs
         SID_free(SID_FARG PHK_halo_index);
         PHK_halo_index=NULL;
         merge_sort(PHK_halo,(size_t)n_halos_allocate,&PHK_halo_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);

         // Read an old decomposition if one exists
         if(ADaPS_exist((plist->data),"PHK_min_local_halos")){
            SID_log("Using the last domain decomposition...",SID_LOG_OPEN|SID_LOG_TIMER);
            PHK_min_local=((int *)ADaPS_fetch((plist->data),"PHK_min_local_halos"))[0];
            PHK_max_local=((int *)ADaPS_fetch((plist->data),"PHK_max_local_halos"))[0];
            // Count how many halos will end-up on each rank
            int i_rank;
            for(i_rank=0;i_rank<SID.n_proc;i_rank++){
               size_t n_halos_bcast;
               int    PHK_min_bcast;
               int    PHK_max_bcast;
               PHK_min_bcast=PHK_min_local;
               PHK_max_bcast=PHK_max_local;
               SID_Bcast(&PHK_min_bcast,sizeof(int),i_rank,SID.COMM_WORLD);
               SID_Bcast(&PHK_max_bcast,sizeof(int),i_rank,SID.COMM_WORLD);
               for(i_halo=0,n_halos_bcast=0;i_halo<n_halos_allocate;i_halo++){
                  if(PHK_halo[i_halo]>=(size_t)PHK_min_bcast && PHK_halo[i_halo]<=(size_t)PHK_max_bcast) n_halos_bcast++;
               }
               SID_Allreduce(SID_IN_PLACE,&n_halos_bcast,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
               if(SID.My_rank==i_rank)
                  n_halos_local=n_halos_bcast;
            }
            SID_log("Done.",SID_LOG_CLOSE);
         }
         // Create a new decomposition
         else{
            // Decide on the key ranges for each core
            SID_log("Performing domain decomposition...",SID_LOG_OPEN|SID_LOG_TIMER);
            int    i_rank;
            int    PHK_start;
            size_t i_halo_start;
            size_t i_halo_stop;
            size_t i_halo_target;
            size_t n_halos_remaining;
            int    n_ranks_remaining;
            PHK_start        =0;
            i_halo_start     =0;
            n_halos_remaining=n_halos;
            for(i_rank=0,n_ranks_remaining=SID.n_proc;i_rank<SID.n_proc;i_rank++,n_ranks_remaining--){
               // If all halos have already been allocated, give a dummy range of keys for this rank ...
               if(n_halos_remaining==0){
                  if(i_rank==0 && i_rank==SID.My_rank){
                     PHK_min_local=(int)0;
                     PHK_max_local=(int)PHK_N_KEYS_3D(n_bits_PHK)-1;
                  }
                  else if(i_rank==SID.My_rank){
                     PHK_min_local=(int)PHK_N_KEYS_3D(n_bits_PHK)-1+i_rank;
                     PHK_max_local=(int)PHK_N_KEYS_3D(n_bits_PHK)-1+i_rank;
                  }
                  n_halos_local=0;
               }
               // ... else we need to continue allocating.
               else{
                  // Set the target halo.  It's PHK will set the upper range of PHKs for this rank.
                  size_t i_halo_step;
                  i_halo_step  =n_halos_remaining/n_ranks_remaining;
                  i_halo_step  =MAX(i_halo_step,1);
                  i_halo_target=i_halo_start+i_halo_step-1;
                  i_halo_target=MIN(i_halo_target,n_halos-1);

                  // Find the rank that has the target rank and it's PHK
                  int     flag_found_target=FALSE;
                  int     PHK_Bcast;
                  int     flag_check;
                  PHK_Bcast =0; 
                  for(i_halo=0;i_halo<n_halos_allocate && !flag_found_target;i_halo++){
                     if(PHK_halo_rank[i_halo]==i_halo_target){
                        flag_found_target=TRUE;
                        PHK_Bcast        =(int)PHK_halo[i_halo];
                        break;
                     }
                  }
                  SID_Allreduce(SID_IN_PLACE,&PHK_Bcast,1,SID_INT,SID_MAX,SID.COMM_WORLD);

                  // Sanity check
                  SID_Allreduce(&flag_found_target,&flag_check,1,SID_INT,SID_SUM,SID.COMM_WORLD);
                  if(flag_check!=1)
                     SID_trap_error("The desired halo (%zd) was found by %d processes during the PHK domain decomposition.",
                                    ERROR_LOGIC,i_halo_target,flag_check);

                  // Determine the global index of the last halo that has the same key as the target.
                  size_t i_halo_highest=0;
                  i_halo=0;
                  if(n_halos_allocate>1){
                     if(PHK_halo[PHK_halo_index[i_halo]]<=PHK_Bcast){
                        while(PHK_halo[PHK_halo_index[i_halo+1]]<=PHK_Bcast){
                           i_halo++;
                           if(i_halo>=(n_halos_allocate-1)) break;
                        }
                        i_halo_highest=PHK_halo_rank[PHK_halo_index[i_halo]];
                     }
                  }
                  SID_Allreduce(&i_halo_highest,&i_halo_stop,1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);

                  // Set the starting halo index of the next rank and 
                  //   the number of halos on this rank and remaining to be allocated.
                  n_halos_remaining-=(i_halo_stop-i_halo_start+1);

                  // Set the i_rank'th's PHK range
                  if(i_rank==SID.My_rank){
                     n_halos_local=i_halo_stop-i_halo_start+1;
                     PHK_min_local=PHK_start;
                     if(n_halos_remaining==0)
                        PHK_max_local=PHK_N_KEYS_3D(n_bits_PHK)-1;
                     else
                        PHK_max_local=PHK_Bcast;
                  }
                  i_halo_start=i_halo_stop+1;
                  PHK_start   =PHK_Bcast+1; // This will be the starting PHK of the next rank
               } // if halos remaining
            } // i_rank
            SID_log("Done.",SID_LOG_CLOSE);

            // Store some stuff
            ADaPS_store(&(plist->data),&PHK_min_local,"PHK_min_local_halos",ADaPS_SCALAR_INT);
            ADaPS_store(&(plist->data),&PHK_max_local,"PHK_max_local_halos",ADaPS_SCALAR_INT);
         }

         // Change n_halos_allocate for the final read
         n_halos_allocate=n_halos_local;
      
         // Clean-up
         SID_free(SID_FARG PHK_halo);
         SID_free(SID_FARG PHK_halo_rank);
         SID_free(SID_FARG PHK_halo_index);
      }
      // ... else, if there is just one core, give it the whole range of keys.
      else{
         n_halos_allocate=n_halos;
         PHK_min_local   =0;
         PHK_max_local   =PHK_N_KEYS_3D(n_bits_PHK)-1;
      }

      // Ok, each core now knows it's key range and the number of halos it's going to read. Allocate RAM ...
      x_halos   =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
      y_halos   =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
      z_halos   =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
      vx_halos  =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
      vy_halos  =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
      vz_halos  =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_halos_allocate);
      read_index=(size_t  *)SID_malloc(sizeof(size_t) *n_halos_allocate);
      PHK_halo  =(size_t  *)SID_malloc(sizeof(size_t) *n_halos_allocate);

      // Determine which PHKs are on the boundary of this rank's domain
      PHK_t *keys_boundary  =NULL;
      int    n_keys_boundary=0;
      compute_PHK_boundary_keys(n_bits_PHK,PHK_min_local,PHK_max_local,PHK_width,&n_keys_boundary,&keys_boundary);

      // ... read the file once to determine boundary/interior member counts ...
      size_t n_boundary=0;
      size_t n_interior=0;
      if(SID.n_proc>1)
         SID_log("Performing boundary/interior counts...",SID_LOG_OPEN|SID_LOG_TIMER);
      for(i_halo=0,n_halos_local=0;i_halo<n_halos;i_halo++){
         grab_next_line_data(fp_in,&line,&line_length);
         grab_real(line,x_column,&x_in);
         grab_real(line,y_column,&y_in);
         grab_real(line,z_column,&z_in);
         if(flag_add_zspace_x){
            grab_real(line,vx_column,&vx_in);
            x_in+=(GBPREAL)(1e3*h_Hubble*((double)vx_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
         }
         else if(flag_add_zspace_y){
            grab_real(line,vy_column,&vy_in);
            y_in+=(GBPREAL)(1e3*h_Hubble*((double)vy_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
         }
         else if(flag_add_zspace_z){
            grab_real(line,vz_column,&vz_in);
            z_in+=(GBPREAL)(1e3*h_Hubble*((double)vz_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
         }
         force_periodic(&x_in,0.,box_size);
         force_periodic(&y_in,0.,box_size);
         force_periodic(&z_in,0.,box_size);

         // Compute key
         PHK_t PHK_i;
         PHK_i=compute_PHK_from_Cartesian(n_bits_PHK,3,(double)x_in/box_size,(double)y_in/box_size,(double)z_in/box_size);

         // Assign to the right rank
         if(PHK_i>=(PHK_t)PHK_min_local && PHK_i<=(PHK_t)PHK_max_local){
            if(n_halos_local>=n_halos_allocate)
               SID_trap_error("Trying to read past storage allocation in PHK doain decomposition.",ERROR_LOGIC);
            if(is_a_member(&PHK_i,keys_boundary,n_keys_boundary,SID_PHK_T))
               n_boundary++;
            else
               n_interior++;
            n_halos_local++;
         }
      }
      rewind(fp_in);

      // Sanity check
      if(n_halos_local!=n_halos_allocate)
         SID_trap_error("The number of halos counted does not correspond to what was allocated (ie. %zd!=%zd).",ERROR_LOGIC,n_halos_local,n_halos_allocate);

      // Report decomposition results
      if(SID.n_proc>1){
         int i_rank;
         SID_log("Results of domain decomposition:",SID_LOG_OPEN);
         for(i_rank=0;i_rank<SID.n_proc;i_rank++){
            size_t n_halos_report;
            size_t n_boundary_report;
            int    PHK_min_report;
            int    PHK_max_report;
            n_halos_report   =n_halos_allocate;
            n_boundary_report=n_boundary;
            PHK_min_report   =PHK_min_local;
            PHK_max_report   =PHK_max_local;
            SID_Bcast(&n_halos_report,   sizeof(size_t),i_rank,SID.COMM_WORLD);
            SID_Bcast(&n_boundary_report,sizeof(size_t),i_rank,SID.COMM_WORLD);
            SID_Bcast(&PHK_min_report,   sizeof(int),   i_rank,SID.COMM_WORLD);
            SID_Bcast(&PHK_max_report,   sizeof(int),   i_rank,SID.COMM_WORLD);
            SID_log("Rank #%03d: n_objects=%6zd n_boundary=%6zd PHK=%d->%d",SID_LOG_COMMENT,i_rank,n_halos_report,n_boundary_report,PHK_min_report,PHK_max_report);
         }
         SID_log("",SID_LOG_SILENT_CLOSE);
      }

      // ... and finally, read everything into the arrays, placing things on the boundary of the 
      //     domain first and things in the interior at the end.  This makes boundary exchanges
      //     easier later.
      SID_log("Performing final read...",SID_LOG_OPEN|SID_LOG_TIMER);
      size_t i_boundary=0;
      size_t i_interior=n_boundary;
      for(i_halo=0,n_halos_local=0;i_halo<n_halos;i_halo++){
         grab_next_line_data(fp_in,&line,&line_length);
         grab_real(line,x_column,&x_in);
         grab_real(line,y_column,&y_in);
         grab_real(line,z_column,&z_in);
         if(flag_add_zspace_x){
            grab_real(line,vx_column,&vx_in);
            x_in+=(GBPREAL)(1e3*h_Hubble*((double)vx_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
         }
         else if(flag_add_zspace_y){
            grab_real(line,vy_column,&vy_in);
            y_in+=(GBPREAL)(1e3*h_Hubble*((double)vy_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
         }
         else if(flag_add_zspace_z){
            grab_real(line,vz_column,&vz_in);
            z_in+=(GBPREAL)(1e3*h_Hubble*((double)vz_in)/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
         }
         force_periodic(&x_in,0.,box_size);
         force_periodic(&y_in,0.,box_size);
         force_periodic(&z_in,0.,box_size);

         // Compute key
         PHK_t PHK_i;
         PHK_i=compute_PHK_from_Cartesian(n_bits_PHK,3,(double)x_in/box_size,(double)y_in/box_size,(double)z_in/box_size);

         // Assign to the right rank
         if(PHK_i>=(PHK_t)PHK_min_local && PHK_i<=(PHK_t)PHK_max_local){
            size_t i_store;
            if(n_halos_local>=n_halos_allocate)
               SID_trap_error("Trying to read past storage allocation in PHK doain decomposition.",ERROR_LOGIC);
            if(is_a_member(&PHK_i,keys_boundary,n_keys_boundary,SID_PHK_T))
               i_store=i_boundary++;
            else
               i_store=i_interior++;
            grab_real(line,vx_column,&vx_in);
            grab_real(line,vy_column,&vy_in);
            grab_real(line,vz_column,&vz_in);
            x_halos[i_store]   =x_in;
            y_halos[i_store]   =y_in;
            z_halos[i_store]   =z_in;
            vx_halos[i_store]  =vx_in;
            vy_halos[i_store]  =vy_in;
            vz_halos[i_store]  =vz_in;
            read_index[i_store]=i_halo;
            PHK_halo[i_store]  =(size_t)PHK_i;

            // Keep track of ranges for a sanity check
            x_min =MIN(x_min, x_halos[i_store]);
            x_max =MAX(x_max, x_halos[i_store]);
            y_min =MIN(y_min, y_halos[i_store]);
            y_max =MAX(y_max, y_halos[i_store]);
            z_min =MIN(z_min, z_halos[i_store]);
            z_max =MAX(z_max, z_halos[i_store]);
            vx_min=MIN(vx_min,vx_halos[i_store]);
            vx_max=MAX(vx_max,vx_halos[i_store]);
            vy_min=MIN(vy_min,vy_halos[i_store]);
            vy_max=MAX(vy_max,vy_halos[i_store]);
            vz_min=MIN(vz_min,vz_halos[i_store]);
            vz_max=MAX(vz_max,vz_halos[i_store]);

            // Count the number of halos read locally
            n_halos_local++;
         }
      }
      SID_log("Done.",SID_LOG_CLOSE);

      // Sanity check
      if(n_halos_local!=n_halos_allocate)
         SID_trap_error("The number of halos read does not correspond to what was allocated (ie. %zd!=%zd).",ERROR_LOGIC,n_halos_local,n_halos_allocate);

      // Sort the local PHKs
      size_t *PHK_index;
      size_t *PHK_index_boundary;
      merge_sort(PHK_halo,n_halos_allocate,&PHK_index,         SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
      merge_sort(PHK_halo,n_boundary,      &PHK_index_boundary,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);

      // Store read indices, PHKs and their sort indices
      ADaPS_store(&(plist->data),(void *)read_index,        "read_index_halos",        ADaPS_DEFAULT);
      ADaPS_store(&(plist->data),(void *)PHK_halo,          "PHK_halos",               ADaPS_DEFAULT);
      ADaPS_store(&(plist->data),(void *)PHK_index,         "PHK_index_halos",         ADaPS_DEFAULT);
      ADaPS_store(&(plist->data),(void *)PHK_index_boundary,"PHK_index_boundary_halos",ADaPS_DEFAULT);
      ADaPS_store(&(plist->data),(void *)(&PHK_min_local),  "PHK_min_local_halos",     ADaPS_SCALAR_INT);
      ADaPS_store(&(plist->data),(void *)(&PHK_max_local),  "PHK_max_local_halos",     ADaPS_SCALAR_INT);
      ADaPS_store(&(plist->data),(void *)(&n_boundary),     "n_boundary_halos",        ADaPS_SCALAR_SIZE_T);

      // Clean-up
      SID_free(SID_FARG keys_boundary);

      if(SID.n_proc>1)
         SID_log("Done.",SID_LOG_CLOSE);

   } // PHK decomposition
   fclose(fp_in);

   // Some sanity checks on what we have read
   size_t n_halos_read;
   SID_Allreduce(&n_halos_local,&n_halos_read,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
   if(n_halos_read!=n_halos)
      SID_trap_error("The correct number of halos was not read (ie. %lld!=%lld).  There must be a coordinate/box size problem.",
                     ERROR_LOGIC,n_halos_read,n_halos);
   SID_Allreduce(SID_IN_PLACE,&x_min, 1,SID_REAL,SID_MIN,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&x_max, 1,SID_REAL,SID_MAX,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&y_min, 1,SID_REAL,SID_MIN,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&y_max, 1,SID_REAL,SID_MAX,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&z_min, 1,SID_REAL,SID_MIN,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&z_max, 1,SID_REAL,SID_MAX,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&vx_min,1,SID_REAL,SID_MIN,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&vx_max,1,SID_REAL,SID_MAX,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&vy_min,1,SID_REAL,SID_MIN,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&vy_max,1,SID_REAL,SID_MAX,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&vz_min,1,SID_REAL,SID_MIN,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,&vz_max,1,SID_REAL,SID_MAX,SID.COMM_WORLD);
   SID_log("x_range =%le->%le",SID_LOG_COMMENT,x_min, x_max);
   SID_log("y_range =%le->%le",SID_LOG_COMMENT,y_min, y_max);
   SID_log("z_range =%le->%le",SID_LOG_COMMENT,z_min, z_max);
   SID_log("vx_range=%le->%le",SID_LOG_COMMENT,vx_min,vx_max);
   SID_log("vy_range=%le->%le",SID_LOG_COMMENT,vy_min,vy_max);
   SID_log("vz_range=%le->%le",SID_LOG_COMMENT,vz_min,vz_max);

   // Store Arrays
   ADaPS_store(&(plist->data),(void *)&n_halos,      "n_all_halos", ADaPS_SCALAR_SIZE_T);
   ADaPS_store(&(plist->data),(void *)&n_halos_local,"n_halos",     ADaPS_SCALAR_SIZE_T);
   ADaPS_store(&(plist->data),(void *)x_halos,       "x_halos",     ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)y_halos,       "y_halos",     ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)z_halos,       "z_halos",     ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)vx_halos,      "vx_halos",    ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)vy_halos,      "vy_halos",    ADaPS_DEFAULT);
   ADaPS_store(&(plist->data),(void *)vz_halos,      "vz_halos",    ADaPS_DEFAULT);

   // Clean-up
   SID_free(SID_FARG line);
   va_end(vargs);
   SID_log("Done.",SID_LOG_CLOSE);
}

