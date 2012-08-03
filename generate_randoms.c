#include <gbpLib.h>
#include <gbpClustering.h>

void generate_randoms(cfunc_info *cfunc,plist_info *plist,char *species_name,char *random_name){

   // Fetch the number of objects and set the number of random
   size_t n_random_target;
   size_t n_species;
   size_t n_species_local;
   n_species      =((size_t *)ADaPS_fetch(plist->data,"n_all_%s",species_name))[0];
   n_species_local=((size_t *)ADaPS_fetch(plist->data,"n_%s",    species_name))[0];
   n_random_target=(size_t)cfunc->F_random*n_species;
   SID_log("Generating %d random halos...",SID_LOG_OPEN,n_random_target);

   // We are assuming PHK decomposition has already been performed
   int PHK_min_local;
   int PHK_max_local;
   PHK_min_local=((int *)ADaPS_fetch(plist->data,"PHK_min_local_%s",species_name))[0];
   PHK_max_local=((int *)ADaPS_fetch(plist->data,"PHK_max_local_%s",species_name))[0];

   // Get the local domain's boundary keys
   int    n_bits_PHK;
   int    n_keys_boundary;
   PHK_t *keys_boundary;
   compute_PHK_boundary_keys(cfunc->n_bits_PHK,PHK_min_local,PHK_max_local,cfunc->PHK_width,&n_keys_boundary,&keys_boundary);

   // Initialize the random number generator
   RNG_info  RNG;
   int       seed=1327621;
   init_RNG(&seed,&RNG,RNG_GLOBAL);

   // Generate the random and determine how many belong to this rank's boundary/interior
   GBPREAL x_i;
   GBPREAL y_i;
   GBPREAL z_i;
   PHK_t   PHK_i;
   size_t  i_random;
   size_t  n_random       =0;
   size_t  n_random_local=0;
   size_t  n_boundary     =0;
   for(i_random=0;i_random<n_random_target;i_random++){
      x_i=random_number(&RNG);force_periodic(&x_i,0.,1.);
      y_i=random_number(&RNG);force_periodic(&y_i,0.,1.);
      z_i=random_number(&RNG);force_periodic(&z_i,0.,1.);
      PHK_i=compute_PHK_from_Cartesian(cfunc->n_bits_PHK,3,(double)x_i,(double)y_i,(double)z_i);
      if(PHK_i>=(PHK_t)PHK_min_local && PHK_i<=(PHK_t)PHK_max_local){
         n_random_local++;
         if(is_a_member(&PHK_i,keys_boundary,n_keys_boundary,SID_PHK_T))
            n_boundary++;
      }
   }
   SID_Allreduce(&n_random_local,&n_random,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);

   // Sanity check
   if(n_random!=n_random_target)
      SID_trap_error("There's a count mismatch while generating random (ie. %zd!=%zd).",ERROR_LOGIC,n_random,n_random_target);

   // Create the arrays
   GBPREAL *x_random;
   GBPREAL *y_random;
   GBPREAL *z_random;
   size_t  *PHK_random;
   x_random  =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_random_local);
   y_random  =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_random_local);
   z_random  =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_random_local);
   PHK_random=(size_t  *)SID_malloc(sizeof(size_t) *n_random_local);

   // Reinitialize the RNG so we generate the same stream as above
   init_RNG(&seed,&RNG,RNG_GLOBAL);

   // Generate the random again and store them in the arrays
   size_t j_random=0;
   size_t i_store;
   size_t i_boundary=0;
   size_t i_interior=n_boundary;
   for(i_random=0;i_random<n_random_target;i_random++){
      x_i  =random_number(&RNG);force_periodic(&x_i,0.,1.);
      y_i  =random_number(&RNG);force_periodic(&y_i,0.,1.);
      z_i  =random_number(&RNG);force_periodic(&z_i,0.,1.);
      PHK_i=compute_PHK_from_Cartesian(cfunc->n_bits_PHK,3,(double)x_i,(double)y_i,(double)z_i);
      if(PHK_i>=(PHK_t)PHK_min_local && PHK_i<=(PHK_t)PHK_max_local){
         int i_store;
         if(is_a_member(&PHK_i,keys_boundary,n_keys_boundary,SID_PHK_T))
            i_store=i_boundary++;
         else
            i_store=i_interior++;
         x_random[i_store]  =x_i*cfunc->box_size;
         y_random[i_store]  =y_i*cfunc->box_size;
         z_random[i_store]  =z_i*cfunc->box_size;
         PHK_random[i_store]=PHK_i;
      }
   }

   // Report decomposition results
   size_t n_species_report;
   size_t n_random_report;
   size_t n_boundary_report;
   if(SID.n_proc>1){
      int i_rank;
      SID_log("Results of domain decomposition:",SID_LOG_OPEN);
      for(i_rank=0;i_rank<SID.n_proc;i_rank++){
         n_species_report =n_species_local;
         n_random_report  =n_boundary;
         n_boundary_report=n_random_local;
         SID_Bcast(&n_species_report, sizeof(size_t),i_rank,SID.COMM_WORLD);
         SID_Bcast(&n_random_report,  sizeof(size_t),i_rank,SID.COMM_WORLD);
         SID_Bcast(&n_boundary_report,sizeof(size_t),i_rank,SID.COMM_WORLD);
         SID_log("Rank #%03d: n_species=%6zd n_random=%6zd n_boundary=%zd",SID_LOG_COMMENT,i_rank,n_species_report,n_random_report,n_boundary_report);
      }
      SID_log("",SID_LOG_SILENT_CLOSE);
   }

   // Tell cfunc to compute the RR's next time.  This only needs to be done once.
   cfunc->flag_compute_RR=TRUE;

   // Sort PHKs
   size_t *PHK_index;
   size_t *PHK_index_boundary;
   merge_sort(PHK_random,n_random_local,&PHK_index,         SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
   merge_sort(PHK_random,n_boundary,    &PHK_index_boundary,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);

   // Store the arrays
   ADaPS_store(&(plist->data),(void *)&n_random,         "n_all_%s",             ADaPS_SCALAR_SIZE_T,random_name);
   ADaPS_store(&(plist->data),(void *)&n_random_local,   "n_%s",                 ADaPS_SCALAR_SIZE_T,random_name);
   ADaPS_store(&(plist->data),(void *)&n_boundary,       "n_boundary_%s",        ADaPS_SCALAR_SIZE_T,random_name);
   ADaPS_store(&(plist->data),(void *)x_random,          "x_%s",                 ADaPS_DEFAULT,      random_name);
   ADaPS_store(&(plist->data),(void *)y_random,          "y_%s",                 ADaPS_DEFAULT,      random_name);
   ADaPS_store(&(plist->data),(void *)z_random,          "z_%s",                 ADaPS_DEFAULT,      random_name);
   ADaPS_store(&(plist->data),(void *)PHK_random,        "PHK_%s",               ADaPS_DEFAULT,      random_name);
   ADaPS_store(&(plist->data),(void *)PHK_index,         "PHK_index_%s",         ADaPS_DEFAULT,      random_name);
   ADaPS_store(&(plist->data),(void *)PHK_index_boundary,"PHK_index_boundary_%s",ADaPS_DEFAULT,      random_name);

   SID_log("Done.",SID_LOG_CLOSE);
}

