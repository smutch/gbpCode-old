#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

typedef struct select_group_params_local select_group_params_local;
struct select_group_params_local{
   int     selection_index;
   int     flag_long_ids_store;
   size_t  n_ids_list;
   void   *ids_list;
   char   *val_list;
};

int count_group_ids_local(int                i_group,
                          int                j_subgroup,
                          int                i_subgroup,
                          int                flag_long_ids,
                          process_halo_info *group_i,
                          process_halo_info *subgroup_i,
                          void              *params);
int count_group_ids_local(int                i_group,
                          int                j_subgroup,
                          int                i_subgroup,
                          int                flag_long_ids,
                          process_halo_info *group_i,
                          process_halo_info *subgroup_i,
                          void              *params){
   ((select_group_params_local *)params)->n_ids_list+=group_i->n_particles;
}

int count_subgroup_ids_local(int                i_group,
                             int                j_subgroup,
                             int                i_subgroup,
                             int                flag_long_ids,
                             process_halo_info *group_i,
                             process_halo_info *subgroup_i,
                             void              *params);
int count_subgroup_ids_local(int                i_group,
                             int                j_subgroup,
                             int                i_subgroup,
                             int                flag_long_ids,
                             process_halo_info *group_i,
                             process_halo_info *subgroup_i,
                             void              *params){
   ((select_group_params_local *)params)->n_ids_list+=subgroup_i->n_particles;
}

int add_group_to_ids_list_local(int                i_group,
                                int                j_subgroup,
                                int                i_subgroup,
                                int                flag_long_ids,
                                process_halo_info *group_i,
                                process_halo_info *subgroup_i,
                                void              *params_in);
int add_group_to_ids_list_local(int                i_group,
                                int                j_subgroup,
                                int                i_subgroup,
                                int                flag_long_ids_in,
                                process_halo_info *group_i,
                                process_halo_info *subgroup_i,
                                void              *params_in){
   select_group_params_local *params=(select_group_params_local *)params_in;
   int     flag_long_ids_store=params->flag_long_ids_store;
   size_t *ids_l=(size_t *)params->ids_list;
   int    *ids_i=(int    *)params->ids_list;
   if(flag_long_ids_in){
      if(flag_long_ids_in!=flag_long_ids_store){
         for(int i_particle=0;i_particle<group_i->n_particles;i_particle++)
            ids_i[params->n_ids_list+i_particle]=(int)(((size_t *)group_i->ids)[i_particle]);
      }
      else
         memcpy(&(ids_l[params->n_ids_list]),group_i->ids,group_i->n_particles*sizeof(size_t));
   }
   else{
      if(flag_long_ids_in!=flag_long_ids_store){
         for(int i_particle=0;i_particle<group_i->n_particles;i_particle++)
            ids_l[params->n_ids_list+i_particle]=(size_t)(((int *)group_i->ids)[i_particle]);
      }
      else
         memcpy(&(ids_i[params->n_ids_list]),group_i->ids,group_i->n_particles*sizeof(int));
   }
   params->n_ids_list+=group_i->n_particles;
}

int add_subgroup_to_ids_list_local(int                i_group,
                                   int                j_subgroup,
                                   int                i_subgroup,
                                   int                flag_long_ids,
                                   process_halo_info *group_i,
                                   process_halo_info *subgroup_i,
                                   void              *params_in);
int add_subgroup_to_ids_list_local(int                i_group,
                                   int                j_subgroup,
                                   int                i_subgroup,
                                   int                flag_long_ids_in,
                                   process_halo_info *group_i,
                                   process_halo_info *subgroup_i,
                                   void              *params_in){
   select_group_params_local *params=(select_group_params_local *)params_in;
   int     flag_long_ids_store=params->flag_long_ids_store;
   size_t *ids_l=(size_t *)params->ids_list;
   int    *ids_i=(int    *)params->ids_list;
   if(flag_long_ids_in){
      if(flag_long_ids_in!=flag_long_ids_store){
         for(int i_particle=0;i_particle<subgroup_i->n_particles;i_particle++)
            ids_i[params->n_ids_list+i_particle]=(int)(((size_t *)subgroup_i->ids)[i_particle]);
      }
      else
         memcpy(&(ids_l[params->n_ids_list]),subgroup_i->ids,subgroup_i->n_particles*sizeof(size_t));
   }
   else{
      if(flag_long_ids_in!=flag_long_ids_store){
         for(int i_particle=0;i_particle<subgroup_i->n_particles;i_particle++)
            ids_l[params->n_ids_list+i_particle]=(size_t)(((int *)subgroup_i->ids)[i_particle]);
      }
      else
         memcpy(&(ids_i[params->n_ids_list]),subgroup_i->ids,subgroup_i->n_particles*sizeof(int));
   }
   params->n_ids_list+=subgroup_i->n_particles;
}

int select_group_index_local(int                i_group,
                             int                j_subgroup,
                             int                i_subgroup,
                             int                flag_long_ids,
                             process_halo_info *group_i,
                             process_halo_info *subgroup_i,
                             void              *params);
int select_group_index_local(int                i_group,
                             int                j_subgroup,
                             int                i_subgroup,
                             int                flag_long_ids,
                             process_halo_info *group_i,
                             process_halo_info *subgroup_i,
                             void              *params){
   return(i_group==((select_group_params_local *)params)->selection_index && j_subgroup==0);
}

int select_subgroup_index_local(int                i_group,
                                int                j_subgroup,
                                int                i_subgroup,
                                int                flag_long_ids,
                                process_halo_info *group_i,
                                process_halo_info *subgroup_i,
                                void              *params);
int select_subgroup_index_local(int                i_group,
                                int                j_subgroup,
                                int                i_subgroup,
                                int                flag_long_ids,
                                process_halo_info *group_i,
                                process_halo_info *subgroup_i,
                                void              *params){
   return(i_subgroup==((select_group_params_local *)params)->selection_index);
}

void make_ids_list(render_info *render,
                   int    i_snap,
                   int    select_function(int                i_group,
                                          int                j_subgroup,
                                          int                i_subgroup,
                                          int                flag_long_ids,
                                          process_halo_info *group_i,
                                          process_halo_info *subgroup_i,
                                          void              *params),
                   int    action_function(int                i_group,
                                          int                j_subgroup,
                                          int                i_subgroup,
                                          int                flag_long_ids,
                                          process_halo_info *group_i,
                                          process_halo_info *subgroup_i,
                                          void              *params),
                   void   **ids_list,
                   char   **val_list,
                   size_t   n_ids_list,
                   int      mode,
                   int      selection);
void make_ids_list(render_info *render,
                   int    i_snap,
                   int    select_function(int                i_group,
                                          int                j_subgroup,
                                          int                i_subgroup,
                                          int                flag_long_ids,
                                          process_halo_info *group_i,
                                          process_halo_info *subgroup_i,
                                          void              *params),
                   int    count_function( int                i_group,
                                          int                j_subgroup,
                                          int                i_subgroup,
                                          int                flag_long_ids,
                                          process_halo_info *group_i,
                                          process_halo_info *subgroup_i,
                                          void              *params),
                   int    action_function(int                i_group,
                                          int                j_subgroup,
                                          int                i_subgroup,
                                          int                flag_long_ids,
                                          process_halo_info *group_i,
                                          process_halo_info *subgroup_i,
                                          void              *params),
                   void   **ids_list,
                   char   **val_list,
                   size_t  *n_ids_list,
                   int      mode,
                   int      selection){

   // Initialize some stuff
   int i_read=render->snap_list[i_snap];
   int flag_long_ids=ADaPS_exist(render->plist_list[i_snap]->data,"flag_LONGIDs");
   select_group_params_local params;
   params.flag_long_ids_store=flag_long_ids;
   params.selection_index    =selection;
   params.n_ids_list         =0;

   // Count the number of particles we will have
   process_SSimPL_halos(render->filename_SSimPL_dir,
                        render->filename_halo_type,
                        render->filename_tree_version,
                        i_read,0,
                        select_function,
                        count_function,
                        &params);

   // Allocate the array for the list
   (*n_ids_list)=params.n_ids_list;
   switch(flag_long_ids){
      case TRUE:
         (*ids_list) =SID_malloc(sizeof(size_t)*(*n_ids_list));
         break;
      case FALSE:
         (*ids_list) =SID_malloc(sizeof(int)*(*n_ids_list));
         break;
   }

   // Allocate a value array
   (*val_list)=NULL;

   // Populate the id list
   params.ids_list  =(*ids_list);
   params.val_list  =(*val_list);
   params.n_ids_list=0;
   process_SSimPL_halos(render->filename_SSimPL_dir,
                        render->filename_halo_type,
                        render->filename_tree_version,
                        i_read,1,
                        select_function,
                        action_function,
                        &params);

   // Sort the list in place
   switch(flag_long_ids){
      case TRUE:
         merge_sort((*ids_list),(*n_ids_list),NULL,SID_SIZE_T,SORT_INPLACE_ONLY,TRUE);
         break;
      case FALSE:
         merge_sort((*ids_list),(*n_ids_list),NULL,SID_INT,   SORT_INPLACE_ONLY,TRUE);
         break;
   }
}

void apply_mark_list(render_info   *render,
                     mark_arg_info *arg,
                     int            i_snap,
                     int            i_type,
                     void          *ids_list,
                     char          *value_list,
                     size_t         n_ids_list,
                     char          *mark,
                     size_t        *mark_count);
void apply_mark_list(render_info   *render,
                     mark_arg_info *arg,
                     int            i_snap,
                     int            i_type,
                     void          *ids_list,
                     char          *value_list,
                     size_t         n_ids_list,
                     char          *mark,
                     size_t        *mark_count){

   // Fetch the IDs and make sure we have sort indices for them
   int     flag_long_ids  =ADaPS_exist(render->plist_list[i_snap]->data,"flag_LONGIDs");
   size_t  n_species_local=((size_t *)ADaPS_fetch(render->plist_list[i_snap]->data,"n_%s", render->plist_list[i_snap]->species[i_type]))[0];
   size_t *ids            = (size_t *)ADaPS_fetch(render->plist_list[i_snap]->data,"id_%s",render->plist_list[i_snap]->species[i_type]);
   size_t *ids_index;
   if(!ADaPS_exist(render->plist_list[i_snap]->data,"id_%s_index",render->plist_list[i_snap]->species[i_type])){
      merge_sort(ids,n_species_local,&ids_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
      ADaPS_store(&(render->plist_list[i_snap]->data),ids_index,"id_%s_index",ADaPS_DEFAULT,render->plist_list[i_snap]->species[i_type]);
   }
   else
      ids_index=(size_t *)ADaPS_fetch(render->plist_list[i_snap]->data,"id_%s_index",render->plist_list[i_snap]->species[i_type]);
   // Perform particle marking
   size_t i_particle,j_particle;
   for(i_particle=j_particle=0;i_particle<n_ids_list && j_particle<n_species_local;i_particle++){
      size_t id_i;
      switch(flag_long_ids){
         case TRUE:
            id_i=(size_t)(((size_t *)ids_list)[i_particle]);
            break;
         case FALSE:
            id_i=(size_t)(((int    *)ids_list)[i_particle]);
            break;
      }
      while(id_i>ids[ids_index[j_particle]] && j_particle<(n_species_local-1)) j_particle++;
      if   (id_i>ids[ids_index[j_particle]])                                   j_particle++;
      if(j_particle>=n_species_local) break;
      if(id_i==ids[ids_index[j_particle]]){
         switch(value_list==NULL){
            case TRUE:
               mark[ids_index[j_particle]]=arg->value;
               break;
            case FALSE:
               mark[ids_index[j_particle]]=value_list[i_particle];
               break;
         }
         (*mark_count)++;
      }
   }
}

void execute_marking_argument_local(render_info *render,mark_arg_info *arg);
void execute_marking_argument_local(render_info *render,mark_arg_info *arg){
   if(!strcmp(arg->type,"all"))
      SID_log("Marking %s particles...",SID_LOG_OPEN|SID_LOG_TIMER,arg->species);
   else if(!strcmp(arg->type,"sphere"))
      SID_log("Marking spherical volume...",SID_LOG_OPEN|SID_LOG_TIMER);
   else if(!strcmp(arg->type,"group_index"))
      SID_log("Marking group halo (index=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,arg->ival[0]);
   else if(!strcmp(arg->type,"subgroup_index"))
      SID_log("Marking subgroup halo (index=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,arg->ival[0]);
   else
      SID_trap_error("Invalid selection type {%s} in perform_marking().",ERROR_LOGIC,arg->type);
   for(int i_snap=0;i_snap<render->n_interpolate;i_snap++){
      if(render->n_interpolate>1)
         SID_log("Processing snapshot %d...",SID_LOG_OPEN|SID_LOG_TIMER);
      plist_info *plist=render->plist_list[i_snap];
      size_t mark_count =0;
      size_t n_particles=0;
      for(int i_type=0;i_type<N_GADGET_TYPE;i_type++){
         char species_name[256];
         strcpy(species_name,plist->species[i_type]);
         if(!strcmp(arg->species,"all") || !strcmp(arg->species,species_name)){
            if(ADaPS_exist(plist->data,"n_%s",species_name)){
               // Fetch or create the mark array
               char *mark;
               size_t n_species_local=((size_t *)ADaPS_fetch(plist->data,"n_%s",species_name))[0];
               if(!ADaPS_exist(plist->data,"mark_%s",species_name)){
                  mark=(char *)SID_calloc(sizeof(char)*n_species_local);
                  ADaPS_store(&(plist->data),mark,"mark_%s",ADaPS_DEFAULT,species_name);
               }
               else
                  mark=(char *)ADaPS_fetch(plist->data,"mark_%s",species_name);
               // Perform marking logic here ...
               // ... set all particles to a value.
               if(!strcmp(arg->type,"all")){
                  for(size_t i_particle=0;i_particle<n_species_local;i_particle++)
                     mark[i_particle]=arg->value;
                  mark_count +=n_species_local;
               }
               // ... set all particles in a sphere to a value.
               else if(!strcmp(arg->type,"sphere")){
                  double x_cen,y_cen,z_cen,r2;
                  GBPREAL *x=NULL;
                  GBPREAL *y=NULL;
                  GBPREAL *z=NULL;
                  if(!render->camera->flag_velocity_space){
                     x_cen=arg->dval[0]*M_PER_MPC/render->h_Hubble;
                     y_cen=arg->dval[1]*M_PER_MPC/render->h_Hubble;
                     z_cen=arg->dval[2]*M_PER_MPC/render->h_Hubble;
                     r2   =arg->dval[3]*arg->dval[3]*pow(M_PER_MPC/render->h_Hubble,2.);
                     x=(GBPREAL *)ADaPS_fetch(plist->data,"x_%s",species_name);
                     y=(GBPREAL *)ADaPS_fetch(plist->data,"y_%s",species_name);
                     z=(GBPREAL *)ADaPS_fetch(plist->data,"z_%s",species_name);
                  }
                  else{
                     x_cen=arg->dval[0]*1e3;
                     y_cen=arg->dval[1]*1e3;
                     z_cen=arg->dval[2]*1e3;
                     r2   =arg->dval[3]*arg->dval[3]*pow(1e3,2.);
                     x=(GBPREAL *)ADaPS_fetch(plist->data,"vx_%s",species_name);
                     y=(GBPREAL *)ADaPS_fetch(plist->data,"vy_%s",species_name);
                     z=(GBPREAL *)ADaPS_fetch(plist->data,"vz_%s",species_name);
                  }
                  for(size_t i_particle=0;i_particle<n_species_local;i_particle++){
                     double r2_i;
                     r2_i =(double)(x[i_particle]-x_cen)*(double)(x[i_particle]-x_cen);
                     r2_i+=(double)(y[i_particle]-y_cen)*(double)(y[i_particle]-y_cen);
                     r2_i+=(double)(z[i_particle]-z_cen)*(double)(z[i_particle]-z_cen);
                     if(r2_i<=r2){
                        mark[i_particle]=arg->value;
                        mark_count++;
                     }
                  }
               }
               // ... set all particles in a group or subgroup to a value.
               else if(!strcmp(arg->type,"group_index") ||
                       !strcmp(arg->type,"subgroup_index")){
                  int selection_index=arg->ival[0];

                  // Set a flag to indicate whether we are processing a group or subgroup
                  int flag_process_group=TRUE;
                  if(!strcmp(arg->type,"subgroup_index"))
                     flag_process_group=FALSE;

                  // Decide which selection function to use
                  int (*select_function)(int                i_group,
                                         int                j_subgroup,
                                         int                i_subgroup,
                                         int                flag_long_ids,
                                         process_halo_info *group_i,
                                         process_halo_info *subgroup_i,
                                         void              *params);
                  int (*count_function) (int                i_group,
                                         int                j_subgroup,
                                         int                i_subgroup,
                                         int                flag_long_ids,
                                         process_halo_info *group_i,
                                         process_halo_info *subgroup_i,
                                         void              *params);
                  int (*action_function)(int                i_group,
                                         int                j_subgroup,
                                         int                i_subgroup,
                                         int                flag_long_ids,
                                         process_halo_info *group_i,
                                         process_halo_info *subgroup_i,
                                         void              *params);
                  if(flag_process_group){
                     select_function=select_group_index_local;
                     count_function =count_group_ids_local;
                     action_function=add_group_to_ids_list_local;
                  }
                  else{
                     select_function=select_subgroup_index_local;
                     count_function =count_subgroup_ids_local;
                     action_function=add_subgroup_to_ids_list_local;
                  }

                  // Create list of particles to mark
                  void   *ids_list=NULL;
                  char   *val_list=NULL;
                  size_t  n_ids_list=0;
                  make_ids_list(render,
                                i_snap,
                                select_function,
                                count_function,
                                action_function,
                                &ids_list,
                                &val_list,
                                &n_ids_list,
                                0,
                                selection_index);

                  // Mark the particles in the list
                  apply_mark_list(render,arg,i_snap,i_type,ids_list,val_list,n_ids_list,mark,&mark_count);

                  // Free the list
                  SID_free(SID_FARG ids_list);
                  SID_free(SID_FARG val_list);
               }
               else
                  SID_trap_error("Invalid selection type {%s} in perform_marking().",ERROR_LOGIC,arg->type);
               n_particles+=n_species_local;
            }
         }
      }
      SID_log("(%lld of %lld particles marked)...",SID_LOG_CONTINUE,mark_count,n_particles);
      if(render->n_interpolate>1)
         SID_log("Done.",SID_LOG_CLOSE);
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

void perform_marking(render_info *render){
   mark_arg_info *current_arg=render->mark_arg_first;
   if(current_arg!=NULL){
      SID_log("Performing particle marking...",SID_LOG_OPEN);
      while(current_arg!=NULL){
         // Perform marking
         execute_marking_argument_local(render,current_arg);
         current_arg=current_arg->next;
      }
      SID_log("Done.",SID_LOG_CLOSE);
   }
}

