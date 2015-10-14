#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <gbpLib.h>
#include <gbpSPH.h>

// In the future ... receive a list of variables to make this adaptive

void write_gadget_ascii(char       *filename_out,
	                    plist_info *plist,
                        size_t     *write_order_indices){
  FILE   *fp;
  double  h_Hubble;
  GBPREAL   *x;
  GBPREAL   *y;
  GBPREAL   *z;
  GBPREAL   *vx;
  GBPREAL   *vy;
  GBPREAL   *vz;
  size_t    *id;
  size_t  n_p;
  int     i_rank;
  int     i_species;

  SID_log("Writing ascii file {%s}...",SID_LOG_OPEN,filename_out);

  // Process each rank in turn
  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
    if(i_rank==SID.My_rank){

      h_Hubble=((double *)ADaPS_fetch(plist->data,"h_Hubble"))[0];

      // Open file
      if(i_rank==MASTER_RANK){
        fp=fopen(filename_out,"w");
        fprintf(fp,"#Columns:\n");
        fprintf(fp,"#  1) Gadget particle type\n");
        fprintf(fp,"#  2) x   [Mpc/h]\n");
        fprintf(fp,"#  3) y   [Mpc/h]\n");
        fprintf(fp,"#  4) z   [Mpc/h]\n");
        fprintf(fp,"#  5) v_x [km/s]\n");
        fprintf(fp,"#  6) v_y [km/s]\n");
        fprintf(fp,"#  7) v_z [km/s]\n");
        fprintf(fp,"#  8) id\n");
      }
      else
        fp=fopen(filename_out,"a");

      // Write this rank's particles
      int i_p=0;
      for(i_species=0;i_species<N_GADGET_TYPE;i_species++){
        if(ADaPS_exist(plist->data,"n_%s",plist->species[i_species])){
          n_p=((size_t *)ADaPS_fetch(plist->data,"n_%s",plist->species[i_species]))[0];
          if(n_p>0){
             x =(GBPREAL *)ADaPS_fetch(plist->data,"x_%s", plist->species[i_species]);
             y =(GBPREAL *)ADaPS_fetch(plist->data,"y_%s", plist->species[i_species]);
             z =(GBPREAL *)ADaPS_fetch(plist->data,"z_%s", plist->species[i_species]);
             vx=(GBPREAL *)ADaPS_fetch(plist->data,"vx_%s",plist->species[i_species]);
             vy=(GBPREAL *)ADaPS_fetch(plist->data,"vy_%s",plist->species[i_species]);
             vz=(GBPREAL *)ADaPS_fetch(plist->data,"vz_%s",plist->species[i_species]);
             id=(size_t  *)ADaPS_fetch(plist->data,"id_%s",plist->species[i_species]);
             for(int j_p=0;j_p<n_p;i_p++,j_p++){
               size_t k_p=j_p;
               if(write_order_indices!=NULL)
                  k_p=write_order_indices[i_p];
               fprintf(fp,"%1d %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %7zd\n",
                       i_species,
                       (double)(x[k_p]/(M_PER_MPC/h_Hubble)),
                       (double)(y[k_p]/(M_PER_MPC/h_Hubble)),
                       (double)(z[k_p]/(M_PER_MPC/h_Hubble)),
                       (double)(vx[k_p]*1e-3),
                       (double)(vy[k_p]*1e-3),
                       (double)(vz[k_p]*1e-3),
                       id[k_p]);
             }
          }
        }
      }
      fclose(fp);
    }
    SID_Barrier(SID.COMM_WORLD);
  }
  SID_log("Done.",SID_LOG_CLOSE);
}

