#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>

void read_gadget_binary_header(char        *filename_root_in,
                               int          snapshot_number,
                               plist_info  *plist){
  FILE   *fp;
  char  **pname;
  size_t  i,j,k; 
  size_t  n_of_type[N_GADGET_TYPE];
  unsigned int n_of_type_tmp[N_GADGET_TYPE];
  size_t  n_particles;
  double  mass_array[N_GADGET_TYPE];
  int     unused[256];
  int     n_all_tmp[N_GADGET_TYPE];
  size_t  n_all[N_GADGET_TYPE];
  int     n_files;
  int     flag_metals;
  int     flag_ages;
  int     flag_entropyICs;
  double  h_Hubble;
  double  box_size;
  double  d_value;
  size_t  s_value;
  double *d_array;
  double *d1_array;
  double *d2_array;
  double *d3_array;
  float   f_temp;
  float   f1_temp;
  float   f2_temp;
  float   f3_temp;
  int     i_value;
  int     i_temp;
  int     i1_temp;
  int     i2_temp;
  char    filename[256];
  int     n_type_used;
  int     record_length;
  int     n_return;
  int     s_load;
  int     flag_alloc_d1_array;
  int     i_file;
  int     flag_filefound=FALSE;
  int     flag_multifile=FALSE;
  int     flag_file_type=0;

  gadget_read_info fp_gadget;
  flag_filefound=init_gadget_read(filename_root_in,snapshot_number,&fp_gadget);

  // A file was found ... read the header
  if(flag_filefound){
    flag_multifile=fp_gadget.flag_multifile;
    flag_file_type=fp_gadget.flag_file_type;
    set_gadget_filename(&fp_gadget,0,filename);
    fp=fopen(filename,"r");

    // Grab species names from the plist header
    pname=plist->species;

    // Read header stuff...

    //  ...header record length...
    fread_verify(&record_length,4,1,fp);
    s_load=0;

    //  ...particle numbers for each species...
    n_particles =0;
    n_return    =fread_verify(n_of_type_tmp,sizeof(unsigned int),N_GADGET_TYPE,fp); 
    s_load     +=n_return*sizeof(unsigned int);
    for(i=0,n_type_used=0;i<N_GADGET_TYPE;i++) {
      n_of_type[i]=(size_t)n_of_type_tmp[i];
      if(n_of_type[i]>0){
        ADaPS_store(&(plist->data),(void *)(&(n_of_type[i])),"n_%s",ADaPS_SCALAR_SIZE_T,pname[i]);
        n_particles+=n_of_type[i];
        n_type_used++;
      }
      else
        n_of_type[i]=0;
    }
    ADaPS_store(&(plist->data),(void *)(&n_particles),"n_particles",ADaPS_SCALAR_SIZE_T);

    //  ...particle mass array...
    n_return =fread_verify(mass_array,sizeof(double),N_GADGET_TYPE,fp); 
    s_load  +=n_return*sizeof(double);

    //  ...expansion factor (or time)...
    n_return=fread_verify(&d_value,sizeof(double),1,fp);
    s_load+=n_return*sizeof(double);
    ADaPS_store(&(plist->data),(void *)(&d_value),"expansion_factor",ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&(plist->data),(void *)(&d_value),"time",            ADaPS_SCALAR_DOUBLE);

    //  ...redshift...
    n_return=fread_verify(&d_value,sizeof(double),1,fp);
    s_load+=n_return*sizeof(double);
    ADaPS_store(&(plist->data),(void *)(&d_value),"redshift",ADaPS_SCALAR_DOUBLE);

    //  ...some flags...
    n_return=fread_verify(&i_value,sizeof(int),1,fp);
    s_load+=n_return*sizeof(int);
    ADaPS_store(&(plist->data),(void *)(&i_value),"flag_Sfr",ADaPS_SCALAR_INT);
    n_return=fread_verify(&i_value,sizeof(int),1,fp);
    s_load+=n_return*sizeof(int);
    ADaPS_store(&(plist->data),(void *)(&i_value),"flag_feedback",ADaPS_SCALAR_INT);

    //  ...number of particles in all files...
    n_return=fread_verify(n_all_tmp,sizeof(int),N_GADGET_TYPE,fp);
    for(i=0;i<N_GADGET_TYPE;i++)
      n_all[i]=(size_t)n_all_tmp[i];
    s_load+=n_return*sizeof(int);

    //  ...another flag...
    n_return=fread_verify(&i_value,sizeof(int),1,fp);
    s_load+=n_return*sizeof(int);
    ADaPS_store(&(plist->data),(void *)(&i_value),"flag_cooling",ADaPS_SCALAR_INT);

    //  ...number of files per snapshot...
    n_return=fread_verify(&i_value,sizeof(int),1,fp);
    s_load+=n_return*sizeof(int);
    ADaPS_store(&(plist->data),(void *)(&i_value),"n_files",ADaPS_SCALAR_INT);

    //  ...box size (n.b.: we need h_Hubble before we can store box_size...do so later)...
    n_return=fread_verify(&box_size,sizeof(double),1,fp);
    s_load+=n_return*sizeof(double);

    // Cosmology

    //  ...Omega_o...
    n_return=fread_verify(&d_value,sizeof(double),1,fp);
    s_load+=n_return*sizeof(double);
    ADaPS_store(&(plist->data),(void *)(&d_value),"Omega_M",ADaPS_SCALAR_DOUBLE);

    //  ...Omega_Lambda...
    n_return=fread_verify(&d_value, sizeof(double),1,fp);
    s_load+=n_return*sizeof(double);
    ADaPS_store(&(plist->data),(void *)(&d_value),"Omega_Lambda",ADaPS_SCALAR_DOUBLE);

    //  ...Hubble parameter...
    n_return=fread_verify(&h_Hubble,sizeof(double),1,fp);
    s_load+=n_return*sizeof(double);
    box_size*=GADGET_LENGTH/h_Hubble;
    ADaPS_store(&(plist->data),(void *)(&box_size),"box_size",ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&(plist->data),(void *)(&h_Hubble),"h_Hubble",ADaPS_SCALAR_DOUBLE);

    // Last few things in the header
    n_return=fread_verify(&flag_metals,sizeof(int),1,fp);
    s_load+=n_return*sizeof(int);
    n_return=fread_verify(&flag_ages,sizeof(int),1,fp);
    s_load+=n_return*sizeof(int);
    n_return=fread_verify(&n_all_tmp,sizeof(int),N_GADGET_TYPE,fp);
    s_load+=n_return*sizeof(int);
    n_return=fread_verify(&flag_entropyICs,sizeof(int),1,fp);
    s_load+=n_return*sizeof(int);
    for(i=0;i<N_GADGET_TYPE;i++)
      n_all[i]+=(((size_t)n_all_tmp[i]) << 32);

    // Store n_all array
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_all[i]>0){
        s_value=n_all[i];
        ADaPS_store(&(plist->data),(void *)(&s_value),"n_all_%s",ADaPS_SCALAR_SIZE_T,pname[i]);
      }
    }

    // Store mass_array
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_all[i]>0){
        mass_array[i]*=GADGET_MASS/h_Hubble;
        d_value=mass_array[i];
        ADaPS_store(&(plist->data),(void *)(&d_value),"mass_array_%s",ADaPS_SCALAR_DOUBLE,pname[i]);
      }
    }

    fclose(fp);
  }
  else
    fprintf(stderr,"Error opening file {%s}.\n",filename_root_in);

};
