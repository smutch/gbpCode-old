#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void read_AHF_for_trees(char       *filename_root,
												plist_info *plist,
												char       *catalog_name,
												int         flag_binary){
	int     n_halos;
	size_t  n_halos_dim;
	float  *x_cen;
	float  *y_cen;
	float  *z_cen;
	float  *vx_cen;
	float  *vy_cen;
	float  *vz_cen;
	float  *M_vir;
	float  *V_max;
	float  *sigma_V;
	float  *spin_x;
	float  *spin_y;
	float  *spin_z;
	float   L;
	float   lambda;
	int     i_halo;
	FILE   *fp;
	char   *line;
	int     line_length=0;
	char    parm_txt[256];
	char    filename_cat[256];
	
	read_AHF_groups(filename_root,
		plist,
		catalog_name,
		flag_binary);
	/*
	sprintf(filename_cat,"%s.AHF_halos",filename_root);
	fprintf(stderr,"Reading halo statistics for {%s}...",filename_root);
	if((fp=fopen(filename_cat,"r"))!=NULL){
			n_halos=(size_t)count_lines_data(fp);
			x_cen  =(float *)malloc(sizeof(float)*n_halos);
			y_cen  =(float *)malloc(sizeof(float)*n_halos);
			z_cen  =(float *)malloc(sizeof(float)*n_halos);
			vx_cen =(float *)malloc(sizeof(float)*n_halos);
			vy_cen =(float *)malloc(sizeof(float)*n_halos);
			vz_cen =(float *)malloc(sizeof(float)*n_halos);
			M_vir  =(float *)malloc(sizeof(float)*n_halos);
			V_max  =(float *)malloc(sizeof(float)*n_halos);
			sigma_V=(float *)malloc(sizeof(float)*n_halos);
			spin_x =(float *)malloc(sizeof(float)*n_halos);
			spin_y =(float *)malloc(sizeof(float)*n_halos);
			spin_z =(float *)malloc(sizeof(float)*n_halos);
			for(i_halo=0;i_halo<n_halos;i_halo++){
				grab_next_line(fp,&line,&line_length);
				grab_float(line,3, &(x_cen[i_halo]));			
				grab_float(line,4, &(y_cen[i_halo]));			
				grab_float(line,5, &(z_cen[i_halo]));			
				grab_float(line,6, &(vx_cen[i_halo]));			
				grab_float(line,7, &(vy_cen[i_halo]));			
				grab_float(line,8, &(vz_cen[i_halo]));			
				grab_float(line,9, &(M_vir[i_halo]));
				grab_float(line,11,&(V_max[i_halo]));			
				grab_float(line,13,&(sigma_V[i_halo]));			
				grab_float(line,14,&lambda);			
				grab_float(line,15,&(spin_x[i_halo]));			
				grab_float(line,16,&(spin_y[i_halo]));			
				grab_float(line,17,&(spin_z[i_halo]));			
				L=sqrt(spin_x[i_halo]*spin_x[i_halo]+
							 spin_y[i_halo]*spin_y[i_halo]+
							 spin_z[i_halo]*spin_z[i_halo]);
				spin_x[i_halo]*=lambda/L;
				spin_y[i_halo]*=lambda/L;
				spin_z[i_halo]*=lambda/L;				
			}
	}
	fprintf(stderr,"Done.\n");
	
	// Store results
	n_halos_dim=(size_t)n_halos;
	fprintf(stderr,"Storing halo information for {%s}...",filename_root);
	sprintf(parm_txt,"n_halos_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(&n_halos),
		SID_INT,
		0,
		NULL);
	sprintf(parm_txt,"x_cen_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(x_cen),
		SID_FLOAT,
		1,
		&n_halos_dim);
	free(x_cen);
	sprintf(parm_txt,"y_cen_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(y_cen),
		SID_FLOAT,
		1,
		&n_halos_dim);
	free(y_cen);
	sprintf(parm_txt,"z_cen_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(z_cen),
		SID_FLOAT,
		1,
		&n_halos_dim);
	free(z_cen);
	sprintf(parm_txt,"vx_cen_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(vx_cen),
		SID_FLOAT,
		1,
		&n_halos_dim);
	free(vx_cen);
	sprintf(parm_txt,"vy_cen_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(vy_cen),
		SID_FLOAT,
		1,
		&n_halos_dim);
	free(vy_cen);
	sprintf(parm_txt,"vz_cen_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(vz_cen),
		SID_FLOAT,
		1,
		&n_halos_dim);
	free(vz_cen);
	sprintf(parm_txt,"M_vir_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(M_vir),
		SID_FLOAT,
		1,
		&n_halos_dim);
	free(M_vir);
	sprintf(parm_txt,"V_max_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(V_max),
		SID_FLOAT,
		1,
		&n_halos_dim);
	free(V_max);
	sprintf(parm_txt,"sigma_V_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(sigma_V),
		SID_FLOAT,
		1,
		&n_halos_dim);
	free(sigma_V);
	sprintf(parm_txt,"spin_x_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(spin_x),
		SID_FLOAT,
		1,
		&n_halos_dim);
	free(spin_x);
	sprintf(parm_txt,"spin_y_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(spin_y),
		SID_FLOAT,
		1,
		&n_halos_dim);
	free(spin_y);
	sprintf(parm_txt,"spin_z_%s",catalog_name);
	SID_store(&(plist->data),
		parm_txt,
		(void *)(spin_z),
		SID_FLOAT,
		1,
		&n_halos_dim);
	free(spin_z);
	fprintf(stderr,"Done.\n");
	*/
}
