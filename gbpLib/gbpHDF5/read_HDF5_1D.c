#include <gbpLib.h>
#include <gbpHDF5.h>

//######################################
// For reading 1D arrays from hdf5 files
//######################################
// arrayName - The name of the dataset to be read
// filename  - The location of the hdf5 file
// data_out  - A 1D float array (not initialised)
// themin	 - Function sets this value to the minimum value of the dataset
// themax	 - Function sets this value to the maximum value of the dataset
// count	 - How many values are present in the dataset
//
int read_HDF5_1D(char *arrayName, char *filename, float **data_out, float *themin, float *themax, int *count)
{
 hid_t hdf5_file, hdf5_headergrp, hdf5_gasgrp;
 hid_t hdf5_dataspace_in_file, hdf5_dataset;
 hid_t hdf5_datatype = 0, hdf5_dataspace_in_memory = 0;
 
 char buf[256];
 int rank, i, j;
 hsize_t dims[2];
 size_t n_values;
 size_t nrow;
 herr_t status;

 fprintf(stderr, "Finished declaring variables\n");

 /* open a random HDF5 file and choose the gas group */
 hdf5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

 //-----------------------------------
 fprintf(stderr, "Opened hdf5 file\n");

 hdf5_headergrp = H5Gopen(hdf5_file, "Structure", H5P_DEFAULT);

 //------------------------------------
 fprintf(stderr, "Opened Particles\n");

 hdf5_gasgrp = H5Gopen(hdf5_headergrp, "Subhalo", H5P_DEFAULT);

 //-----------------------------------------
 fprintf(stderr, "Opened /Structure/Subhalo\n");
 
  /* All datasets are floats */
 hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);

 /* Define the dataset of interest */
// buf = (*arrayName);
 sprintf(buf, "%s", arrayName);
//
//
//
//
//
 

  /* Open the dataset of interest */
 hdf5_dataset = H5Dopen(hdf5_gasgrp, buf, H5P_DEFAULT);

 fprintf(stderr, "%s data set opened\n",buf);

  /* Now we can either define the file size by hand (shown later) or as in this section we will query HDF5 */
 hdf5_dataspace_in_file = H5Dget_space(hdf5_dataset);

 
 // Only position is 3D, i.e. rank = 2. Rest is rank = 1
 rank = H5Sget_simple_extent_ndims(hdf5_dataspace_in_file); 
// Only position has any extra arrays in the 2 dimension, i.e. dims[1] = 3 for x,y,z while rest is dims[1] = 1
// All arrays have dims[0] = nparticles
 status = H5Sget_simple_extent_dims(hdf5_dataspace_in_file, dims, NULL);  
 
/* Malloc data_out */
	if(( *data_out =(float *)malloc(dims[0]*sizeof(float))) == NULL) {
			fprintf(stderr, "Malloc of %s array failed\n", buf);
			exit(-1);
	}

  /* Position has x,y,z while other arrays are just vectors */
 dims[1] = 1; 

/* Data space in file */
 hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
   
 /* Space in memory */ 
 hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);

 status = H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, *data_out);

 
 // Here I'm just going to print out the array. Kinda pointless.
  /* print it by rows */
 n_values = (size_t)(dims[0] * dims[1]);

 fprintf(stderr,"dims[0]=%ld\n", (long)dims[0]);
 fprintf(stderr,"dims[1]=%ld\n", (long)dims[1]);
 fprintf(stderr,"%ld nrow\n",    (long)dims[1]);

 nrow = (size_t)dims[1];

 fprintf(stderr, "n_values=%ld\n", (long)n_values); 
 fprintf(stderr, "data_out[10]=%f\n", (*data_out)[10]); 

 for (i=0; i<n_values/nrow; i++ )
 {
  	for (j=0; j<nrow; j++){
                if(themin!=NULL)
		  *themin = MIN(*themin, (*data_out)[i*nrow+j]);
                if(themax!=NULL)
		  *themax = MAX(*themax, (*data_out)[i*nrow+j]);
		(*count)++;
	}
 }

 /* Close the dataset and dataspace links in reverse order */ 
 H5Sclose(hdf5_dataspace_in_memory);
 H5Sclose(hdf5_dataspace_in_file);
 H5Dclose(hdf5_dataset);


/* Close the groups and file in the reverse order */ 
 H5Tclose(hdf5_datatype);
 H5Gclose(hdf5_gasgrp);
 H5Gclose(hdf5_headergrp);
 H5Fclose(hdf5_file);

 return 1;
}
			  
