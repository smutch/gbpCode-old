#include <stdio.h>
#include <limits.h>
#include <common.h>

size_t read_binary_by_rank(void    *array,
			   int      type,
			   size_t   n_read,
			   int      My_rank,
			   FILE    *fp){
  int     i;
  int     i_rank;
  int     j_rank;
  size_t  n_remaining_local;
  size_t  n_remaining_read;
  size_t  n_remaining;
  size_t  n_buffer;
  size_t  n_buffer_max;
  size_t  n_actual;
  size_t  index=0;
  size_t  bytes_per_read;
  void   *buffer=NULL;
  int     flag_continue;
  int     master_rank;
#ifdef USE_MPI
  MPI_Datatype mpi_type;
  MPI_Status   status;
#endif
  int flag_error=FALSE;


  switch(type){
  case ADaM_INT:
    bytes_per_read=sizeof(int);
#ifdef USE_MPI
    mpi_type      =MPI_INT;
#endif
    break;
  case ADaM_LONG:
    bytes_per_read=sizeof(long int);
#ifdef USE_MPI
    mpi_type      =MPI_LONG;
#endif
    break;
  case ADaM_FLOAT:
    bytes_per_read=sizeof(float);
#ifdef USE_MPI
    mpi_type      =MPI_FLOAT;
#endif
    break;
  case ADaM_DOUBLE:
    bytes_per_read=sizeof(double);
#ifdef USE_MPI
    mpi_type      =MPI_DOUBLE;
#endif
    break;
  default:
    // Unknown data type
    flag_error=TRUE;
    break;
  }

#ifdef USE_MPI
  MPI_Allreduce(&My_rank,&master_rank,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD);
#else
  master_rank=MASTER_RANK;
#endif

  printf("SID.My_rank=%d My_rank=%d master_rank=%d\n",SID.My_rank,My_rank,master_rank);

  // BUFSIZ is defined as "the optimal read size for this platform"
  n_buffer_max     =BUFSIZ/bytes_per_read;
  n_remaining_local=n_read;
  n_remaining      =LONG_MAX;
  buffer           =malloc(n_buffer_max*bytes_per_read);

  // Read for each rank in turn
  for(i_rank=0;i_rank<SID.n_proc && !flag_error;i_rank++){
#ifdef USE_MPI
    if(i_rank!=master_rank){
      if(My_rank==master_rank)
	MPI_Send(&n_remaining_local,
		 1, 
		 MPI_SIZE_T,
		 i_rank,   
		 0,
		 MPI_COMM_WORLD);
      else if(My_rank==i_rank)
	MPI_Recv(&n_remaining_read,
		 1,
		 MPI_SIZE_T,
		 master_rank,
		 MPI_ANY_TAG,
		 MPI_COMM_WORLD, 
		 &status);
    }
#endif
    if(i_rank==master_rank)
      n_remaining_read=n_remaining_local;

    // Loop until n_read for this rank has been read
    flag_continue=TRUE;
    while(flag_continue && !flag_error){

      // Load buffer
      if(i_rank==master_rank){
	n_buffer=MIN(n_buffer_max,n_remaining_read);
	if(n_buffer<=0)
	  flag_continue=FALSE;
	else{
	  n_actual=fread(buffer,n_buffer,bytes_per_read,fp);
	  if(n_actual!=n_buffer)
	    flag_error=TRUE;
	  n_remaining_read-=n_actual;
	}
      }

      // Make sure something was read ...
#ifdef USE_MPI
      if(i_rank!=master_rank){
	if(My_rank==master_rank)
	  MPI_Send(&n_actual,
		   1, 
		   MPI_SIZE_T,
		   i_rank,   
		   0,
		   MPI_COMM_WORLD);
	else if(My_rank==i_rank)
	  MPI_Recv(&n_actual,  
		   1,
		   MPI_SIZE_T,
		   master_rank,
		   MPI_ANY_TAG,
		   MPI_COMM_WORLD, 
		   &status);
      }
#endif

      // ... if it was then send the buffer to i_rank ...
      if(n_actual>0){
	if(i_rank!=master_rank){
	  if(My_rank==master_rank)
	    MPI_Send(&buffer,
		     n_actual, 
		     mpi_type,
		     i_rank,   
		     0,
		     MPI_COMM_WORLD);
	  else if(My_rank==i_rank)
	    MPI_Recv(&buffer,  
		     n_actual,
		     mpi_type,
		     master_rank,
		     MPI_ANY_TAG,
		     MPI_COMM_WORLD, 
		     &status);
	}
	if(My_rank==master_rank){
	  // ... and then dump it to the array
	  switch(type){
	  case ADaM_INT:
	    for(i=0;i<n_actual;i++)
	      ((int *)array)[index++]=((int *)buffer)[i];
	    break;
	  case ADaM_LONG:
	    for(i=0;i<n_actual;i++)
	      ((long *)array)[index++]=((long *)buffer)[i];
	    break;
	  case ADaM_FLOAT:
	    for(i=0;i<n_actual;i++)
	      ((float *)array)[index++]=((float *)buffer)[i];
	    break;
	  case ADaM_DOUBLE:
	    for(i=0;i<n_actual;i++)
	      ((double *)array)[index++]=((double *)buffer)[i];
	    break;
	  }
	}
      }
    } // while
    if(i_rank==master_rank)
      n_remaining_local=n_remaining_read;

    // Update this rank's n_remaining
#ifdef USE_MPI
    if(i_rank!=master_rank){
      if(My_rank==master_rank)
	MPI_Send(&n_remaining_read,
		 1, 
		 MPI_SIZE_T,
		 i_rank,   
		 0,
		 MPI_COMM_WORLD);
      else if(My_rank==i_rank)
	MPI_Recv(&n_remaining_local,  
		 1,
		 MPI_SIZE_T,
		 master_rank,
		 MPI_ANY_TAG,
		 MPI_COMM_WORLD, 
		 &status);
    }
    MPI_Bcast(&flag_error,1,MPI_INT,master_rank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

#ifdef USE_MPI
  MPI_Allreduce(&n_remaining_local,&n_remaining,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
#else
  n_remaining=n_remaining_local;
#endif

  // Clean-up
  if(buffer!=NULL)
    free(buffer);

  return(n_remaining);
}


