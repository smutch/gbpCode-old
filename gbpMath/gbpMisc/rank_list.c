#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <common.h>
#include <ADaM.h>

void sort(double  *sval,
	  size_t   nval,
	  int      type,
	  int      flag_local,
	  size_t **index){
  int     i_rank;
  int     j_rank;
  int     rank_send_to;
  int     rank_receive_from;
  size_t  i,ir,j,l,k;
  size_t  nval_max;
  size_t  nval_tmp;
  size_t  nval_all;
  size_t *index_tmp;
  size_t *index_tmp_tmp;
  size_t *increment;
  void   *sval_tmp;
  double *sval_dtmp;

  // Rank local array in ascending order; must
  //   use merge_sort (despite extra RAM use) 
  //   because it is "stable"
  merge_sort(sval,
	     &index_tmp,
	     nval,
             type,
             SORT_COMPUTE_RANK,
             SORT_COMPUTE_INPLACE);
  heap_sort(index_tmp,
	    nval,
	    index,
	    SORT_COMPUTE_INDEX,
            SORT_COMPUTE_INPLACE);

for(i_rank=0;i_rank<SID.n_proc;i_rank++){
if(SID.My_rank==i_rank){
  for(i=0;i<nval;i++)
    fprintf(stderr,"A%d: %lf %lld\n",i_rank,sval[i],(*index)[i]);
}
fflush(stderr);MPI_Barrier(MPI_COMM_WORLD);
}

  // If this is being run in parallel, consider
  //   the arrays on the other ranks as well
  //   if flag_local is set to SORT_GLOBAL
#ifdef USE_MPI
  if(flag_local==SORT_GLOBAL && SID.n_proc>1){
    MPI_Allreduce(&nval,&nval_max,1,MPI_SIZE_T,MPI_MAX,MPI_COMM_WORLD);
    switch(type){
    case ADaM_INT:
      sval_tmp=(void *)malloc(sizeof(int)*(nval+nval_max));
      break;
    case ADaM_LONG:
      sval_tmp=(void *)malloc(sizeof(long)*(nval+nval_max));
      break;
    case ADaM_SIZE_T:
      sval_tmp=(void *)malloc(sizeof(size_t)*(nval+nval_max));
      break;
    case ADaM_DOUBLE:
      sval_tmp=(void *)malloc(sizeof(double)*(nval+nval_max));
      break;
    case ADaM_FLOAT:
      sval_tmp=(void *)malloc(sizeof(float)*(nval+nval_max));
      break;
    }
    index_tmp=(size_t *)malloc(sizeof(size_t)*(nval+nval_max));
    increment=(size_t *)malloc(sizeof(size_t)*nval);
    for(i=0;i<nval;i++) 
      increment[i]=0;
    rank_send_to     =SID.My_rank;
    rank_receive_from=SID.My_rank;
    for(i_rank=0;i_rank<SID.n_proc-1;i_rank++){
      rank_send_to++;
      rank_receive_from--;
      if(rank_receive_from<0)      rank_receive_from+=SID.n_proc;
      if(rank_send_to>=SID.n_proc) rank_send_to     -=SID.n_proc;
      MPI_Send(&nval,
	       1, 
	       MPI_SIZE_T,
	       rank_send_to,   
	       0,
	       MPI_COMM_WORLD);
      MPI_Recv(&nval_tmp,  
	       1,
	       MPI_SIZE_T,
	       rank_receive_from,
	       MPI_ANY_TAG,
	       MPI_COMM_WORLD, 
	       MPI_STATUS_IGNORE);
      MPI_Barrier(MPI_COMM_WORLD);
      nval_all=nval+nval_tmp;
      // Ranks < My_rank are placed first in the new data list
      //  so that duplicate entries are placed in a unique order
      switch(type){
      case ADaM_INT:
	MPI_Send(sval,
		 nval, 
		 MPI_INTEGER,
		 rank_send_to,   
		 0,
		 MPI_COMM_WORLD);
	if(rank_receive_from<SID.My_rank){
	  MPI_Recv(sval_tmp,  
		   nval_tmp,
		   MPI_INTEGER,
		   rank_receive_from,
		   MPI_ANY_TAG,
		   MPI_COMM_WORLD, 
		   MPI_STATUS_IGNORE);
	  for(i=nval_tmp;i<nval_all;i++)
	    ((int *)sval_tmp)[i]=((int *)sval)[i-nval_tmp];
	}
	else{
	  MPI_Recv(&(((int *)sval_tmp)[nval]),  
		   nval_tmp,
		   MPI_INTEGER,
		   rank_receive_from,
		   MPI_ANY_TAG,
		   MPI_COMM_WORLD, 
		   MPI_STATUS_IGNORE);
	  for(i=0;i<nval;i++)
	    ((int *)sval_tmp)[i]=((int *)sval)[i];
	}
	break;
      case ADaM_LONG:
	MPI_Send(sval,
		 nval, 
		 MPI_LONG,
		 rank_send_to,   
		 0,
		 MPI_COMM_WORLD);
	if(rank_receive_from<SID.My_rank){
	  MPI_Recv(sval_tmp,  
		   nval_tmp,
		   MPI_LONG,
		   rank_receive_from,
		   MPI_ANY_TAG,
		   MPI_COMM_WORLD, 
		   MPI_STATUS_IGNORE);
	  for(i=nval_tmp;i<nval_all;i++)
	    ((long *)sval_tmp)[i]=((long *)sval)[i-nval_tmp];
	}
	else{
	  MPI_Recv(&(((long *)sval_tmp)[nval]),  
		   nval_tmp,
		   MPI_LONG,
		   rank_receive_from,
		   MPI_ANY_TAG,
		   MPI_COMM_WORLD, 
		   MPI_STATUS_IGNORE);
	  for(i=0;i<nval;i++)
	    ((long *)sval_tmp)[i]=((long *)sval)[i];
	}
	break;
      case ADaM_SIZE_T:
	MPI_Send(sval,
		 nval, 
		 MPI_SIZE_T,
		 rank_send_to,   
		 0,
		 MPI_COMM_WORLD);
	if(rank_receive_from<SID.My_rank){
	  MPI_Recv(sval_tmp,  
		   nval_tmp,
		   MPI_SIZE_T,
		   rank_receive_from,
		   MPI_ANY_TAG,
		   MPI_COMM_WORLD, 
		   MPI_STATUS_IGNORE);
	  for(i=nval_tmp;i<nval_all;i++)
	    ((size_t *)sval_tmp)[i]=((size_t *)sval)[i-nval_tmp];
	}
	else{
	  MPI_Recv(&(((size_t *)sval_tmp)[nval]),  
		   nval_tmp,
		   MPI_SIZE_T,
		   rank_receive_from,
		   MPI_ANY_TAG,
		   MPI_COMM_WORLD, 
		   MPI_STATUS_IGNORE);
	  for(i=0;i<nval;i++)
	    ((size_t *)sval_tmp)[i]=((size_t *)sval)[i];
	}
	break;
      case ADaM_DOUBLE:
	MPI_Send(sval,
		 nval, 
		 MPI_DOUBLE,
		 rank_send_to,   
		 0,
		 MPI_COMM_WORLD);
	if(rank_receive_from<SID.My_rank){
	  MPI_Recv(sval_tmp,  
		   nval_tmp,
		   MPI_DOUBLE,
		   rank_receive_from,
		   MPI_ANY_TAG,
		   MPI_COMM_WORLD, 
		   MPI_STATUS_IGNORE);
	  for(i=nval_tmp;i<nval_all;i++)
	    ((double *)sval_tmp)[i]=((double *)sval)[i-nval_tmp];
	}
	else{
	  MPI_Recv(&(((double *)sval_tmp)[nval]),  
		   nval_tmp,
		   MPI_DOUBLE,
		   rank_receive_from,
		   MPI_ANY_TAG,
		   MPI_COMM_WORLD, 
		   MPI_STATUS_IGNORE);
	  for(i=0;i<nval;i++)
	    ((double *)sval_tmp)[i]=((double *)sval)[i];
	}
	break;
      case ADaM_FLOAT:
	MPI_Send(sval,
		 nval, 
		 MPI_FLOAT,
		 rank_send_to,   
		 0,
		 MPI_COMM_WORLD);
	if(rank_receive_from<SID.My_rank){
	  MPI_Recv(sval_tmp,  
		   nval_tmp,
		   MPI_FLOAT,
		   rank_receive_from,
		   MPI_ANY_TAG,
		   MPI_COMM_WORLD, 
		   MPI_STATUS_IGNORE);
	  for(i=nval_tmp;i<nval_all;i++)
	    ((float *)sval_tmp)[i]=((float *)sval)[i-nval_tmp];
	}
	else{
	  MPI_Recv(&(((float *)sval_tmp)[nval]),  
		   nval_tmp,
		   MPI_FLOAT,
		   rank_receive_from,
		   MPI_ANY_TAG,
		   MPI_COMM_WORLD, 
		   MPI_STATUS_IGNORE);
	  for(i=0;i<nval;i++)
	    ((float *)sval_tmp)[i]=((float *)sval)[i];
	}
	break;
      }
      MPI_Barrier(MPI_COMM_WORLD);

      /*
      sort((double *)sval_tmp,
	     nval_all,
	     &index_tmp,
	     FALSE);
      */
  sval_dtmp=(double *)malloc(sizeof(double)*nval_all);
  for(i=0;i<nval_all;i++) sval_dtmp[i]=((double *)sval_tmp)[i];
  merge_sort((double *)sval_dtmp,
	     &index_tmp,
	     nval_all);
  free(sval_dtmp);

      sort_s(index_tmp,
	     nval_all,
	     &index_tmp_tmp,
	     TRUE);
      free(index_tmp);

      // Figure-out where the local values lie in the new sort
      //   and increment rank for each if any values from
      //   receive_from_rank's array are less than it
      if(rank_receive_from<SID.My_rank){
	for(i=0,j=0;i<nval;i++)
	  increment[i]+=index_tmp_tmp[i+nval_tmp]-(*index)[i];
for(j_rank=0;j_rank<SID.n_proc;j_rank++){
if(SID.My_rank==j_rank){
  for(i=0;i<nval_all;i++){
    if(i>nval_tmp)
      fprintf(stderr,"B%2d->%2d: %lf %lld %lld\n",rank_receive_from,j_rank,((double *)sval_tmp)[i],index_tmp_tmp[i],increment[i-nval_tmp]);
    else
      fprintf(stderr,"B%2d->%2d: %lf %lld\n",rank_receive_from,j_rank,((double *)sval_tmp)[i],index_tmp_tmp[i]);
  }
}
fflush(stderr);MPI_Barrier(MPI_COMM_WORLD);
}
      }
      else{
	for(i=0,j=0;i<nval;i++)
	  increment[i]+=index_tmp_tmp[i]-(*index)[i];
for(j_rank=0;j_rank<SID.n_proc;j_rank++){
if(SID.My_rank==j_rank){
  for(i=0;i<nval_all;i++){
    if(i<nval)
      fprintf(stderr,"B%2d->%2d: %lf %lld %lld\n",rank_receive_from,j_rank,((double *)sval_tmp)[i],index_tmp_tmp[i],increment[i]);
    else
      fprintf(stderr,"B%2d->%2d: %lf %lld\n",rank_receive_from,j_rank,((double *)sval_tmp)[i],index_tmp_tmp[i]);
  }
}
fflush(stderr);MPI_Barrier(MPI_COMM_WORLD);
}
      }
      free(index_tmp_tmp);
    }
    for(i=0;i<nval;i++)
      (*index)[i]+=increment[i];

    free(increment);
    free(sval_tmp);
  }
#endif
for(i_rank=0;i_rank<SID.n_proc;i_rank++){
  if(i_rank==SID.My_rank){
    for(i=0;i<nval;i++)
      fprintf(stderr,"C%d: %lf %lld\n",SID.My_rank,sval[i],(*index)[i],increment[i]);
  }
  fflush(stderr);MPI_Barrier(MPI_COMM_WORLD);
}

}
