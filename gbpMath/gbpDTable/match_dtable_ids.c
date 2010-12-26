#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <common.h>
#include <ADaM.h>
void match_dtable_ids(dtable *table1,
		      dtable *table2,
		      char    var1[DCOLUMN_NAME_LENGTH],
		      char    var2[DCOLUMN_NAME_LENGTH],
		      char   *match_name,
		      int    *n_matches){
  int      i,j;
  int     *results;
  long    *id_1;
  long    *id_2;
  int     *index_1;
  int     *index_2;
  int      n_data_1;
  int      n_data_2;
  dcolumn *current_column;
  int      flag;

  /* Get ids for table 1*/
  if(exist_dcolumn(table1,var1)){
    n_data_1=table1->n_data;
    id_1    =(long *)fetch_dcolumn(table1,var1);
    flag    =TRUE;
  }
  else{
    fprintf(stderr,"WARNING: column {%s} does not exist in table 1!\n",var1);
    flag=FALSE;
  }

  /* Get ids for table 2*/
  if(exist_dcolumn(table2,var2)){
    n_data_2=table2->n_data;
    id_2    =(long *)fetch_dcolumn(table2,var2);
    flag   *=TRUE;
  }
  else{
    fprintf(stderr,"WARNING: column {%s} does not exist in table 2!\n",var2);
    flag=FALSE;
  }
	    
  /* Proceed if data is found ... */
  (*n_matches)=0;
  if(flag){

    /* Generate sort indices */
    sort_l(id_1,n_data_1,&index_1,FALSE);
    sort_l(id_2,n_data_2,&index_2,FALSE);
    
    /* Generate matches */
    results=(int *)malloc(sizeof(int)*n_data_1);
    for(i=0,j=0;i<n_data_1 && j<n_data_2;i++){
      while(id_2[index_2[j]]<id_1[index_1[i]] && j<n_data_2-1 ) j++;
      if(id_2[index_2[j]]==id_1[index_1[i]]){
	results[index_1[i]]=index_2[j];
	(*n_matches)++;
      }
      else
	results[index_1[i]]=0;
    }
    
    /* Store match column */
    add_dcolumn(table1,
		match_name,
		ADaM_INT,
		(void *)results,
		n_data_1);
    
    free(results);
    free(index_1);
    free(index_2);
  }
}
