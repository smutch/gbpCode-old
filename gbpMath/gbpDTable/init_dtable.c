#include <stdio.h>
#include <common.h>
void init_dtable(dtable *table){
  table->n_data   =0;
  table->n_columns=0;
  table->remarks  =NULL;
  table->columns  =NULL;
}
