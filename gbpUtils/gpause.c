#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[]){
  long dt_wait;
  long t_start;
  long t,dt;

  if(argc!=2){
    fprintf(stderr,"\nSyntax: %s t_to_wait_in_seconds\n\n",argv[0]);
    return(100);
  }
  dt_wait=(long)atoi(argv[1]);

  t_start=(long)time(NULL);
  t=t_start;
  while(dt_wait>(t-t_start)) t=(long)time(NULL);

  return(0);
}

