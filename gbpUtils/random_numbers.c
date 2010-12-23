/********************************************************/
/*                                                      */
/*                    random_number.c                   */
/*                    ---------------                   */
/*                                                      */
/********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float ran1(long *seed);
#define  CLEAN_EXIT          0
#define  SYNTAX_ERROR       -1
#ifndef  TRUE
#define  TRUE 1
#endif
#ifndef  FALSE
#define  FALSE 0
#endif
int main(int argc, char *argv[]){
    char   filename[256];
    int    i,j,n_l,n_c;
    long   seed;
    FILE  *fp;

    // Check syntax and read-in input parameters 
    if(argc!=3){
        printf("\n  Syntax: %s number_of_lines number_of_columns\n",
            argv[0]);
        return(SYNTAX_ERROR);
    }
    n_l=(int)atoi(argv[1]);
    n_c=(int)atoi(argv[2]);

    // Generate and print results 
    seed=-1*(long)times(NULL);
    for(i=0;i<n_l;i++){
       for(j=0;j<n_c;j++)
         printf("%f ",ran1(&seed));
       printf("\n");
    }
    return(CLEAN_EXIT);
}
