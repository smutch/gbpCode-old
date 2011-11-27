/********************************************************/
/*                                                      */
/*                    random_number.c                   */
/*                    ---------------                   */
/*                                                      */
/********************************************************/
#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
int main(int argc, char *argv[]){
    char      filename[256];
    int       i,j,n_l,n_c;
    int       seed;
    FILE     *fp;
    RNG_info  RNG;

    // Check syntax and read-in input parameters 
    if(argc!=3){
        printf("\n  Syntax: %s number_of_lines number_of_columns\n",
            argv[0]);
        return(ERROR_SYNTAX);
    }
    n_l=(int)atoi(argv[1]);
    n_c=(int)atoi(argv[2]);

    // Generate and print results 
    init_RNG(&seed,&RNG,RNG_DEFAULT);
    for(i=0;i<n_l;i++){
       for(j=0;j<n_c;j++)
         printf("%f ",random_number(&RNG));
       printf("\n");
    }
    free_RNG(&RNG);
    return(ERROR_NONE);
}
