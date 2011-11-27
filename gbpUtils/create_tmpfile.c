/********************************************************/
/*                                                      */
/*                   create_tmpfile.c                   */
/*                   ----------------                   */
/*                                                      */
/*  This program creates a new, random tmpfile in the   */
/*  present directory, checking to make sure that it    */
/*  does not already exist.                             */
/*                                                      */
/*  SYNTAX: find_tempfile_name                          */
/*  ------                                              */
/*                                                      */
/*  Created by:   Greg Poole                            */
/*  Last updated: Feb 10/2005                           */
/*                                                      */
/********************************************************/
#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpMath.h>
#include <gbpLib.h>

int main(int argc, char *argv[]){
    char      filename[256];
    int       flag;
    int       irand;
    int       seed;
    FILE     *fp;
    RNG_info  RNG;

    /*********************************************/
    /* Check syntax and read-in input parameters */
    /*********************************************/
    if(argc!=1){
        printf("\n  Syntax: %s\n",
            argv[0]);
        return(ERROR_SYNTAX);
    }

    // Initialize random number gerator
    init_seed_from_clock(&seed);
    init_RNG(&seed,&RNG,RNG_DEFAULT);
   
    flag=TRUE;
    while(flag){
      irand=(int)(32768.0*random_number(&RNG));
      sprintf(filename,".tmpfile%d",irand);
      if((fp=fopen(filename,"r"))==NULL)
	flag=FALSE;
    }
    fp=fopen(filename,"w");fclose(fp);
    free_RNG(&RNG);

    /******************/
    /* Report results */
    /******************/
    printf("%s\n",filename);

    return(ERROR_NONE);
}
