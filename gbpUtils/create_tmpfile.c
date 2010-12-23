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
    int    flag;
    int    irand;
    float  frand;
    long   seed;
    FILE  *fp;

    /*********************************************/
    /* Check syntax and read-in input parameters */
    /*********************************************/
    if(argc!=1){
        printf("\n  Syntax: %s\n",
            argv[0]);
        return(SYNTAX_ERROR);
    }

    flag=TRUE;
    seed=-1*(long)times(NULL);
    while(flag){
      frand=ran1(&seed);
      irand=(int)(32768.0*frand);
      sprintf(filename,".tmpfile%d",irand);
      if((fp=fopen(filename,"r"))==NULL)
	flag=FALSE;
    }
    fp=fopen(filename,"w");close(fp);

    /******************/
    /* Report results */
    /******************/
    printf("%s\n",filename);

    return CLEAN_EXIT;
}
