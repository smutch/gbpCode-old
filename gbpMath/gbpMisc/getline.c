#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <common.h>
int getline(char **string,int *n,FILE *fp){
	size_t  size;
	size_t  size_old;
	size_t  len  = 0;
	size_t  last = 0;
        int     i;
	int     start_flag=TRUE;
	size=(size_t)(*n);
	if(size<=0){
		size   =BUFSIZ; /* BUFSIZ is defined as "the optimal read size for this platform" */
		*string=NULL;
	}
	do {
		size_old=size;
		if(!start_flag)
			size+=BUFSIZ;
		*string=(char *)realloc(*string,size); // realloc(NULL,n) is the same as malloc(n)
                if(start_flag)
                   strcpy(*string,"\0");
		start_flag=FALSE;
		fgets(*string+len,BUFSIZ,fp);
		len =strlen(*string)-1; // Ignore trailing '\0's
		last=len-1; /* Remove trailing '\0's */
	} while (!feof(fp) && (*string)[len]!='\n');
	*n=size;
	return((int)len);
}

