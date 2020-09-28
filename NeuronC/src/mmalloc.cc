/* segment emalloc in program nc */

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
char *malloc(unsigned int size);
void free(void *ptr);

#ifdef __cplusplus
}
#endif

#include "ncio.h"

/*------------------------------------------------*/

char *zmalloc(unsigned int n)	

/* original malloc */

{
	char *p, *tp;
	unsigned int i;

	p = malloc((size_t)n);	

	if (p==(char*)NULL) 
	  fprintf(stderr,"emalloc: out of memory.\n");
	else {
	  for (tp=p,i=0; i<n; i++)		/* zero the allocated space */
	   	*tp++ = 0; 
	}
	 return p;

}

/*---------------------------------------------------*/

void efree(void *pnt)

{ 
  free (pnt);
}

/*------------------------------------------------*/

char *emalloc(unsigned int n)

{
   return (zmalloc(n));		/* allocate from regular heap */
}

/*------------------------------------------------*/

char *smalloc(unsigned int n)

{
   return (zmalloc(n));		/* allocate from "save heap" */
}

/*------------------------------------------------*/

/*
char *zzsmalloc(unsigned int n)
*/

/* allocate memory for space that won't be erased much */

/*  Original simpler version. Not used. */

/* 
{
	int i,size,nn;
	char *ospnt, *tp,*x;
	static int smfree=0;
	static char *spnt=(char*)NULL;

#define max(x, y)       (((x) < (y)) ? (y) : (x))

#define PGRAN (sizeof(char*))
#define SBLKSIZ 2000

   if (n <= 0) return (char*)NULL;
   nn = n/PGRAN;
   if (nn*PGRAN!=n) n = int(nn+1)*PGRAN;
   if (n > smfree) {
      size = max(n,SBLKSIZ);
      spnt = malloc((size_t)size);	
      if (!spnt) {
	   fprintf(stderr,"smalloc: out of memory.\n");
           return (char *)NULL;
      }
      else {
	for (tp=spnt,i=0; i<size; i++)
	   	*tp++ = 0; 
        smfree = size;
      }
   }
   ospnt = spnt;
   spnt += n;
   smfree -= n;
   return (ospnt);
}
*/

/*------------------------------------------------*/

/*
#define MODSIZ 100000000

static char modarr[MODSIZ]={0};
static char *modfree = &modarr[0];

 char *zmalloc(unsigned int len, int h)
*/

/* return pointer to newly allocated space.
   Works like malloc() and emalloc().
*/

/*	version which doesn't use library version of malloc()  */

/*
{
     char *tfree;

    if (modfree+len >= modarr+MODSIZ) {
	fprintf (stderr,"out of space\n");
        return ((char*)NULL);
    }
   else {
    tfree = modfree;
    modfree += len;
    return tfree;
   }
}
*/

/*------------------------------------------------*/

