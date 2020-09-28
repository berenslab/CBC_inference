/* segment convarr */

/* subroutines to read and write from large array
*/

#include <sys/types.h>
#include "ncsub.h" 
#include "stim.h" 

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#ifdef __cplusplus
}
#endif

#include "ncio.h"

#if (CONVSIZE > 256)		/* if 32 bit compiler */

#define XMAX 512			/* make 512 x 512 array */
#define YMAX 512
#define NUMCONV 12			/* max number of conv arrays */

int convsize=XMAX;			/* size of convolution array */
static int cxmax=XMAX;			/* individual x and y dimensions */
static int cymax=YMAX;


static double *arrtab[NUMCONV];

#ifdef __cplusplus
extern "C" {
#endif
  void free (void *ptr);
  void *malloc (size_t size);
#ifdef __cplusplus
}
#endif


/*------------------------------------*/

char *xsmalloc(unsigned int n)	/* check return from malloc() */
	           
{
	char *p;
	int i;

	p = (char *)malloc((size_t)n);
	if (p == 0)
		ncfprintf(stderr,"Out of memory.\n");
	else
	for (i=0; i<n; i++)		/* clear array */
	  *(p+i) = 0;
	return p;
}

/*------------------------------------*/
static int initcalled=0;

int initc(int num, int size)
{
   unsigned int asize;
   int arr;
   extern int convsize;

  if (initcalled) return(0);
  else initcalled = 1;

  convsize=cxmax = cymax = size;
  asize = sizeof(double) * cxmax * cymax; 
  if (num >= NUMCONV) num = NUMCONV-1;
  if (num < 0) num = 0;
  for (arr=0; arr<num; arr++)
      if ((arrtab[arr] = (double *)xsmalloc(asize))==NULL) {
         return (0);
      }
  return (arr);   /* return number of arrays made */
}

/*------------------------------------*/

void delc(int num)
{
   int arr;

  initcalled = 0;
  if (num >= NUMCONV) num = NUMCONV-1;
  if (num < 0) num = 0;
  for (arr=0; arr<num; arr++)
      if (arrtab[arr]) free (arrtab[arr]);
}

/*------------------------------------*/

double readc(int arr, int x, int y)
{
   int offset;
   double *conarr;

 if (x < 0) x = 0;
 else if (x >= cxmax) x = cxmax-1;
 if (y < 0) y = 0;
 else if (y >= cymax) y = cymax-1;
 conarr = arrtab[arr];
 offset = y * cxmax + x;
 return ((double)*(conarr+offset)); 
}

/*------------------------------------*/

void writec(int arr, int x, int y, double val)
{
   int offset;
   double *conarr;

 if (x < 0) x = 0;
 else if (x >= cxmax) x = cxmax-1;
 if (y < 0) y = 0;
 else if (y >= cymax) y = cymax-1;
 conarr = arrtab[arr];
 offset = y * cxmax + x;
 *(conarr+offset) = val; 
}

/*------------------------------------*/

void addc(int arr, int x, int y, double val)
{
   int offset;
   double *conarr;

 if (x < 0) x = 0;
 else if (x >= cxmax) x = cxmax-1;
 if (y < 0) y = 0;
 else if (y >= cxmax) y = cxmax-1;
 conarr = arrtab[arr];
 offset = y * cxmax + x;
 *(conarr+offset) += val; 
}


/*------------------------------------*/

void mulc(int arr, int x, int y, double val)
{
   int offset;
   double *conarr;

 if (x < 0) x = 0;
 else if (x >= cxmax) x = cxmax-1;
 if (y < 0) y = 0;
 else if (y >= cxmax) y = cxmax-1;
 conarr = arrtab[arr];
 offset = y * cxmax + x;
 *(conarr+offset) *= val; 
}


/*------------------------------------*/

void copyc(int arr1, int arr2)
                 
/* copy arr1 into arr2 */

{
   double *conarr1,*conarr2,*arrend;

   conarr1 = arrtab[arr1];
   conarr2 = arrtab[arr2];
   arrend = conarr1 + (cxmax*cymax);
   for ( ; conarr1<arrend; ) 
       *conarr2++ = *conarr1++; 
}

/*------------------------------------*/

void setc(int arr, double val)

/* set arr to val */

{
   double *conarr,*arrend;

   conarr = arrtab[arr];
   arrend = conarr + (cxmax*cymax);
   for ( ; conarr<arrend; ) *conarr++ = val; 
}

/*------------------------------------*/

#else		/* 16 bit compiler, not 32bit */


#define NUMCONV 10			/* max number of conv arrays */
#define CONVSIZ 8			/* number of sub-arrays per conv arr */

#define YMASK 0x1f 
#define XMAX 256			/* make 256 x 256 array */
#define YMAX 256

static int cxmax=XMAX;
static int cymax=YMAX;

double *arrtab[NUMCONV][CONVSIZ];

/*------------------------------------*/

char *xsmalloc(unsigned int n)	/* check return from malloc() */
	           
{
	char *p;
	int i;

	p = (char *)malloc((size_t)n);
	if (p == 0)
		ncfprintf(stderr,"out of memory\n");
	else
	for (i=0; i<n; i++)	/* clear array */
	  *(p+i) = 0;
	return p;
}

/*------------------------------------*/

void initc(int num, int size)
{
   unsigned int asize;
   int arr,i;

  size = 32768; 
  asize = size; 
  if (num >= NUMCONV) num = NUMCONV-1;
  if (num < 0) num = 0;
  for (arr=0; arr<num; arr++)
    for (i=0; i<CONVSIZ; i++)
      if ((arrtab[arr][i] = (double *)xsmalloc(asize))==NULL) {
         return (0);
      }
  return (arr);		/* return number of arrays made */
}

/*------------------------------------*/

void delc(int num)
{
   int arr,i;

  if (num >= NUMCONV) num = NUMCONV-1;
  if (num < 0) num = 0;
  for (arr=0; arr<num; arr++)
    for (i=0; i<CONVSIZ; i++)
      if (arrtab[arr][i]) free (arrtab[arr][i]);
}

/*------------------------------------*/

double readc(int arr, int x, int y)
{
   int offset,subarr;
   double *conarr;

 if (x < 0) x = 0;
 else if (x >= cxmax) x = cxmax-1;
 if (y < 0) y = 0;
 else if (y >= cymax) y = cymax-1;
 offset = ((y & YMASK) << 8) | x;
 subarr = y >> 5;
 conarr = arrtab[arr][subarr];
 return ((double)*(conarr+offset)); 
}

/*------------------------------------*/

void writec(int arr, int x, int y, double val)
{
   int offset,subarr;
   double *conarr;

 if (x < 0) x = 0;
 else if (x >= cxmax) x = cxmax-1;
 if (y < 0) y = 0;
 else if (y >= cymax) y = cymax-1;
 offset = ((y & YMASK) << 8) | x;
 subarr = y >> 5;
 conarr = arrtab[arr][subarr];
 *(conarr+offset) = val; 
}

/*------------------------------------*/

void addc(int arr, int x, int y, double val)
{
   int offset,subarr;
   double *conarr;

 if (x < 0) x = 0;
 else if (x >= cxmax) x = cxmax-1;
 if (y < 0) y = 0;
 else if (y >= cxmax) y = cxmax-1;
 offset = ((y & YMASK) << 8) | x;
 subarr = y >> 5;
 conarr = arrtab[arr][subarr];
 *(conarr+offset) += val; 
}


/*------------------------------------*/

void mulc(int arr, int x, int y, double val)
{
   int offset,subarr;
   double *conarr;

 if (x < 0) x = 0;
 else if (x >= cxmax) x = cxmax-1;
 if (y < 0) y = 0;
 else if (y >= cxmax) y = cxmax-1;
 offset = ((y & YMASK) << 8) | x;
 subarr = y >> 5;
 conarr = arrtab[arr][subarr];
 *(conarr+offset) *= val; 
}


/*------------------------------------*/

copyc(int arr1, int arr2)
                 

/* copy arr1 into arr2 */

{
   int x,y;

   for (x=0; x<cxmax; x++)
     for (y=0; y<cymax; y++)
        writec (arr2, x, y, readc(arr1, x, y)); 		/* */
}

/*------------------------------------*/

void setc(int arr, double val)

/* set arr to val */

{
   int x,y;

   for (x=0; x<cxmax; x++)
     for (y=0; y<cymax; y++)
        writec (arr, x, y, val); 		/* */
}

/*------------------------------------*/

#endif			/* end of not BIT32 */

