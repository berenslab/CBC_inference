/* segment emalloc in program nc */

/* Simple standalone memory management */

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
char *malloc(unsigned int size);
//void free(void *ptr);

#ifdef __cplusplus
}
#endif

#define max(x, y)       (((x) < (y)) ? (y) : (x))

#include "ncio.h"

char *zmalloc(unsigned int size);

/*------------------------------------------------*/

char *xxemalloc(unsigned int n)	

/* original malloc */
	           
{
	char *p, *tp;
	int i;

	p = malloc((size_t)n);	

	if (p==(char*)NULL) 
	  ncfprintf(stderr,"emalloc: out of memory.\n");
	else {
	  for (tp=p,i=0; i<n; i++)		/* zero the allocated space */
	   	*tp++ = 0; 
	}
	 return p;

}

/*---------------------------------------------------*/

#define MHASHSIZ1 509             /* prime number to make better hash */
#define MHASHSIZ2 2281   
#define MHASHSIZ3 8273  
#define MHASHSIZ4 32213
#define MHASHSIZ5 129559

#define MHASHSIZ MHASHSIZ3
#define SHASHSIZ MHASHSIZ3

#define MBLKSIZ 2000		/* size of memory blocks for reg heap 0 */
#define SBLKSIZ 2000 		/* size of memory blocks for "save heap" 1 */

#define MAXBLKNUM  2147483000

/* The idea here is to have several heaps and to remember which
   heap a memory space was allocated on.  Then when it is time
   to free the memory, it can be correctly freed from the proper heap.
   This is a simple way to allow similar data types to share blocks,
   allowing for "memory compaction".
*/

#define NHEAP 2

struct memblock;

struct memspace {

        memblock *mblk;
        int size;
};

struct memblock {

        int blknum;
        int nspc;
        int heap;
        int size;
        memblock *hnext;
        char *mpnt;
};

static memblock *mhashtab[NHEAP][MHASHSIZ] = {0};   /* memory hash table */
static memblock *memfree[NHEAP] = {0};		/* list of memblocks to re-use*/
static int mblknum[NHEAP] = {0};		/* unique block number */
static memblock *curblock[NHEAP]={0}; 		/* current block */

void delhashmem (memblock *pnt, int heap);

/*---------------------------------------------------*/


void minithash(void)

{
   int i;

   for (i=0; i<MHASHSIZ; i++) {
     mhashtab[0][i] = (memblock *)NULL;
   }
   for (i=0; i<SHASHSIZ; i++) {
     mhashtab[1][i] = (memblock *)NULL;
   }
   for (i=0; i<NHEAP; i++) {
     mblknum[i] = 0;
     curblock[i] = (memblock *)NULL;
     memfree[i] = (memblock *)NULL;
   }
}

/*---------------------------------------------------*/

int memhash(int num, int heap)

/* make disordered index of memblock num */

{
   unsigned long int i;

   i = num;
   if (i<0) i = -i;
   switch (heap) {
     default:
     case 0: return (i % MHASHSIZ); break;
     case 1: return (i % SHASHSIZ); break;
   }
}

/*---------------------------------------------------*/

void mfree (memblock *newmem)
              
/* Install a memblock in free list.  */
/*  Uses "direct chaining" with field "hnext". */

{
     int i, heap;
     memblock *mpnt,*mlast;

   heap = newmem->heap;
   delhashmem(newmem,heap);			/* remove from hash list */
   newmem->hnext = (memblock*)NULL;
   if (!(mpnt=memfree[heap])) {		/* install first in list ? */
      memfree[heap] = newmem;
   }
   else {				/* otherwise, go to end of list */
      for (; mpnt; mpnt=mpnt->hnext)	/* find last element */
          mlast = mpnt;    
      mlast->hnext = newmem;		/* put new symbol at end of free list */
   }
 /*ncfprintf (stderr,"mfree: freeing unused block * %d.\n",newmem->blknum); /**/
}

/*---------------------------------------------------*/

void minstall (memblock *newmem, int heap)
              
/* Install a memblock in hash table */
/*  Uses "direct chaining" with field "hnext". */

{
    int i;
    memblock *mpnt,*mlast;
    static int minitfl=0;

   if (!newmem) return;

   if (!minitfl) {                      /* initialize table once at start */
        minitfl = 1;
        minithash();
   }
   newmem->hnext = 0;
   i=memhash(newmem->blknum,heap); 	/* initial index into mhashtab */ 
   if (!(mpnt=mhashtab[heap][i])) {     /* install directly in table */
      mhashtab[heap][i] = newmem;
   }
   else {                               /* otherwise, go to end of list */
      for (; mpnt; mpnt=mpnt->hnext)    /* find last element */
          mlast = mpnt;    
      mlast->hnext = newmem;            /* put new symbol at end of list */
   }
}
/*---------------------------------------------------*/

void delhashmem (memblock *pnt, int heap)

/* delete a block of memory from the lookup table */

{
   memblock *mpnt,*mlast;
   int i,found;

  if (!pnt) return;
  
  i = memhash(pnt->blknum, heap);
  mlast = (memblock *)NULL; 
  for (found=0,mpnt = mhashtab[heap][i]; mpnt; mlast=mpnt, mpnt=mpnt->hnext) {
    if (mpnt==pnt) {
       found = 1;
       break;
    }
  }
  if (found) {
     if (mlast) mlast->hnext = mpnt->hnext;     /* delete hash pointer */
     else mhashtab[heap][i] = mpnt->hnext;
  }
}

/*------------------------------------------------*/

void efree(void *pnt)

/* Free a memspace of memory.  When a block of freed memspaces is all */
/*  free, reallocate that block to allow memory re-use. */
/*  No need to zero the old memory, as it will be zeroed in */
/*  emalloc() below. */

{
     int blk;
     memspace *mpnt;
     memblock *mblk;

  mpnt = (memspace *)(((char*)pnt)-sizeof(memspace));
  mblk = mpnt->mblk;

  /* if (mpnt->size > 100) ncfprintf (stderr,"efree %d\n",mpnt->size); /* */
  /* ncfprintf (stderr,"efree: block %d heap %d spc %d\n",
			mblk->blknum,mblk->heap,mblk->nspc); /* */
  if (mblk->nspc > 0) mblk->nspc--;
  else {
     ncfprintf (stderr,"efree: error, can't free unused space.\n");
     return;
  }
  if (mblk!=curblock[mblk->heap])
     if (mblk->nspc==0) mfree (mblk);

  /* ncfprintf (stderr,"efree: block %d free unused memspace %d.\n",
				mblk->blknum,mblk->nspc); */
}

/*------------------------------------------------*/

#define PGRAN (sizeof(char*))

char *zmalloc(unsigned int n, int heap)
	           
{
	int i,found,msize,size,nn;
	char *tp,*x;
	memblock *mpnt,*mlast,*curblk;
	memspace *newspc;

   if (n <= 0) return ((char*)NULL);
   if (heap >= NHEAP) heap - NHEAP-1;
   nn = n/PGRAN;
   if (nn*PGRAN!=n) n = int(nn+1)*PGRAN;
   n += sizeof(memspace);		/* n is size of entire new mspace */
   curblk = curblock[heap];
   if (!curblk || (n > (curblk->size - (curblk->mpnt - ((char*)curblk))))) {
      msize = n+sizeof(memblock);

				     /* search for block that is big enough */
      mlast = (memblock*)NULL;
      for (found=0,mpnt=memfree[heap]; mpnt; mlast=mpnt,mpnt=mpnt->hnext) {
	 if (mpnt->size >= msize) {
		found=1;
		break;
	  }
      }
      if (found) { 				  /* remove from free list */
	  if (mlast) mlast->hnext  = mpnt->hnext;  
	  else       memfree[heap] = mpnt->hnext;
	  curblk = mpnt;   
	  curblk->hnext = (memblock *)NULL;
	  size = mpnt->size;
      /* ncfprintf(stderr,"emalloc: reallocating old block %d heap %d.\n",
				curblk->blknum,heap);/**/
      } 
      else {
	 switch (heap) {
	  default: 
          case 0: size = max(msize,MBLKSIZ); break;
          case 1: size = max(msize,SBLKSIZ); break;
	 }
	 curblk = (memblock*)malloc((size_t)size); /* allocate another */

        /* ncfprintf(stderr,"emalloc: allocating new block of memory.\n"); /* */
        if (!curblk) { ncfprintf(stderr,"emalloc: out of memory.\n");
           return (char *)NULL;
        }
        else {
	  if (++mblknum[heap] > MAXBLKNUM) {
		ncfprintf (stderr,"zmalloc: too many blocks.\n");
		mblknum[heap] = MAXBLKNUM;
	        return ((char *)NULL);
	  }
	  curblk->blknum = mblknum[heap];
	  curblk->nspc = 0;
	  curblk->heap = heap;
	  curblk->size = size;
 	  /* ncfprintf(stderr,"emalloc: allocating new block %d on heap %d.\n",
			mblknum[heap],heap);/**/
        }
     }	/* else not found */ 

     curblk->mpnt = ((char*)curblk) + sizeof(memblock);
     minstall(curblk,heap);	   	/* stow away the block */
   }
   newspc = (memspace *)curblk->mpnt;
   newspc->mblk = curblk;
   newspc->size = n;
   for (i=sizeof(memspace),tp=((char*)newspc)+i; i<n; i++) /* erase */ 
	   	*tp++ = 0; 
   curblk->mpnt += n;
   curblk->nspc++;
   if (curblk->nspc > MAXBLKNUM) curblk->nspc = MAXBLKNUM;
   curblock[heap] = curblk;
   return (((char*)newspc) + sizeof(memspace));
}

/*------------------------------------------------*/

char *emalloc(unsigned int n)

{
   zmalloc(n,0);		/* allocate from regular heap */
}

/*------------------------------------------------*/

char *smalloc(unsigned int n)

{
   zmalloc(n,1);		/* allocate from "save heap" */
}

/*------------------------------------------------*/

char *zzsmalloc(unsigned int n)

/* allocate memory for space that won't be erased much */

/*  Original simpler version. Not used. */
 
{
	int i,size,nn;
	char *ospnt, *tp,*x;
	static int smfree=0;
	static char *spnt=(char*)NULL;

   if (n <= 0) return (char*)NULL;
   nn = n/PGRAN;
   if (nn*PGRAN!=n) n = int(nn+1)*PGRAN;
   if (n > smfree) {			/* allocate another block */
      size = max(n,SBLKSIZ);
      spnt = malloc((size_t)size);	
      /* ncfprintf(stderr,"smalloc: allocating block of memory %d.\n",size); /* */
      if (!spnt) {
	   ncfprintf(stderr,"smalloc: out of memory.\n");
           return (char *)NULL;
      }
      else {
	for (tp=spnt,i=0; i<size; i++)		/* zero the allocated space */
	   	*tp++ = 0; 
        smfree = size;
      }
   }
   ospnt = spnt;
   spnt += n;
   smfree -= n;
   return (ospnt);
}

/*------------------------------------------------*/

void zzsfree(void *pnt)

/* deallocate memory for space that won't be erased much */

/*  Original version. Not used. */
{
return;
}

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
	ncfprintf (stderr,"out of space\n");
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

