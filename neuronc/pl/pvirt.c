/* Segment PVIRT in program mprinttf.c */

/* virtual memory subroutines */

#include <stdio.h>

/*----------------------------------*/

#include <fcntl.h>
#include <unistd.h>

#define BBITS 5			/* bits in one side of square block */
#define BSIZE (1 << (BBITS*2)) 	/* size of blocks in bytes */
#define BLKBUF 256		/* number of blocks in buffer space */

#define MAXMBITS 10		/* largest size arrray for direct mem arr */
#define MMAX  ((1 << MAXMBITS)-1)	/* limit for (x,y) accessing memory */

int bkaddr[BLKBUF] = {0};		/* loc of most recently used blocks */
int bkwrit[BLKBUF] = {0};		/* their write flags */
int bkpri[BLKBUF] = {0};		/* their most recently used pri */

static int lowmask = 0;		/* mask for addressing into blocks */
static int highmask = 0;
static int hbits = 0;		/* number of high bits used for finding block */
static int mmax = 0;		/* max x or y */


/* char mem[BLKBUF][BSIZE] = {0};	/* the old way: virtual buffer */

static char *mem = 0;		/* pointer to allocated memory */
static int memfl = 0;		/* = 1 -> memory has been allocated */
static int vmemfl = 0;		/* = 1 -> pixel array is entirely in memory */
static int filsiz = 0;		/* = size of file (size of array) */
static int filblks = 0;		/* = size of file in blocks */

static char *vnam = 0;
static int fd = 0;		/* the file descriptor */

char *malloc (unsigned int n);
void free ();
char *meminit(int bits, int backgr);
void mwrite (int x, int y, int val);
int mread (int x, int y);


/****************************************/

virtinit(int bits, int backgr)

{
    static int i,j,ptr,blksize;
    static char tempfil[] = "/usr/tmp/virtXXXXXX";
    char *mktemp();

mmax = (1 << bits) - 1;

if (bits <= MAXMBITS) {		/* Pixel array is 1 MB or less */
    vmeminit(bits, backgr);
    vmemfl = 1;
}
else {				/* array is greater than 1 MB */
   vmemfl = 0;

   filsiz = (1 << (bits * 2));
  
  if (! vnam) {
    vnam = mktemp (tempfil);
    fd=open (vnam,O_RDWR | O_CREAT | O_TRUNC,0666);	/* make file */
  }

  blksize = BLKBUF*BSIZE;
  mem = meminit (blksize,backgr);	/* allocate mem and set to "backgr" */
  if (!mem) return;

  filblks = filsiz / BSIZE;
  for (i=0; i<filblks; i++) {
    lseek (fd, i*BSIZE, SEEK_SET);
    write (fd,mem,BSIZE);				/* clear file */
  }

  for (i=0; i<BLKBUF; i++)				/* clear tables */
   {
    bkaddr[i] = -1;
    bkwrit[i] = 0;
    bkpri[i] = BLKBUF - i;			/* set up initial priority */
   }

						/* set up address masks */
  lowmask = ~(~0 << BBITS);			/* 0x001f */
  hbits = bits - BBITS;
  highmask =  ~((~0 << bits) | lowmask);

 }  	/* else virtual memory */

}

/*--------------------------------------*/

virtend()
/* erase file when done */

{
  if (fd) {
      close(fd);
      unlink (vnam);
      fd = 0;
  }
  if (memfl) {
     if (mem) free (mem);
     memfl = 0;
  }
}

/*--------------------------------------*/

int vread (x,y)
      int x,y;

/* Read one byte from virtual memory */
/*  The idea here is to use the lower bits of both
    x and y as the address in the block.  If
    adjacent pixels are accessed sequentially
    then this reduces the number of new block  
    reads.
*/

{
   static int word,block,rblk; 

 if (vmemfl) return (mread(x,y));
 else {
  word  = ((y & lowmask) << BBITS)  + (x & lowmask); /* word is 5 lower bits */

  block=((y & highmask)<<(hbits-BBITS)) + (x >> BBITS);/* block is upper bits*/

  rblk = getblk (block);
  return (mem[rblk*BSIZE+word]);
 }
}

/*--------------------------------------*/
 
vwrite (x,y,val)
    int x,y,val;

/* write one byte to virtual memory */

{
   static int word,block,wblk; 

 if (vmemfl) mwrite(x,y,val);
 else {
  word  = ((y & lowmask) << BBITS)  + (x & lowmask); /* word is 5 lower bits */

  block=((y & highmask)<<(hbits-BBITS)) + (x >> BBITS);/* block is upper bits*/

  wblk = getblk (block);
  bkwrit[wblk] = 1;
  mem[wblk*BSIZE + word] = val;
 }
}

/*--------------------------------------*/

int getblk (req)
      int req;

/* Map the requested block into a smaller number of in-memory blocks.
   Return a memory block that contains the requested block's data.
*/

{
   int i,j,block,lru,lrupri;
   static int found;
   static int lastblk=0;

 for (i=lastblk,j=0,found=0; j<BLKBUF; j++)	/* valid memory image? */
   {
    if (bkaddr[i] == req)
     {
      found = 1;
      break;
     } 
    if (++i >= BLKBUF) i = 0;
   }

 if (found) block = i;				/* found? use this block */

 else   	/* not found */			/* no, then find LRU block */
  {
   for (lrupri=0,i=0; i<BLKBUF; i++)
     if (bkpri[i] > lrupri) {			/* check for LRU */
       lrupri = bkpri[i];
       lru = i;
       if (lrupri == BLKBUF-1) break;
     }
   block = lru;					/* use this block */
   for (i=0; i<BLKBUF; i++)			/* make new priorities */
    bkpri[i]++;
   bkpri[block] = 0;				/* this is now  MRU block */
   if (bkwrit[block]) {				/* write old block */ 
     lseek (fd, bkaddr[block]*BSIZE, SEEK_SET);
     write (fd,mem+block*BSIZE,BSIZE);
   }
   lseek (fd, req*BSIZE, SEEK_SET);
   read (fd,mem+block*BSIZE,BSIZE);
   bkaddr[block] = req;				/* set the new address */
   bkwrit[block] = 0;				/* reset the write flag */
  }

 lastblk = block;
 return (block);
}
 

/****************************************/

				/* Here we're using memory not buffered. */

char *malloc(unsigned int n);
void free (char *p);

static int nbits=0;			/* size of array actually allocated */

static int vsize;

/*--------------------------------------*/

vmeminit(int xbits,int backgr)

{
 nbits = xbits;

 if (nbits > 12) nbits = 12;		/* limit array to 4096 x 4096 */
 if (nbits < 4) nbits = 4;

 vsize = 1 << (nbits * 2);		/* square array */
 mem = meminit(vsize,backgr);
}

/*--------------------------------------*/

int mread (x,y)
      int x,y;

/* Read one byte from virtual memory */


{
   long int byte;

  if (!mem) return 0;
  else {
   if (x > mmax) x = mmax; 
   if (x < 0)    x = 0;
   if (y > mmax) y = mmax; 
   if (y < 0)    y = 0;
   byte  = (y << nbits)  + x; 
   return (mem[byte]);
  }
}

/*--------------------------------------*/
 
void mwrite (x,y,val)
    int x,y,val;

/* write one byte to virtual memory */

{
    long int byte; 


  if (!mem) return;
  else {
   if (x > mmax) x = mmax; 
   if (x < 0)    x = 0;
   if (y > mmax) y = mmax; 
   if (y < 0)    y = 0;
   byte  = (y << nbits)  + x; 
   mem[byte] = val;
  }
}

/*--------------------------------------*/

putarr()

{

   int i,rblk,x,y;
 
  if (vmemfl) {
    for (i=0; i<vsize; i++) 
      putc(mem[i],stdout);
  }
  else {
   if (!fd) return;
   for (y=0; y<=mmax; y++) 
    for (x=0; x<=mmax; x++) 
      putc(vread(x,y),stdout);
  }

}

/*--------------------------------------*/

char *meminit(int msize, int backgr)

{
    static int i,j;
    char *p;

  if ((p = malloc (msize))==NULL) {
     fprintf (stderr,
	"pvirt: meminit: %g MB of memory cannot be allocated.\n", 
		(double)msize/1048576.);
     return NULL;
  }
  else {
    memfl = 1;
    for (i=0; i<msize; i++)
       p[i] = backgr;   
  }

return (p);
}

