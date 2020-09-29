/* Segment HVIRT in program hide */

/* virtual memory subroutines */

#include <stdio.h>
#include "mdef.h"

/*----------------------------------*/

#ifndef HIDE9


#define VSIZE 256 		/* size of buffer in words */

#ifdef HIDE10
#define BLKBUF 60		/* number of blocks in buffer space */
#define FILSIZ 256				/* # of blocks/file */
#endif

#ifdef HIDE11
#define BLKBUF 40 		/* number of blocks in buffer space */
#define FILSIZ 1024 				/* # of blocks/file */
#endif

#ifdef HIDE12
#define BLKBUF 100 		/* number of blocks in buffer space */
#define FILSIZ 4096 				/* # of blocks/file */
#endif

short int virt[BLKBUF][VSIZE] = {0};		/* virtual buffer */

int bkaddr[BLKBUF] = {0};		/* loc of most recently used blocks */
int bkwrit[BLKBUF] = {0};		/* their write flags */
int bkpri[BLKBUF] = {0};		/* their most recently used pri */

pfil virtfil = {0};		/* virtual memory file */
static TEXT *vnam = 0;

static short int masktab[16] = {1,2,4,8,0x10,0x20,0x40,0x80,0x100,0x200,
			0x400,0x800,0x1000,0x2000,0x4000,0x8000}; 

/****************************************/

virtinit()

begin
    static int i,j,k;
    static TEXT template[] = {"/usr/tmp/virtXXXXXX"};
    extern TEXT *vnam;
    char *mktemp();

  if (! vnam)
   begin
    vnam = mktemp (template);
    k=prvcreate (vnam,0,0,&virtfil,DATFIL,RW);		/* make file */
   end

  for (i=0; i<BLKBUF; i++)
   for (j=0; j<VSIZE; j++)				/* clear buffer */
     virt[i][j] = 0;   

  for (i=0; i<FILSIZ; i++)
    pwrite (&virtfil,i,virt,1);				/* clear file */

  for (i=0; i<BLKBUF; i++)				/* clear tables */
   begin
    bkaddr[i] = -1;
    bkwrit[i] = 0;
    bkpri[i] = BLKBUF - i;			/* set up initial priority */
   end

end

/****************************************/

virtend()
/* erase file when done */

begin
   extern TEXT *vnam;

  prclose(&virtfil);
  unlink (vnam);
end

/****************************************/

int vread (x,y)
      register int x,y;

/* read one bit from virtual memory */


begin
   static int bit,word,block,rblk; 

  bit   =   y & 0X00F;
  word  = ((y & 0X030) << 2)  + (x & 0X03F); 

#ifdef HIDE10
  block = ((y & 0X3C0) >> 2) + ((x & 0X3C0) >> 6);
#endif
#ifdef HIDE11
  block = ((y & 0X7C0) >> 1) + ((x & 0X7C0) >> 6);
#endif
#ifdef HIDE12
  block = ((y & 0XFC0) ) + ((x & 0XFC0) >> 6);
#endif

  rblk = getblk (block);
  return (virt[rblk][word] & masktab[bit] ? 1 : 0);
end

/****************************************/
 
vwrite (x,y)
    register int x,y;

/* write one bit to virtual memory */

begin
    static int bit,word,block,wblk; 

  bit   =   y & 0X00F;
  word  = ((y & 0X030) << 2)  + (x & 0X03F); 

#ifdef HIDE10
  block = ((y & 0X3C0) >> 2) + ((x & 0X3C0) >> 6);
#endif
#ifdef HIDE11
  block = ((y & 0X7C0) >> 1) + ((x & 0X7C0) >> 6);
#endif
#ifdef HIDE12
  block = ((y & 0XFC0) ) + ((x & 0XFC0) >> 6);
#endif

  wblk = getblk (block);
  bkwrit[wblk] = T;
  virt[wblk][word] |= masktab[bit];
end

/****************************************/

vline (x,y1,y2)
      int x,y1,y2;

/* write a vertical line to virtual memory */

begin
   static int top,bot;
   static int lo,hi,y,xwrd,xblk;
   static int word,block,wblk;

 top = max (y1,y2);
 bot = min (y1,y2);
 y = bot;
 if (top-bot >= 32)
  begin
   lo = (bot + 16) & 0XFFF0;		/* find low 16 bit boundary */
   hi = top & 0XFFF0;			/* find high 16 bit boundary */
   for (; y<lo; y++)
    vwrite (x,y);
    xwrd = x & 0X03F;

#ifdef HIDE10
   xblk = (x & 0X3C0) >> 6;
#endif
#ifdef HIDE11
   xblk = (x & 0X7C0) >> 6;
#endif
#ifdef HIDE12
   xblk = (x & 0XFC0) >> 6;
#endif

   for (y=lo; y<hi; y+=16)		/* write whole words */
    begin
     word  = ((y & 0X030) << 2) + xwrd;

#ifdef HIDE10
     block = ((y & 0X3C0) >> 2) + xblk;
#endif
#ifdef HIDE11
     block = ((y & 0X7C0) >> 1) + xblk;
#endif
#ifdef HIDE12
     block = ((y & 0XFC0) ) + xblk;
#endif

     wblk = getblk (block);
     bkwrit[wblk] = T;
     virt[wblk][word] = 0XFFFF;		/* set all bits in word */
    end
  end     /* if */

  for ( ; y<=top; y++)			/* set any extra bits on top */
  vwrite (x,y);

end

/****************************************/

int getblk (req)
      int req;

begin
   int i,j,block,lru,lrupri;
   static BOOL found;
   static int lastblk=0;

 for (i=lastblk,j=0,found=F; j<BLKBUF; j++)	/* valid memory image? */
   begin
    if (bkaddr[i] == req)
     begin
      found = T;
      break;
     end 
    if (++i >= BLKBUF) i = 0;
   end

 if (found) block = i;				/* found? use this block */

 else   	/* not found */			/* no, then find LRU block */
  begin
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
   if (bkwrit[block])				/* write old block */ 
     pwrite (&virtfil,bkaddr[block],virt[block],1);
   pread (&virtfil,req,virt[block],1);			/* Read the new block */
   bkaddr[block] = req;				/* set the new address */
   bkwrit[block] = F;				/* reset the write flag */
  end

 lastblk = block;
 return (block);
end

 
#endif			/* HIDE10 or HIDE11 (not HIDE9) */

/****************************************/


#ifdef HIDE9

#define HIDSIZ 512
#define VSIZE 16384

short virt[VSIZE] = {0};	/* virtual buffer */

static short masktab[16] = {1,2,4,8,0x10,0x20,0x40,0x80,0x100,0x200,0x400,
		         0x800,0x1000,0x2000,0x4000,0x8000}; 

/****************************************/

virtinit()

#define FILSIZ 256 				/* # of blocks/file */

begin
    static int i,j,k;


  for (i=0; i<VSIZE; i++)
     virt[i] = 0;   

end

/****************************************/

virtend()

begin
end

/****************************************/

int vread (x,y)
      register int x,y;

/* read one bit from virtual memory */


begin
   static int bit,word,block,rblk; 

  bit   =   y & 0X00F;
  word  = ((y & 0X3F0) << 5)  + x; 
  return (virt[word] & masktab[bit] ? 1 : 0);
end

/****************************************/
 
vwrite (x,y)
    register int x,y;

/* write one bit to virtual memory */

begin
    static int bit,word,block,wblk; 

  bit   =   y & 0X00F;
  word  = ((y & 0X3F0) << 5)  + x; 
  virt[word] |= masktab[bit];
end

/****************************************/

vline (x,y1,y2)
      int x,y1,y2;

/* write a vertical line to virtual memory */

begin
   static int top,bot;
   static int lo,hi,y,xwrd,xblk;
   static int word,block,wblk;

 top = max (y1,y2);
 bot = min (y1,y2);
 y = bot;
 if (top-bot >= 32)
  begin
   lo = (bot + 16) & 0XFFF0;		/* find low 16 bit boundary */
   hi = top & 0XFFF0;			/* find high 16 bit boundary */
   for (; y<lo; y++)
    vwrite (x,y);
   for (y=lo; y<hi; y+=16)		/* write whole words */
    begin
     word  = ((y & 0X3F0) << 5) + x;
     virt[word] = 0XFFFF;		/* set all bits in word */
    end
  end     /* if */

  for ( ; y<=top; y++)			/* set any extra bits on top */
  vwrite (x,y);

end

#endif		/* HIDE9 */

/****************************************/

