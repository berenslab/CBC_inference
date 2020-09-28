/*  Segment MPRINTH in program MONT */

/* System-dependent graphics and text I/O;

   For use with Hercules type graphics cards where text
   can't be displayed at the same time as graphics.

*/

/*  Latest mod  22-Feb-87	R.G.Smith       */

#include <stdio.h>
#include <fcntl.h>
/* #include <rtx.h> */

#include "mdef.h"
#include "colors.h"
#include "grdef.h"

#define XHI 479
#define XHHI 719 
#define YHI 347

#define XMAX 79
#define YMAX 24 

#define XSCALE 31
#define YSCALE 47

#define STDSCL (STDSIZ/YHI) 		/* scale from 348 to STDSIZ */

#define CWIDTH  (6  * STDSCL)		/* char box size for enh graph */
#define CHEIGHT (14 * STDSCL)

extern int scrcolor;
extern int scrrev;			/* =1 -> video screen reversed */

static int xpos,ypos;			/* location of last move for moveh */ 
static double cwidth = CWIDTH;
static double cheight = CHEIGHT;

static int xscale = XSCALE;		/* scale for x on screen */
static int yscale = YSCALE;		/* scale for y on screen */
static int xhi = XHI;			/* scale for x on screen */
static int yhi = YHI;			/* scale for y on screen */
static int x,y;			/* point values for line drawing */

/********************************************/

thinit()

{
   unsigned int i;

  gmodeh(); 
/*  sethcol(1); */
end

/*****************************/

thexit()

{
  tmodeh();
end

/*****************************/

int sethcol (color)
   int color;

{
   static int oldcolor;

  oldcolor = scrcolor;
  scrcolor = color;
  return (oldcolor);
end


/********************************************/

xyconh (x,y)
   int x,y;

/* X,Y cursor addressing for console */
/* x = columns, y = rows; (0,0) is at bottom left */

{
 x = min(x,XMAX);
 y = min(y,YMAX);
 x = max(x,0);
 y = max(y,0);
 xyvidh ((int)(x * cwidth), (int)(y * cheight));
end

/********************************************/

xyvidh(x,y)
   int x,y;

/* X,Y cursor addressing for alphanumerics */
/* x,y are graphics coords (STDSIZ, STDSIZ) */

{
  moveh (x,y);
  setorg (x,y);
}

/********************************************/

hdot(x,y,color)
   int x,y,color;
{
 return;
/*  moveh (x,y);
  drawh (x,y);
  setorg (x,y);
*/
}

/********************************************/

double xscalh(size)
   double size;

{
   extern int xscale,xhi;
   int temp;

  temp = xscale;
/*  xscale = ((double)XSCALE) / size; */
/*  xhi = (STDSIZ / xscale) - 2; */
  return(((double)XSCALE) / temp);
}

/********************************************/

double yscalh(size)
   double size;

{
   extern int yscale,yhi;
   int temp;

  temp = yscale;
/*  yscale = ((double)YSCALE) / size; */
/*  yhi = (STDSIZ / yscale) - 2; */
  return(((double)YSCALE) / temp);
}

/********************************************/

int thxsize()

{
 return ((STDSIZ / yhi)*xhi);
}

/********************************************/

int thysize()

{
 return (STDSIZ);
}

/********************************************/

moveh(x,y)
  int x,y;

/* Device independent routine for Hercules board.
   Screen size is 348 * 480.
   Physical resolution is 720 in X. This must be multiplied by 3/2 to make 
   square pixels; this gives an overall X resolution of 480.
   Scale by 31 ((2/3)16384/348) in X, 47 (16384/348) in Y.
*/

{
  xpos = x / xscale;
  ypos = y / yscale;
  if (! scrrev) ypos = yhi - ypos; 
}

/*--------------------------------------*/

drawh(x,y)
  int x,y;

/* Device independent routine for Hercules board.
   Screen size is 348 * 480.
   Physical resolution is 720 in X. This must be multiplied by 3/2 to make 
   square pixels; this gives an overall X resolution of 480.
   Scale by 31 ((2/3)16384/348) in X, 47 (16384/348) in Y.
*/

{
   static int lindat[5] = {0};
   static int *linpnt = lindat;

  x /= xscale;
  y /= yscale;
  if (! scrrev) y = yhi - y; 

  if (x > XHHI) x = XHHI;
  if (x < 0)   x = 0;
  if (y > yhi) y = yhi;
  if (y < 0)   y = 0;
  if (xpos > XHHI) xpos = XHHI;
  if (xpos < 0)   xpos = 0;
  if (ypos > yhi) ypos = yhi;
  if (ypos < 0)   ypos = 0;

  linsegh(xpos,ypos,x,y);
  xpos = x;
  ypos = y;
}

/**************************************************/

/* Driver for Hercules card */

/* Draws line segments in color on Hercules Hi-Res board 
*/



#define outb(a,b) rtio_outb(a,b) /* */
#define inb(a)    rtio_inb(a)    /* */

/********************************************/

/*   Hercules hi-res graphics board driver */

/*   Simple driver to start Hi-Res (720 * 348) mode */

#define VIDBUF 0xb8000		/* start of video buffer seg */
#define INDEX  0x3b4		/* index register */
#define CNTRL  0x3b8		/* control register */
#define DSTAT  0x3ba		/* display status register */
#define SOFTSW 0x3bf		/* softw switch port,0=text,1=page1,3=both g */
#define VERTRET 0x80		/* vertical retrace bit in display status reg */

#define VIDBUFSIZE 0x8000 	/* size of one page in graphics screen */

#define SCRNON 8
#define GRPH   0x82
#define TEXT   0x20

static char gtable[12] = {		/* contents for 6845 registers */
			0x35,0x2d,0x2e,0x07,
			0x5b,0x02,0x57,0x57,
			0x02,0x03,0x00,0x00 };

static char ttable[12] = {0x61,0x50,0x52,0x0f,
	 		0x19,0x06,0x19,0x19,
	 		0x02,0x0d,0x0b,0x0c };


static char *vptr;			/* pointer to video buffer */
extern char *malloc();			/* memory allocation from system */

/********************************************/

linsegh (x1,y1,x2,y2)
   int x1,y1,x2,y2;
   
/* Bresenham's Algorithm;
   see Newman and Sproul, 2nd Ed.  for explanation.

 */

{
        static int deltax,deltay,dx,dy;
        static int incr1,incr2,d,xend,yend;
        static int xinc,yinc;

 deltax = x2 - x1;
 deltay = y2 - y1;
 dx = abs(deltax); dy = abs(deltay);

 if (dy <= dx)
   {
     if (deltax < 0)
       {
         yinc = ((deltay < 0) ? 1 : -1); 
         x = x2; y = y2; xend = x1;
       }
     else   /* deltax >= 0 */
       {
         yinc = ((deltay < 0) ? -1 : 1);
         x = x1; y = y1; xend = x2;
       }
     pointh();
     d = (dy << 1) - dx;
     incr1 = dy << 1;
     incr2 = (dy - dx) << 1;
     while (x < xend)
       {
         x++;
         if (d < 0) d += incr1;
         else     { y += yinc; d += incr2; }
         pointh();
       }
    } 
  else   /* dy > dx */
    {
      if (deltay < 0)
        {
          xinc = ((deltax < 0) ? 1 : -1);
          x = x2; y = y2; yend = y1;
        }
      else
        {
          xinc = ((deltax < 0) ? -1 : 1);
          x = x1; y = y1; yend = y2;
        }
      pointh();
      d = (dx << 1) - dy;
      incr1 = dx << 1;
      incr2 = (dx - dy) << 1;
      while (y < yend)
       {
         y++;
         if (d < 0) d += incr1;
         else     { x += xinc; d += incr2; }
         pointh();
       }
   }  /* else */

} /* linsegh */

/********************************************/

gmodeh()

{
   static int i;
   
  if (!setpriv()) return;
					/* alloc virt space */

 if ((vptr=malloc((VIDBUFSIZE / VMGRAN + 1) * VMGRAN)) == (char *)-1) { 
    fprintf (stderr,"mprinte: can't allocate video memory\n");
    exit(1);
 }
					/* round to mappable boundary */

 vptr = (char *) ((((unsigned long) vptr) + (VMGRAN-1)) & (~(VMGRAN-1)));

					/* map video memory */

 if (rtphys(vptr,VIDBUF,VIDBUFSIZE,RT_PHYS_RW) == -1) {
    fprintf (stderr,"mprinte: can't map video memory\n");
    exit(1);
  }

/*  phys(0,0,(unsigned)0x5c0);	/* set 0xb8000 as extra segment page */

  outb(SOFTSW,3);		/* allow graphics mode */

  outb(CNTRL,GRPH);		/* set to graphics mode with screen off */

  for(i=0; i<12; i++)		/* set params for 720 x 348 graphics mode */
   outib(INDEX,i,gtable[i]);

  outb(CNTRL,GRPH+SCRNON);	/* screen on, page 0 */
}

/********************************************/

tmodeh()

{
   static int i;

  if (!setpriv()) return;

  outb(CNTRL,TEXT);		/* change back to text mode with screen off */

  for(i=0; i<12; i++)		/* set params for text mode */
   outib(INDEX,i,ttable[i]);

  outb(CNTRL,TEXT+SCRNON);	/* screen on, page 0 */

  endpriv();
}

/**************************************************/

clrch()

/* erase screen on Hercules board */

{
  register long i,j;

 for (i=0; i<VIDBUFSIZE; i++)
  *(vptr+i) = 0;			/* send color byte */
 xyconh (0,24);				/* home up cursor */
}

/**************************************************/

pointh ()

{
  register int bit;
  char val;
  int loc;
  static char bittab[] = {1,2,4,8,0x10,0x20,0x40,0x80};

  loc = (((y&3) << 13) + 90 * (y >> 2) + (x>>3));
  bit = 7-(x&7);
  val = *(vptr+loc);					/* latch byte */
  if (scrcolor & 0x80)	      val ^=  bittab[bit];
  else if (scrcolor & 0x0f)   val |=  bittab[bit];
  else			      val &= ~bittab[bit];
  *(vptr+loc) = val;

}

/**************************************************/
