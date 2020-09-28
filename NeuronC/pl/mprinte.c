/*  Segment MPRINTE in program MONT */

/* System-dependent graphics and text I/O;

   Designed for use with Enhanced Graphics Adapter

   Special code for VENIX sys V rel 2.2 has been commented out.


    This code allows several ways of writing to the EGA screen.

	1) direct access to the board with:
		a) I/O privileges given by setpriv() - console(7);
		b) I/O to ports with inb() and outb();
		c) writes to screen with putesb(3C) and clresn(3C). 

	2)direct access to the board with:
		a) direct I/O to ports with "rtio_inb()" and "rtio_outb()"
		b) use with Interactive 386/ix with VENIX RTX

	3) access through /dev/iomem and putesb(3C)

    Only the direct access method is currently used, because
     it is faster and easier. The other methods are commented out.
    See "mprintm.c" for access to EGA through standard VENIX drivers. 
*/

/*  Latest mod  5-Mar-88	R.G.Smith       */

#include <stdio.h>
#include <fcntl.h>
#include "mdef.h"
#include "colors.h"
#include "grdef.h"
/* #include <rtx.h>  */

#define XHHI 639		/* Physical screen coords */
#define XHI 486			/* (640 * 19/25)  Physical screen coords */
#define YHI 349

#define XMAX 79			/* Logical alpha coords */
#define YMAX 24 

#define XSCALE 36
#define YSCALE 47

#define STDSCL  (STDSIZ/YHI)		/* scale from 350 to STDSIZ */

#define CWIDTH  (6  * STDSCL)		/* char matrix size for enh graph */
#define CHEIGHT (14 * STDSCL)

static int xloc,yloc;			/* location of last move */ 
static int xpos,ypos;			/* location of last move for movee */ 
static int cwidth=CWIDTH;
static int cheight=CHEIGHT;

static int xscale = XSCALE;		/* scale for x on screen */
static int yscale = YSCALE;		/* scale for y on screen */
static int xhi = XHI;			/* size for x on screen */
static int yhi = YHI;			/* size for y on screen */

extern int scrcolor;			/* color of screen for enhanced adapt */
extern int scrrev;			/* =1 -> screen reversed vertically */
static int x,y;			/* point value for line drawing */

/********************************************/

putse(str)
   TEXT *str;

/* Display string on Enhanced graphics card overlaid on graphics. */

{
 for (; *str; str++)				/* enhanced color monitor */
  putce(*str);
}

/*-----------------------------------------*/

putce(ch)
   char ch;

/* output dot matrix char onto 
   Enhanced Graphics screen */

{
  switch (ch)
   {
      case 0x08:
	  xloc -= CWIDTH;
	  break;
      case '\r':
	  break;
      case '\n':
	  xloc = 0;
	  yloc -= CHEIGHT;
	  break;
      default:
     /*     scr_chr(ch); */			/* in vidio.s */
          xloc += CWIDTH;
	  break;
       }      /* switch */

 xycon(xloc/CWIDTH,yloc/CHEIGHT);
}

/************************************/
  
teinit()


#define HIRES 0x10		/* Enhanced graphics high res mode (350x640) */
#define ALPHA 3

{

  gmodee();
/*  scrcolor = WHITE; */
}

/*--------------------------------------*/

teexit()

{

 tmodee();
}

/*****************************/

int setecol (color)
   int color;

{
   static int oldcolor;

  oldcolor = scrcolor;
  scrcolor = color;
  return (oldcolor);
}

/*****************************/

xyvide (xpos,ypos)
  int xpos,ypos;

/* Set text cursor (close) to a graphics position.
   Of course, this can't really be done well with
   IBM enhanced graphics board using standard text fonts. 
   Use GRAPHICS mode if accurate positioning is really needed.
 */

{
  xycone (xpos/CWIDTH,ypos/CHEIGHT);
}

/*****************************/

xyvid (xpos,ypos)
    int xpos,ypos;

/* set text cursor to a graphics position in GRAPHICS mode */

{
  movee (xpos,ypos);		/* move graphics loc */
  setorg (xpos,ypos);  		/* set the text loc */
}

/*****************************/

xycone (x,y)
   int x,y;

/* X,Y cursor addressing for alphanumerics on console */
/* x = columns, y = rows; (0,0) is at bottom left */

{
 x = min(x,XMAX);
 y = min(y,YMAX);
 x = max(x,0);
 y = max(y,0);
 xloc = x * CWIDTH;
 yloc = y * CHEIGHT;
/* scr_xy(x,YMAX-y); */
 xyvid(xloc,yloc);
}

/********************************************/

xycon (x,y)
   int x,y;

/* Set text cursor to a col, row position in GRAPHICS mode */

{
 x = min(x,XMAX);
 y = min(y,YMAX);
 x = max(x,0);
 y = max(y,0);
 xloc = x * cwidth;
 yloc = y * cheight;
 xyvid(xloc,yloc);
}

/********************************************/

int texsize()

{
  return ((STDSIZ / yhi) * xhi);
}

/********************************************/

int teysize()

{
  return (STDSIZ);
}

/**************************************************/

edot (x,y,color)

/* Draw one dot only. Not currently used. */

{
return;
/*setecol(color);
movee (x,y);
drawe (x,y);
*/
}

/*--------------------------------------*/

double xscalee(size)
   double size;

{
   extern int xscale,xhi;
   int temp;

  temp = xscale;
/*  xscale = ((double)XSCALE) / size; */
/*  xhi = (STDSIZ / xscale) - 2; */
  return(((double)temp);
}

/*--------------------------------------*/

double yscalee(size)
   double size;

{
   extern int yscale,yhi;
   int temp;

  temp = yscale;
/*  yscale = ((double)YSCALE) / size; */
/*  yhi = (STDSIZ / yscale) - 2; */
  return(((double)temp);
}

/*--------------------------------------*/

movee(x,y)
  int x,y;

/* move to (x,y) routine for Enhanced Graphics.
   square screen size is 350 * 350.
   Scale by 47 (16384/350).
   Scale x by 36 (19/25 * 47); EGA physical screen is 19/25 cm.
*/

{
     extern int scrrev;

  x /= xscale;
  y /= yscale;

  xpos = x;
  ypos = y;
  if (! scrrev) ypos = yhi-ypos;
}

/*--------------------------------------*/

drawe(x,y)
  int x,y;

/* draw to (x,y) routine for Enhanced Graphics.
   square screen size is 350 * 350.
   Scale by 47 (16384/350).
*/

{
   extern int scrrev;
 
  x /= xscale;
  y /= yscale;
  if (! scrrev) y = yhi-y;

  if (x > XHHI) x = XHHI;
  if (x < 0)   x = 0;
  if (y > yhi) y = yhi;
  if (y < 0)   y = 0;
  if (xpos > XHHI) xpos = XHHI;
  if (xpos < 0)   xpos = 0;
  if (ypos > yhi) ypos = yhi;
  if (ypos < 0)   ypos = 0;

  linseg(xpos,ypos,x,y);
  xpos = x;
  ypos = y;
}

/*--------------------------------------*/

/* Driver for EGA card */

/* Draws line segments in color on Enhanced Graphics Adapter
*/

#define GRF 0x3ce

/* #define VIDBUF 0xa0000		/* page zero */
#define VIDBUF 0xa8000			/* page one  */

#define VIDBUFSIZE 0x8000		/* page size */


#define outb(a,b) rtio_outb(a,b) /* */
#define inb(a)    rtio_inb(a)    /* */

/********************************************/

/* simple low-level driver for Enhanced Graphics Adapter */

/*   graphics mode 0x10 register settings */

static char gattrib[] = {0,1,2,3,4,5,6,7,
	 	  0x38,0x39,0x3a,0x3b,
	 	  0x3c,0x3d,0x3e,0x3f,
	  	  1,0,0x0f,0};


static char gseq[] =	  {1,0xf,0,6};

static char gcrt[] =	{0x5b,0x4f,0x53,0x37,
		 0x52,0x00,0x6c,0x1f,
	 	 0,0,0,0,
	 	 0x80,0,0,0,			/* page one */
	/* 	 0x00,0,0,0,			/* page zero */
	 	 0x5e,0x2b,0x5d,0x28,
	 	 0x0f,0x5f,0x0a,0xe3,
	 	 0xff};

static char ggrfc[] = 	 {0,0,0,0,
	  	  0,0,5,0xf,
	 	  0xff};

/* text mode (mode 3) register settings */

static char tattrib[] = {0,1,2,3,4,5,0x14,7,
	 	  0x38,0x39,0x3a,0x3b,
	 	  0x3c,0x3d,0x3e,0x3f,
	 	  8,0,0x0f,0};

static char tseq[] =  	 {1,3,0,3};

static char tcrt[] =    {0x5b,0x4f,0x53,0x37,
	 	  0x51,0x5b,0x6c,0x1f, 
	 	  0,0x0d,0x0b,0x0c,
	 	  0,0,0,0,
	 	  0x5e,0x2b,0x5d,0x28,
	 	  0x08,0x5f,0x0a,0xa3,
	 	  0xff};

static char tgrfc[] = 	 {0,0,0,0,
	 	 0,0x10,0x0e,0,
	 	 0xff};

static char bittab[] = {0x80,0x40,0x20,0x10,8,4,2,1};



static char *vptr;			/* pointer to video buffer */

/**************************************************/

linseg (x1,y1,x2,y2)
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
     pointe();
     d = (dy << 1) - dx;
     incr1 = dy << 1;
     incr2 = (dy - dx) << 1;
     while (x < xend)
       {
         x++;
         if (d < 0) d += incr1;
         else     { y += yinc; d += incr2; }
         pointe();
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
      pointe();
      d = (dx << 1) - dy;
      incr1 = dx << 1;
      incr2 = (dx - dy) << 1;
      while (y < yend)
       {
         y++;
         if (d < 0) d += incr1;
         else     { x += xinc; d += incr2; }
         pointe();
       }
   }  /* else */

} /* drawto */


/**************************************************/

/*	Driver for enhanced graphics adapter */

gmodee()

/* initialize 640 * 350 graphics mode */

{
   static int i,q;
   extern char *malloc();

 if (!setpriv()) return;		/* set I/O privilege */

					/* alloc virt space into this process */

 if ((vptr=malloc(((VIDBUFSIZE / VMGRAN) + 1) * VMGRAN)) == (char *)-1) { 
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

/* phys(0,0,(unsigned)0x540);		/* set es to 0xa8000 (page one) */
/* phys(0,0,(unsigned)0x500);		/* set es to 0xa0000 (page zero) */

  outib(0x3d4,0,0);			/* blank screen */

  inb(0x03da);				/* reset 3c0 (ARX) to "index" state */

  for (i=0; i<20; i++)  {		/* setup attributes */
   outb(0x3c0,i);
   outb(0x3c0,gattrib[i]); 
  }
  outb(0x3c2,0xa7);			/* misc register */

  outib(0x3c4,0,1);			/* halt sequencer */

  for (i=0; i<4; i++) 
   outib(0x3c4,i+1,gseq[i]);		/* program sequencer */

  outib(0x3c4,0,3);			/* start sequencer */

  for (i=0; i<25; i++)			/* program crt controller */
   outib(0x3d4,i,gcrt[i]);

  outb(0x3cc,0);			/* set graphics board position */

  outb(0x3ca,1);

  for (i=0; i<9; i++)			/* set graphics controller regs */
   outib(GRF,i,ggrfc[i]); 

/* vidclr() */

  inb(0x3da);

  outb(0x3c0,0x20);			/* enable display */

  outib(GRF,5,2);			/* set write mode 2 */

  outib(0x3c4,2,0xf);			/* set map mask */

}

/**************************************************/

tmodee()

/* return to hi-res alpha mode */

{
   static int i;

  if (!setpriv()) return;

  outib(0x3d4,0,0);			/* blank screen */
  inb(0x3da);
  
  for (i=0; i<20; i++)  {		/* setup attributes */
   outb(0x3c0,i);
   outb(0x3c0,tattrib[i]); 
  }

  outb(0x3c2,0xa7);			/* misc register */

  outib(0x3c4,0,1);			/* halt sequencer */

  for (i=0; i<4; i++) 
   outib(0x3c4,i+1,tseq[i]);		/* program sequencer */


  outib(0x3c4,0,3);			/* start sequencer */

  for (i=0; i<25; i++)
   outib(0x3d4,i,tcrt[i]);

  outb(0x3cc,0);			/* set graphics board position */

  outb(0x3ca,1);

  for (i=0; i<9; i++)			/* set graphics controller regs */
   outib(GRF,i,tgrfc[i]); 

  inb(0x3da);

  outb(0x3c0,0x20);

  endpriv();

}

/**************************************************/

pointe()

{
    static int bit,ym,loc;
    static char ch;
    extern char bittab[];

  bit = *(bittab + (x & 7)); 			/* find the bit to be written */
  outb(GRF,8);					/* set the graphics index reg */
  outb(GRF+1,bit);				/* set the bit */
  ym = y << 4;					/* yloc = y * 16 */
  loc = (ym + (ym << 2) + (x >> 3));
						/* (y * 80) + (x / 8)   */
  ch = *(vptr+loc);				/* latch byte */
  *(vptr+loc) = scrcolor;			/* send color byte */
}

/**************************************************/

clrce ()

/* clear screen and home up */

{
   long i,j;

 outib(GRF,5,0);			/* set write mode 0 for fast erase */

 outib(GRF,3,0);			/* set to replace mode */

 outib(GRF,8,0xff);			/* set write mask to enable all bits */

 for (i=0; i<VIDBUFSIZE; i++)
   *(vptr+i) = 0;			/* clear video buffer */

 outib(GRF,5,2);			/* set write mode 2 */
 xycon (0,24);				/* home up cursor */

}

/**************************************************/
