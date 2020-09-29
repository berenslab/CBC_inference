/*  Segment MPRINTV in program MONT */

/* System-dependent graphics and text I/O;

   Designed for use with VGA 


	direct access to the board with:
		a) direct I/O to ports with "rtio_inb()" and "rtio_outb()"
		b) use with Interactive 386/ix with VENIX RTX


    See "mprintm.c" for access to VGA through standard VENIX drivers. 
*/

/*  Latest mod  27-July-89	R.G.Smith       */

#include <stdio.h>
#include <fcntl.h>
#include "mdef.h"
#include "vcolors.h"
#include "grdef.h"
/* #include <rtx.h>  */

#define XHHI 639		/* Physical screen coords */
#define XHI 639			/* (640 * 20/27)  Physical screen coords */
#define YHI 479

#define XMAX 79			/* Logical alpha coords */
#define YMAX 24 

#define XSCALE 35
#define YSCALE 35

#define CWIDTH  (.01709 * STDSIZ)	/* char matrix size for enh graph */
#define CHEIGHT (.04 * STDSIZ)

static int xloc,yloc;			/* location of last move */ 
static int xpos,ypos;			/* location of last move for movee */ 
static int cwidth=CWIDTH;
static int cheight=CHEIGHT;

static int xscale = XSCALE;		/* scale for x on screen */
static int yscale = YSCALE;		/* scale for y on screen */
static int xhi = XHI;			/* size for x on screen */
static int yhi = YHI;			/* size for y on screen */
static int x,y;			/* point value for line drawing */

extern int scrcolor;			/* color of screen for enhanced adapt */
extern int scrrev;			/* =1 -> screen reversed vertically */

/********************************************/

putsv(str)
   TEXT *str;

/* Display string on VGA card overlaid on graphics. */

{
 for (; *str; str++)				/* enhanced color monitor */
  putcv(*str);
}

/*-----------------------------------------*/

putcv(ch)
   char ch;

/* output dot matrix char onto VGA screen */

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

 xyconv(xloc/CWIDTH,yloc/CHEIGHT);
}

/************************************/
  
tvinit()

{

  gmodev();
/*  scrcolor = WHITE; */
}

/*--------------------------------------*/

tvexit()

{

 tmodev();
}

/*****************************/

int setvcol (color)
   int color;

{
   static int oldcolor;

  oldcolor = scrcolor;
  scrcolor = color;
  return (oldcolor);
}

/*****************************/

/* xyvidv (xpos,ypos)
  int xpos,ypos;

/* Set text cursor (close) to a graphics position.  */

/* 
{
  xyconv (xpos/CWIDTH,ypos/CHEIGHT);
}
*/

/*****************************/

xyvidv (xpos,ypos)
    int xpos,ypos;

/* set text cursor to a graphics position in GRAPHICS mode */

{
  movev (xpos,ypos);		/* move graphics loc */
  setorg (xpos,ypos);  		/* set the text loc */
}

/*****************************/

/*xyconv (x,y)
   int x,y;

/* X,Y cursor addressing for alphanumerics on console */
/* x = columns, y = rows; (0,0) is at bottom left */

/*{
 x = min(x,XMAX);
 y = min(y,YMAX);
 x = max(x,0);
 y = max(y,0);
 xloc = x * CWIDTH;
 yloc = y * CHEIGHT;
 xyvidv(xloc,yloc);
}
*/

/********************************************/

xyconv (x,y)
   int x,y;

/* Set text cursor to a col, row position in GRAPHICS mode */

{
 x = min(x,XMAX);
 y = min(y,YMAX);
 x = max(x,0);
 y = max(y,0);
 xloc = x * cwidth;
 yloc = y * cheight;
 xyvidv(xloc,yloc);
}

/********************************************/

int tvxsize()

{
  return ((STDSIZ / yhi) * xhi);
}

/********************************************/

int tvysize()

{
  return (STDSIZ);
}

/**************************************************/

vdot (xd,yd)
  int xd,yd;

/* Draw one dot only. */

{
  x = xd / xscale;
  y = yd / yscale;
  if (! scrrev) y = yhi-y;
  pointv();
  xpos = x;
  ypos = y;
}

/*--------------------------------------*/

double xscalev(size)
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

double yscalev(size)
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

movev(x,y)
  int x,y;

/* move to (x,y) routine for VGA.
   square screen size is 480 * 480.
   Scale by 35 (16384/480).
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

drawv(x,y)
  int x,y;

/* draw to (x,y) routine for VGA.
   square screen size is 480 * 480.
   Scale by 35 (16384/480).
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

  linsegv(xpos,ypos,x,y);
  xpos = x;
  ypos = y;
}

/*--------------------------------------*/

/* Driver for VGA card */

/* Draws line segments in color on VGA
*/

#define GRFPG 0x6000			/* graphics page start */
#define GRFPGH (GRFPG / 0x100)		/* high byte of graphics page */

#define DISP_ADDR 0xA0000  /* the physical address of the vga card */
#define DISP_SIZE 0x10000 /* the size of display memory on vga card */
#define GRFSIZE (DISP_SIZE - GRFPG) /* the size of graphics memory used */

#define ARX 0x3c0 /* attribute controller address registers and data
                   registers share same port address- for READS and WRITES */
#define GRX   0x3ce	/* graphics controller index register */
#define GR    0x3cf	/* graphics controller data register */
#define CRX   0x3d4	/* CRT controller index register */
#define SRX   0x3c4	/* sequencer index register */
#define SR    0x3c5	/* sequencer data register */
#define STAT  0x3da	/* input status #1 register controls attribute addr
                 multiplexer, selects attribute address or data to a attr reg */
#define MISC 0x3c2	/* miscellaneous output register */
#define PALL 0x20	/* pallette address source bit in ARX */


#define outb(a,b) rtio_outb(a,b) /* */
#define inb(a)    rtio_inb(a)    /* */

/********************************************/

/* simple low-level driver for VGA */

/*   graphics mode 0x12 register settings */

static char gattrib[] = 
		  {BLACK,BLUE,GREEN,CYAN,RED,MAGENTA,BROWN,WHITE,
		  IWHITE,IBLUE,IGREEN,ICYAN,IRED,IMAGENTA,YELLOW,
		  BWHITE, 0x01,0x00,0x0f,0x00};


static char gseq[] =	  {0,1,0xf,0,6};

static char gcrt[] =	{0x5f,0x4f,0x50,0x82,
		 0x54,0x80,0x0b,0x3e,
	 	 0,0x40,0,0,
	 	 GRFPGH,0,0,0,			/* page zero */
	 	 0xea,0x8c,0xdf,0x28,
	 	 0x00,0xe7,0x04,0xe3,
	 	 0xff};

static char ggrfc[] = {0,0,0,0,
	  	  0,2,5,0xf,
	 	  0xff};

/* text mode (mode 3+) register settings */

static char tattrib[] = {0,1,2,3,4,5,0x14,7,
	 	  0x38,0x39,0x3a,0x3b,
	 	  0x3c,0x3d,0x3e,0x3f,
	 	  0x0c,0,0x0f,8};

static char tseq[] = {0,0,3,0,2};

static char tcrt[] =    {0x5f,0x4f,0x50,0x82,
	 	  0x55,0x81,0xbf,0x1f, 
	 	  0,0x4f,0,0x0e,
	 	  0,0,0,0,
	 	  0x9c,0x8e,0x8f,0x28,
	 	  0x1f,0x96,0xb9,0xa3,
	 	  0xff};

static char tgrfc[] =  {0,0,0,0,
	 	 0,0x10,0x0e,0,
	 	 0xff};


static char *gptr;			/* pointer to graphics page  */
static char *ptr,*vptr;			/* vga pointers  */
static char vgaset=0;			/* 1 => vga mapped already */

/**************************************************/

linsegv (x1,y1,x2,y2)
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
     pointv();
     d = (dy << 1) - dx;
     incr1 = dy << 1;
     incr2 = (dy - dx) << 1;
     while (x < xend)
       {
         x++;
         if (d < 0) d += incr1;
         else     { y += yinc; d += incr2; }
         pointv();
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
      pointv();
      d = (dx << 1) - dy;
      incr1 = dx << 1;
      incr2 = (dx - dy) << 1;
      while (y < yend)
       {
         y++;
         if (d < 0) d += incr1;
         else     { x += xinc; d += incr2; }
         pointv();
       }
   }  /* else */

} /* drawto */


/**************************************************/

/*	Driver for VGA */

vgamap()

{
   extern char *malloc();
					/* alloc virt space */

 if ((vptr=malloc((DISP_SIZE / VMGRAN + 1) * VMGRAN)) == NULL) { 
    fprintf (stderr,"mprintv: can't allocate video memory\n");
    exit(1);
 }
					/* round to mappable boundary */

 vptr = (char *) ((((unsigned int) vptr) + (VMGRAN-1)) & (~(VMGRAN-1)));

 gptr = vptr + GRFPG;			/* pointer to graphics page */

					/* map video memory */

 if (rtphys(vptr,DISP_ADDR,DISP_SIZE,RT_PHYS_RW) == -1) {
    fprintf (stderr,"mprintv: can't map video memory\n");
    exit(1);
  }
  vgaset = 1;
}


/**************************************************/

gmodev()

/* initialize 640 * 350 graphics mode */

{
   static int i;

 if (!setpriv()) return;		/* set I/O privilege */

 if (!vgaset) vgamap();		/* map vga graphics memory into user space */
  
  			/* disable the write prot. on CR 0-7 */
  			/* zero the whole register; will get rewritten */
  outib(CRX,0x11,0);			/* blank screen */


 			/* load all registers for vga graphics mode h12 */

  for (i=0; i<25; i++)			/* program crt controller */
   outib(CRX,i,gcrt[i]);



  for (i=0; i<20; i++)  {		/* setup attributes */
    inb(STAT);
    outb(ARX,i);  	  	/* overwrite PALL to access pallett regs */
    outb(ARX,gattrib[i]); 
  }

  inb(STAT);
  outb(ARX,PALL); 			/* reset PALL to use pallet regs */ 

  for (i=0; i<9; i++)			/* set graphics controller regs */
   outib(GRX,i,ggrfc[i]); 

/* zero seq. reset reg. SR0 to allow proper access to SR and MISC regs*/ 

  for (i=0; i<5; i++) 
   outib(SRX,i,gseq[i]);		/* program sequencer */

  outb(MISC,0xe3);			/* misc register */

  outib(SRX,0,3);			/* start sequencer */

}


/**********************************/

tmodev()

/* return to hi-res alpha mode */

{
   static int i;

 if (!setpriv()) return;
 
  			/* disable the write prot. on CR 0-7 */
  			/* zero the whole register; will get rewritten */
  outib(CRX,0x11,0);			/* blank screen */

 			/* load all registers for vga graphics mode h12 */

  for (i=0; i<25; i++)			/* program crt controller */
   outib(CRX,i,tcrt[i]);


  for (i=0; i<20; i++)  {		/* setup attributes */
    inb(STAT);
    outb(ARX,i);  	  	/* overwrite PALL to access pallett regs */
    outb(ARX,tattrib[i]); 
  }

  inb(STAT);
  outb(ARX,PALL); 			/* reset PALL to use pallet regs */ 

  for (i=0; i<9; i++)			/* set graphics controller regs */
   outib(GRX,i,tgrfc[i]); 

/* zero seq. reset reg. SR0 to allow proper access to SR and MISC regs*/ 

  for (i=0; i<5; i++) 
   outib(SRX,i,tseq[i]);		/* program sequencer */

  outb(MISC,0x67);			/* misc register */

  outib(SRX,0,3);			/* start sequencer */

  endpriv();

}

/**************************************************/

pointv()

/* similar to EGA version */
/* writes one pixel to screen in write mode 2, grphcs mode h13 */

{
   extern int x,y;
   int bytepos;
   register char bitmask;
   char *vptr,dummy;
   extern char *gptr;

bitmask = 0x80 >> x%8;
bytepos = y*80 + x/8;

/* use bitmask to disable all but bit interested in */

outb(GRX,0x08);  			/* select bitmask register*/
outb(GRX+1,bitmask);  			/* select bit plane 2 */

vptr = gptr + bytepos;
dummy   = *vptr; 
*vptr= scrcolor;   		/* color passed implicitly from global var */

}
/**************************************************/
 
clrcv ()

/* clear screen and home up */

{
   long i,j;

 outib(GRX,5,0);			/* set write mode 0 for fast erase */

 outib(GRX,3,0);			/* set to replace mode */

 outib(GRX,8,0xff);			/* set write mask to enable all bits */

 for (i=0; i<GRFSIZE; i++)
   *(gptr+i) = 0;			/* clear video buffer */

 outib(GRX,5,2);			/* set write mode 2 */
 xyconv (0,24);				/* home up cursor */

}

/*clrczzz()   /* blanks the display screen in graphics mode */
/*
{
char *ptr;

/* use bitmask to enable all bits */
/*outb(GRX,0x08);  /* select bitmask register*/
/*outb(GRX+1,0xff);  /* select all bits */

/* clear the entire display memory */
/*for(ptr=gptr; ptr<gptr+GRFSIZE; ptr++)
   {
    *ptr=0x00;
   } 

}
*/
/**************************************************/

