/*  Segment MPRINTTF in program MONT */

/* System-dependent graphics and text I/O;

   Designed for use with HP LaserJet II printer.
*/

/*  Latest mod  5-Mar-88	R.G.Smith       */

#include <stdio.h>
#include <fcntl.h>
#include "colors.h"
#include "mdef.h"

#include <signal.h>
void (*istat)();		/* Used by interrupt routine setup below */
				/* see "UNIX Programming, second edition" */

/* typedef int __sighandler_t;	/* */

#define STDSIZ 16384

#define ESC 033
#define FF  014

#define RASRES  280		/* resolution for laserjet (lines/inch) */
/* #define XRAS (8*RASRES)		/* x raster size for array */
/* #define YRAS (10*RASRES)	/* y raster size for array */

#define RASBITS 10		/* basic resolution of array */
#define XRAS (1<<RASBITS)	/* x raster size for array */
#define YRAS (1<<RASBITS) 	/* y raster size for array */

#define XHHI (XRAS-1)		/* Physical screen coords */
#define XHI XRAS		/* Physical screen coords */
#define YHI YRAS

#define XMAX 79			/* Logical alpha coords */
#define YMAX 24 

#define XSCALE ((STDSIZ/XRAS)+1)		/* STDSIZ/XRAS (16384/1200) */
#define YSCALE XSCALE

#define STDSCL  (STDSIZ/YHI)		/* scale from 1500 to STDSIZ */

#define CWIDTH  (6  * STDSCL)		/* char matrix size for enh graph */
#define CHEIGHT (14 * STDSCL)

#define pltout stdout

static int xloc,yloc;			/* location of last move */ 
static int xpos,ypos;			/* location of last move for movee */ 
static int cwidth=CWIDTH;
static int cheight=CHEIGHT;

static int xscale = XSCALE;		/* scale for x on screen */
static int yscale = YSCALE;		/* scale for y on screen */
static int xhi = XHI;			/* size for x on screen */
static int yhi = YHI;			/* size for y on screen */

static int x,y;
static int grflag=0;		/*  =1 -> graphics have used draw() */
static int gmodefl=0;		/*  =1 -> gmode has been called */
static int psize = RASBITS;		/* 10 bits -> 1024 x 1024 (1 MB) array */

extern int scrcolor;			/* color of screen for enhanced adapt */
extern int scrrev;			/* =1 -> screen reversed vertically */
extern int backgr;
extern int bitres;			/* size of pixel array from "vid.c" */

int vread (int x, int y);
void vwrite (int x, int y, int val);
void virtinit(int bits, int backgr);
void virtend();
void putcf(char ch);
void xyvidf (int xpos, int ypos);

/********************************************/

onintr(void)
{
 virtend();
 exit();
}

/**************************************************/

putsf(char *str)
             

/* Display string on graphics. */

{
 for (; *str; str++)				/* enhanced color monitor */
  putcf(*str);
}

/*-----------------------------------------*/

void putcf(char ch)
           

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

 xyconf(xloc/CWIDTH,yloc/CHEIGHT);
}

/************************************/
  
tfinit(void)
{

if (bitres==0) bitres = psize;
if (! gmodefl)			/* if this routine hasn't already been called */
  virtinit(bitres,backgr);	/* initialize the virtual raster array */

xscale = STDSIZ / (1<<bitres);
yscale = xscale;
xhi = yhi = (1 << bitres);

istat = signal (SIGINT, SIG_IGN);	/* save orig staus */
if (istat != SIG_IGN)
  (void) signal (SIGINT,(__sighandler_t)onintr);   /* go to "onintr" at ^C */
 grflag = 0;
 gmodefl = 1;
/* scrcolor = WHITE; */
}

/*--------------------------------------*/

tfexit(void)
{
  virtend();
}

/*--------------------------------------*/


grdumpf(void)

/* output graphics to printer after tfexit() */

{

if (grflag) {				/* if graphics were used */

/* 
fprintf (pltout,"} def\n");
fprintf (pltout,"1 1 scale\n");
fprintf (pltout,"20 20 translate\n");
fprintf (pltout,"%d %d %g bitdump\n",XRAS,YRAS,0.27);
*/

  frastout();				/* send raster data to printer */

/*
fprintf (pltout,"\n");
fprintf (pltout,"showpage\n");
fprintf (pltout,"%%%%Trailer\n");
*/

  grflag=0;
}

if (gmodefl) virtend();	 		/* erase virtual raster file */

gmodefl = 0;				/* reset graphics mode flag */

}

/*****************************/

int setfcol (int color)
{
   static int oldcolor;

  oldcolor = scrcolor;
  scrcolor = color;
  return (oldcolor);
}

/*****************************/

void xyvidf (int xpos, int ypos)

/* set text cursor to a graphics position in GRAPHICS mode */

{
  movef (xpos,ypos);		/* move graphics loc */
  setorg (xpos,ypos);  		/* set the text loc */
}

/*****************************/

xyconf (int x, int y)
           
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
 xyvidf(xloc,yloc);
}

/********************************************/

int tfxsize(void)
{
/*  return ((STDSIZ / yhi) * xhi); */
  return (STDSIZ);
}

/********************************************/

int tfysize(void)
{
/*  return (STDSIZ); */
  return (STDSIZ * 1.25);

}

/**************************************************/

fdot (int x, int y, int color)

/* Draw one dot only. Not currently used. */

{
return;
/*setfcol(color);
movef (x,y);
drawf (x,y);
*/
}

/*--------------------------------------*/

double xscalf(double size)
{
   extern int xhi;
   int temp;

  temp = xscale;
/*  xscale = ((double)XSCALE) / size; */
/*  xhi = (STDSIZ / xscale) - 2; */
  return((double)temp);
}

/*--------------------------------------*/

double yscalf(double size)
{
   extern int yhi;
   int temp;

  temp = yscale;
  yscale = ((double)YSCALE) / size;
/*  yhi = (STDSIZ / yscale) - 2; */
  return((double)temp);
}

/*--------------------------------------*/

movef(int x, int y)

/* move to (x,y) routine for Laserjet II.
   Square screen size is XRAS * XRAS (1200).
   Scale x by 14 (16384/1200).
   Scale y by 14. 
*/

{
     extern int scrrev;

  x /= xscale;
  y /= yscale;

  xpos = x;
  ypos = y;
/*  if (! scrrev) ypos = yhi-ypos; */
}

/*--------------------------------------*/

drawf(int x, int y)

/* draw to (x,y) routine for Laserjet II.
   Square screen size is XRAS * XRAS (1200).
   Scale x by 14 (16384/1200).
   Scale y by 14. 
*/

{
   extern int scrrev;

  x /= xscale;
  y /= yscale;
/*  if (! scrrev) y = yhi-y; */

  if (x > XHHI) x = XHHI;
  if (x < 0)   x = 0;
  if (y > yhi) y = yhi;
  if (y < 0)   y = 0;
  if (xpos > XHHI) xpos = XHHI;
  if (xpos < 0)   xpos = 0;
  if (ypos > yhi) ypos = yhi;
  if (ypos < 0)   ypos = 0;

  linsegf(xpos,ypos,x,y);
  xpos = x;
  ypos = y;
  grflag = 1;			/* set graphics flag when graphics used */
}

/**************************************************/

linsegf (int x1, int y1, int x2, int y2)
                   
   
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
     pointf();
     d = (dy << 1) - dx;
     incr1 = dy << 1;
     incr2 = (dy - dx) << 1;
     while (x < xend)
       {
         x++;
         if (d < 0) d += incr1;
         else     { y += yinc; d += incr2; }
         pointf();
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
      pointf();
      d = (dx << 1) - dy;
      incr1 = dx << 1;
      incr2 = (dx - dy) << 1;
      while (y < yend)
       {
         y++;
         if (d < 0) d += incr1;
         else     { x += xinc; d += incr2; }
         pointf();
       }
   }  /* else */

} /* drawto */

/**************************************************/

gmodef(void)

/* set graphics mode */

{
}

/**************************************************/

tmodef(void)

/* return to text mode */

{
}

/**************************************************/

frastout(void)

/* Put out a raster of color bytes from the page buffer.
*/

{
   int i,j,pixel,size;

 putarr(); 


/*  size = (1 << bitres);

  for (j=0; j<size; j++) {			/* for all rows */
/*    for (i=0; i<size; i++) {
	 pixel = vread(i,j); 
         fputc (pixel, pltout);			/* send binary byte */
/*     }
    }
*/


}

/**************************************************/

pointf(void)
{
  vwrite(x,y,scrcolor);
}

/**************************************************/

clrcf (void)

/* clear screen and home up */

{
}

/**************************************************/
