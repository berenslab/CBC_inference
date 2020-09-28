/*  Segment MPRINTL in program MONT */

/* System-dependent graphics and text I/O;

   Designed for use with HP LaserJet II printer.
*/

/*  Latest mod  5-Mar-88	R.G.Smith       */

#include <stdio.h>
#include <fcntl.h>
#include "mdef.h"
#include "colors.h"

#include <signal.h>
void (*istat)();		/* Used by interrupt routine setup below */
				/* see "UNIX Programming, second edition" */

#define STDSIZ 16384

#define ESC 033
#define FF  014

#define RASRES  280		/* resolution for laserjet (lines/inch) */
#define XRAS (8*RASRES)		/* x raster size for array */
#define YRAS (10*RASRES)	/* y raster size for array */

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

extern int scrcolor;			/* color of screen for enhanced adapt */
extern int scrrev;			/* =1 -> screen reversed vertically */

/********************************************/

onintr(void)
{
 virtend();
 exit();
}

/**************************************************/

putsl(TEXT *str)
             

/* Display string on graphics. */

{
 for (; *str; str++)				/* enhanced color monitor */
  putcl(*str);
}

/*-----------------------------------------*/

putcl(char ch)
           

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

 xyconl(xloc/CWIDTH,yloc/CHEIGHT);
}

/************************************/
  
tlinit(void)
{

if (! gmodefl)			/* if this routine hasn't already been called */
  virtinit();			/* initialize the virtual raster array */

istat = signal (SIGINT, SIG_IGN);	/* save orig staus */
if (istat != SIG_IGN)
  (void) signal (SIGINT,(__sighandler_t)onintr);   /* go to "onintr" at ^C */
 grflag = 0;
 gmodefl = 1;
/* scrcolor = WHITE; */
}

/*--------------------------------------*/

tlexit(void)
{
  virtend();
}

/*--------------------------------------*/


grdump(void)

/* output graphics to printer after tlexit() */

{

if (grflag) {				/* if graphics were used */

fprintf (pltout,"%%!PS-Adobe-1.0\n");
fprintf (pltout,"%%%%BoundingBox: 0 0 %d %d\n",XRAS,YRAS);
fprintf (pltout,"%%%%Creator: sim:rob \n");
fprintf (pltout,"%%%%Title: vid picture \n");
fprintf (pltout,"%%%%CreationDate: Tue Feb 25 18:12:34 1992\n");
fprintf (pltout,"%%%%EndComments\n");
fprintf (pltout,"%%%%Pages: 1\n");
fprintf (pltout,"%%%%EndProlog\n");
fprintf (pltout,"%%%%Page: 1 1\n");
fprintf (pltout,"\n");
fprintf (pltout,"/bitdump %% stk: width, height, iscale\n");
fprintf (pltout,"%% dump a bit image with ll corner at current origin,\n");
fprintf (pltout,"%% scaling by iscale (iscale=1 means 1/300 inch per pixel)\n");
fprintf (pltout,"{\n");
fprintf (pltout,"	%% read arguments\n");
fprintf (pltout,"	/iscale exch def\n");
fprintf (pltout,"	/height exch def\n");
fprintf (pltout,"	/width exch def\n");
fprintf (pltout,"\n");
fprintf (pltout,"	%% scale appropriately\n");
fprintf (pltout,"	width iscale mul height iscale mul scale\n");
fprintf (pltout,"\n");
fprintf (pltout,"	%% allocate space for one scanline of input\n");
fprintf (pltout,"	/picstr %% picstr holds one scan line\n");
fprintf (pltout,"	   %% width of image in bytes = ceiling(width/8)\n");
fprintf (pltout,"	  width 7 add 8 idiv\n");
fprintf (pltout,"		string\n");
fprintf (pltout,"		def\n");
fprintf (pltout,"\n");
fprintf (pltout,"	%% read and dump the image\n");
fprintf (pltout,"	width height 1 [width 0 0 height neg 0 height]\n");
fprintf (pltout,"	{ currentfile picstr readhexstring pop }\n");
fprintf (pltout,"	image\n");
fprintf (pltout,"} def\n");
fprintf (pltout,"1 1 scale\n");
fprintf (pltout,"20 20 translate\n");
fprintf (pltout,"%d %d %g bitdump\n",XRAS,YRAS,0.27);

  rastout();				/* send raster data to printer */

fprintf (pltout,"\n");
fprintf (pltout,"showpage\n");
fprintf (pltout,"%%%%Trailer\n");

  grflag=0;
}

if (gmodefl) virtend();	 		/* erase virtual raster file */

gmodefl = 0;				/* reset graphics mode flag */

}

/*****************************/

int setlcol (int color)
{
   static int oldcolor;

  oldcolor = scrcolor;
  scrcolor = color;
  return (oldcolor);
}

/*****************************/

xyvidl (int xpos, int ypos)
                  

/* set text cursor to a graphics position in GRAPHICS mode */

{
  movel (xpos,ypos);		/* move graphics loc */
  setorg (xpos,ypos);  		/* set the text loc */
}

/*****************************/

xyconl (int x, int y)
           

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
 xyvidl(xloc,yloc);
}

/********************************************/

int tlxsize(void)
{
/*  return ((STDSIZ / yhi) * xhi); */
  return (STDSIZ);
}

/********************************************/

int tlysize(void)
{
/*  return (STDSIZ); */
  return (STDSIZ * 1.25);

}

/**************************************************/

ldot (int x, int y, int color)

/* Draw one dot only. Not currently used. */

{
return;
/*setlcol(color);
movel (x,y);
drawl (x,y);
*/
}

/*--------------------------------------*/

double xscall(double size)
{
   extern int xscale,xhi;
   int temp;

  temp = xscale;
/*  xscale = ((double)XSCALE) / size; */
/*  xhi = (STDSIZ / xscale) - 2; */
  return(((double)temp);
}

/*--------------------------------------*/

double yscall(double size)
{
   extern int yscale,yhi;
   int temp;

  temp = yscale;
  yscale = ((double)YSCALE) / size;
/*  yhi = (STDSIZ / yscale) - 2; */
  return(((double)temp);
}

/*--------------------------------------*/

movel(int x, int y)
          

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
  if (! scrrev) ypos = yhi-ypos;
}

/*--------------------------------------*/

drawl(int x, int y)
          

/* draw to (x,y) routine for Laserjet II.
   Square screen size is XRAS * XRAS (1200).
   Scale x by 14 (16384/1200).
   Scale y by 14. 
*/

{
   static int lindat[5] = {0};
   static int *linpnt = lindat;
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

  linsegl(xpos,ypos,x,y);
  xpos = x;
  ypos = y;
  grflag = 1;			/* set graphics flag when graphics used */
}

/**************************************************/

linsegl (int x1, int y1, int x2, int y2)
                   
   
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
     pointl();
     d = (dy << 1) - dx;
     incr1 = dy << 1;
     incr2 = (dy - dx) << 1;
     while (x < xend)
       {
         x++;
         if (d < 0) d += incr1;
         else     { y += yinc; d += incr2; }
         pointl();
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
      pointl();
      d = (dx << 1) - dy;
      incr1 = dx << 1;
      incr2 = (dx - dy) << 1;
      while (y < yend)
       {
         y++;
         if (d < 0) d += incr1;
         else     { x += xinc; d += incr2; }
         pointl();
       }
   }  /* else */

} /* drawto */

/**************************************************/

gmodel(void)

/* set graphics mode */

{
}

/**************************************************/

tmodel(void)

/* return to text mode */

{
}

/**************************************************/

rastout(void)
{
   int i,j,k,xras;
   short int grbyte;
   static short int bittab[8] = {0x80,0x40,0x20,0x10,8,4,2,1};

  xras = XRAS / 8;

  for (j=0; j<YRAS; j++) {			/* for all rows */
    /* fprintf (pltout,"%dW",xras);		/* row header */
    for (i=0; i<XRAS; i+=8) {
      grbyte = 0;
      for (k=0; k<8; k++) {
	if (vread(j,i+k)) 
           grbyte |= bittab[k];			/* get bit value */
      }
      fprintf (pltout,"%2.2x",~grbyte & 0xff);	/* translate to hex */
    }
    fprintf (pltout,"\n");			/* end of row */
  }
}

/**************************************************/

pointl(void)
{
  vwrite(y,x);
}

/**************************************************/

clrcl (void)

/* clear screen and home up */

{
/*  putc (FF,pltout); */
}

/**************************************************/
