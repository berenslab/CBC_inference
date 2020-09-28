/*  Segment MPRINTT in program MONT */

/* System-dependent graphics and text I/O;

   For use with Textronix 4014 terminals and emulators
*/

/*  Latest mod  18-Jun-86	R.G.Smith       */

#include <stdio.h>
#include "../h/mdef.h"
#include "../h/colors.h"
#include "../h/grdef.h"

#define XHI  512
#define YHI  512

#define XMAX 79
#define YMAX 24 

#define CWIDTH 8			/* char box size for enh graph */
#define CHEIGHT 14

#define ERASE	014
#define ESC	033
#define VFAST	1		/* */
#define GS	035
#define pltout	stdout

static int xloc,yloc;			/* location of last move */ 
static int cwidth=CWIDTH;
static int cheight=CHEIGHT;

extern int scrntyp;			/* monochrome PC vidio screen adapter */

/*-----------------------------------------*/
  
clrctx (void)

/* clear alpha screen and home up */

{
 putc (ESC,pltout);
 putc (ERASE,pltout);
}

/************************************/

ttxinit(void)
{
}

/*****************************/

ttxexit(void)
{
} 

/*****************************/

settxcol (int color)
{
  int oldcolor, temp;

  temp = oldcolor;
  oldcolor = color; 
  return (temp);
}

/*****************************/

xyvidvtx(int x, int y)
{
  movetx (x,y);
}

/*****************************/


xyvidtx(int x, int y)
{
  x = x * cwidth;
  y = y * cheight;
  movetx (x,y);
}

/*****************************/

static char xhi,yhi,ylo;

movetx(int x, int y)
           

/* Device independent routine which scales from
   STDSIZ (16384) to 512 for TX4014 emulators. 
   This code takes advantage of vector addresses. */

{

  x >>= 5;
  y >>= 5;

  putc(GS,pltout);
  putc((y>>5)|040,pltout); 	/* y/32+32 hi y byte */
  putc((y&037)|0140,pltout); 	/* y%32+96 low y byte */
  putc((x>>5)|040,pltout); 	/* x/32+32 hi x byte */
  putc((x&037)|0100,pltout); 	/* x%32+64 low x byte */

#ifdef VFAST
{
  yhi = (y>>5)|040;
  ylo = (y&037)|0140;
  xhi = (x>>5)|040;
}
#endif

} 

/*****************************/

drawtx(int x, int y)
          

#ifndef VFAST
   {
        x >>= 5;
        y >>= 5;

	putc((y>>5)|040,pltout); 	/* y/32+32 hi y byte */
	putc((y&037)|0140,pltout); 	/* y%32+96 low y byte */
	putc((x>>5)|040,pltout); 	/* x/32+32 hi x byte */
	putc((x&037)|0100,pltout); 	/* x%32+64 low x byte */
   }
#else
   {
     register char tempx, tempy;

        x >>= 5;
        y >>= 5;

	tempy= (y>>5)|040;	/* y/32+32 hi y byte */
	if(tempy != yhi)
	   {
		putc(tempy,pltout);
		yhi= tempy;
	   }
	tempy= (y&037)|0140;	/* y%32+96 low y byte */
	tempx= (x>>5)|040;	/* x/32+32 hi x byte */
	if(tempx != xhi)
	   {
		putc(tempy,pltout);
		putc(tempx,pltout);
		ylo= tempy;
		xhi= tempx;
	   }
	  else if (tempy != ylo)
	   {
		putc(tempy,pltout);
		ylo= tempy;
	   }
	putc( (x&037)|0100,pltout);	/* x%32+64 low x byte */
   }
#endif

/***********************************************/

int ttxsize(void)
{
  return(STDSIZ);
}

/***********************************************/

int ttysize(void)
{
  return(STDSIZ);
}

