/*  Segment MPRINTM in program MONT */

/* System-dependent graphics and text I/O;

   For use with the VENIX-supplied drivers for:

		Monochrome Graphics
		Color Graphics
		Enhanced Graphics
*/

/*  Latest mod  5-Mar-88	R.G.Smith       */

#include <stdio.h>
#include "../h/mdef.h"
#include "../h/colors.h"
#include "../h/grdef.h"

#define XMHI STDSIZ
#define YMHI STDSIZ

#define XMAX 79				/* Logical coords for alphanumerics */
#define YMAX 24 

#define CWIDTH  (6  * 48)		/* char box size for enh graph */
#define CHEIGHT (14 * 48)

static int xloc,yloc;			/* location of last move */ 
static int cwidth=CWIDTH;
static int cheight=CHEIGHT;

extern int scrntyp;			/* MONOCHR HERCULES ENHANCED */
extern int scrcolor;			/* color of screen for color boards */

/********************************************/

putcm(ch)
   char ch;

/* Display on primary graphics screen, overlaid on graphics. */

{
   static char str[2]=0;

  str[0]=ch;
  move (xloc,yloc);
  switch (ch)
   begin
      case 0x08:
	  xloc -= CWIDTH;
	  break;
      case '\r':
	  xloc = 0;
	  break;
      case '\n':
	  xloc = 0;
	  yloc -= CHEIGHT;
	  break;
      default:
          label(str);
          xloc += CWIDTH;
	  break;
       end      /* switch */

 move(xloc,yloc);
}

/********************************************/

putsm(str)
   TEXT *str;

/* Display on primary graphics screen, overlaid on graphics. */

{
 for (; *str; str++)
  putcm(*str);
}

/*-----------------------------------------*/
  
clrcm ()

/* clear alpha screen and home up */

begin
  erase();
  xloc = 0;
  yloc = YMHI;
  setorg (xloc,yloc);
end

/************************************/

tminit()

begin
  openpl();
  space(0,0,XMHI,YMHI);
end

/*****************************/

tmexit()

begin
  closepl();
  system ("erspc");
end 

/*****************************/

int setmcol (scolor)
   int scolor;

begin
   static int oldcolor=0, tcolor;

  tcolor = oldcolor;
/*  scr_colr(WHITE,BLACK,0); */
  color(scolor);
  oldcolor = scolor;
  return (tcolor);
end

/*****************************/

xyvidm (xpos,ypos)
  int xpos,ypos;

/* Set text cursor (close) to a graphics position.
   The VENIX plot(3x) drivers are used.  
 */

begin
  move (xpos,ypos);
  xloc = xpos;
  yloc = ypos;
end

/*****************************/

xyvidmg (xpos,ypos)
    int xpos,ypos;

/* set text cursor to a graphics position in GRAPHICS mode */

begin
  move (xpos,ypos);		/* move graphics loc */
  setorg (xpos,ypos);  		/* set the text loc */
end

/*****************************/

xyconm (x,y)
   int x,y;

/* X,Y cursor addressing for alphanumerics on console */
/* x = columns, y = rows; (0,0) is at bottom left */

begin

 x = min(x,XMAX);
 y = min(y,YMAX);
 x = max(x,0);
 y = max(y,0);
 xloc = x * cwidth;
 yloc = y * cheight;
 move(xloc,yloc);
end

/********************************************/

xycong (x,y)
   int x,y;

begin
 x = min(x,XMAX);
 y = min(y,YMAX);
 x = max(x,0);
 y = max(y,0);
 xloc = x * cwidth;
 yloc = y * cheight;
 xyvidmg(x*cwidth,y*cheight);
end

/********************************************/

int tmxsize()

{
  return (STDSIZ);
}

/********************************************/

int tmysize()

{
  return (STDSIZ);
}

/********************************************/



