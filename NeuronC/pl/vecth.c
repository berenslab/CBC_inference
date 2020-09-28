/* segment VECT in program SPIKE */

/* latest mod   Oct 1985     R.G.Smith */

/* subroutines to draw vectors on IBM Enhanced graphics card */

#include "../h/colors.h"

#define VIDMAX 349             /* resolution of screen */
#define VIDSIZ 350             /* size of graphics screen */
#define YSCALEDIV  47          /* divide or matrox from STDSIZ */
#define XSCALEDIV  47 
#define STDSIZ 16384            /* size of screen for move() and draw() */

#define abs(a)  ((a) < 0 ? -(a) : (a))


static int      xpos = 0;
static int      ypos = 0;       /* inverted for matrox card */
static int	xloc = 0;
static int	yloc = 0;

static int      oldx=0;        	/* device indep coords, used in "charvec.c" */
static int      oldy=0;
static int	scrncol=7;	/* initial screen color is white */

static int xscalediv = XSCALEDIV; /* used by fscale below */
static int yscalediv = YSCALEDIV; /* used by fscale below */
static int xoffset= 0;
static int yoffset= 0;

/*------------------------------------------------------------*/

drawto (deltax,deltay,color)
   int deltax,deltay,color;
   
/* Bresenham's Algorithm;
   see Newman and Sproul, 2nd Ed.  for explanation.

   color = 1 -> draw line;
   color = 0 -> erase line;

   device = 0 -> draw on matrox;
   device = 0 -> draw on plotter;

 */

{
        static int dx,dy,incr1,incr2,d,x,y,xend,yend;
        static int xinc,yinc,tempx,tempy;

 deltay = -deltay;
 tempx = xpos; tempy = ypos;
 dx = abs(deltax); dy = abs(deltay);

 if (dy <= dx)
   {
     if (deltax < 0)
       {
         yinc = ((deltay < 0) ? 1 : -1); 
         x = xpos + deltax; y = ypos + deltay; xend = xpos;
       }
     else   /* deltax >= 0 */
       {
         yinc = ((deltay < 0) ? -1 : 1);
         x = xpos; y = ypos; xend = xpos + deltax;
       }
     pointx (x,y,color);
     d = (dy << 1) - dx;
     incr1 = dy << 1;
     incr2 = (dy - dx) << 1;
     while (x < xend)
       {
         x++;
         if (d < 0) d += incr1;
         else     { y += yinc; d += incr2; }
         pointx (x,y,color);
       }
    } 
  else   /* dy > dx */
    {
      if (deltay < 0)
        {
          xinc = ((deltax < 0) ? 1 : -1);
          x = xpos + deltax; y = ypos + deltay; yend = ypos;
        }
      else
        {
          xinc = ((deltax < 0) ? -1 : 1);
          x = xpos; y = ypos; yend = ypos + deltay;
        }
      pointx (x,y,color);
      d = (dx << 1) - dy;
      incr1 = dx << 1;
      incr2 = (dx - dy) << 1;
      while (y < yend)
       {
         y++;
         if (d < 0) d += incr1;
         else     { x += xinc; d += incr2; }
         pointx (x,y,color);
       }
   }  /* else */

  xpos = deltax + tempx; ypos = deltay + tempy; /* inverted for Matrox */
  xloc = xpos; yloc = VIDMAX - ypos;

} /* drawto */

/*------------------------------------*/

clrscr(screen)

        int screen;
{
}

/*-------------------------------------*/

moveto(x,y)
        int x,y;
{
 xpos = x; ypos = VIDMAX - y;
 xloc = x; yloc = y;
}

/*--------------------------------------*/

pointx (x,y,color)
   int x,y,color;
{
/*  scr_don(x,y,color); 	/* for Enhanced color graphics */
  hdot(x,y);			/* for Hercules graphics */

/*  printf ("x %d y %d color %d\n",x,y,color);
 */ 
}

/*--------------------------------------*/

movee(x,y)
  int x,y;

/* Device independent routine.
   input coordinates range from 0 to 16384.
*/

{
  oldx = x;
  oldy = y;
  moveto(x+xoffset, y+yoffset);
}

/*--------------------------------------*/

drawe(x,y)

/* device independent routine */

{
  oldx = x;
  oldy = y;
  drawto(x-xloc+xoffset,y-yloc+yoffset,scrncol);
}

/*--------------------------------------*/

setgcol(color)
  int color;

{
   scrncol = color;
}
/*--------------------------------------*/


/* extern char sc_mode,sc_page,sc_attr; */

vidstart()

#define HIRES 0x10		/* Enhanced graphics high res mode (350x640) */
#define ALPHA 3

{
/*   scr_set();
 /*  scr_mod(HIRES);		/* for enhanced graphics only */
/*   scr_mod(2);
   scr_colr(BLACK,WHITE,BLUE); 
 /*  scr_clr();
   scr_colr(WHITE,BLACK,BLUE);
*/
  gmode();			/* for Hercules card */

/*  printf ("mode %d page %d attr %d\n",sc_mode,sc_page,sc_attr); */

}

/*--------------------------------------*/

vidret()

{
/*  scr_mod(ALPHA);  */
  tmode();
}

