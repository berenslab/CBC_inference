/* Segment VECT in program SPIKE */

/* latest mod   Oct 1985     R.G.Smith */

/* subroutines to draw vectors on IBM Enhanced graphics card */

#include "../h/colors.h"

#define abs(a)  ((a) < 0 ? -(a) : (a))

extern int scrcolor;		/* color of text and vectors on screen */

static int      xpos = 0;
static int      ypos = 0;
static int	xloc = 0;
static int	yloc = 0;

static int xoffset= 0;
static int yoffset= 0;

/*------------------------------------------------------------*/

drawto (deltax,deltay)
   int deltax,deltay;
   
/* Bresenham's Algorithm;
   see Newman and Sproul, 2nd Ed.  for explanation.

   device = 0 -> draw on matrox;
   device = 0 -> draw on plotter;

 */

{
        static int dx,dy,incr1,incr2,d,x,y,xend,yend;
        static int xinc,yinc,tempx,tempy;

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
     pointx (x,y);
     d = (dy << 1) - dx;
     incr1 = dy << 1;
     incr2 = (dy - dx) << 1;
     while (x < xend)
       {
         x++;
         if (d < 0) d += incr1;
         else     { y += yinc; d += incr2; }
         pointx (x,y);
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
      pointx (x,y);
      d = (dx << 1) - dy;
      incr1 = dx << 1;
      incr2 = (dx - dy) << 1;
      while (y < yend)
       {
         y++;
         if (d < 0) d += incr1;
         else     { x += xinc; d += incr2; }
         pointx (x,y);
       }
   }  /* else */

  xpos = deltax + tempx; ypos = deltay + tempy; /* inverted for Matrox */
  xloc = xpos; yloc = ypos;

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
 xpos = x; ypos = y;
 xloc = x; yloc = y;
}

/*--------------------------------------*/

pointx (x,y)
   int x,y;
{
  tdot (x,y,scrcolor);
/*  scr_don(x,y,scrcolor); 	/* for Enhanced color graphics */
/*  hdot(x,y,scrcolor);		/* for Hercules graphics */

}

/*--------------------------------------*/

moveh(x,y)
  int x,y;

/* Device independent routine for Hercules board.
   Screen size is 348 * 480.
   Physical resolution is 720 in X. This must be multiplied by 3/2 to make 
   square pixels; this gives an overall X resolution of 480.
   Scale by 31 ((2/3)16384/348) in X, 47 (16384/348) in Y.
*/

{
  x /= 31;
  y /= 47;
  moveto(x+xoffset, y+yoffset);
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
  x /= 31;
  y /= 47;
  drawto(x-xloc+xoffset,y-yloc+yoffset,scrcolor);
}

/*--------------------------------------*/

setgcol(color)
  int color;

{
   static int tcolor;

   tcolor = scrcolor;
   scrcolor = color;
   return (tcolor);
}
/*--------------------------------------*/

