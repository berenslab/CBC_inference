/* Segment PLOTHP */

/* Subroutines for IBM 7372 color plotter,
to be linked with an outside main(). */

/* Based on Neil Friedman's Fortran
drivers written in Nov, 1978 */

/*  Latest Mod 	26-Sept-86		R.G.Smith	*/

#include <stdio.h>
#include "../h/mdef.h"


/****************************************/

extern psend ();

/* This routine must be in the main segment.
This way, we can direct output wherever we want. */

/* Routine to put char string in output buffer
for plotter.  When buffer is full, it is
dumped to the standard output.  Use this
routine for all graphics commands to plotter.
*/

/* begin
  fprintf (stdout,str);
end */

/****************************************/

/****************************************/

preset ()

/* Find mechanical limit switch on plotter;
reset all plotter graphics variables. */

begin
 psend ("DF;IN;");
end


/****************************************/

pinit()

/* Routine to initialize the plotter.
Initialize the mechanical limits on plotter,
and return current pen to its stall. */

begin
 preset();		/* Reset mechanical limits */
 ppen (0);		/* Return current pen to its stall */
end

/****************************************/

psetup()

begin
 pinit();
end

/****************************************/

pset2()

begin
 pinit();
end

/****************************************/

pveloc (speed)
   int speed;

/* Set max velocity of the plotter pen. */

begin
 if (speed < 1) speed = 1;
 if (speed > 36) speed = 36;
 psend ("VS %d;",speed);
end

/****************************************/

papnt (x,y)
    int x,y;

/* Position the plotter pen at abs postion x,y.  */

begin
  psend ("PU");
  psend ("PA %d %d;",x,y);
end

/****************************************/

pdraw (x,y)
   int x,y;

/* Move to abs loc (x,y) with the pen down. */

begin
  psend ("PD");
  psend ("PA %d %d;",x,y);
end

/****************************************/

pvect (x,y)
  int x,y;

/* Draw a line from current pos. to new
pos. using x,y as delta values. */

begin
  psend ("PD");
  psend ("PR %d %d;",x,y);
end

/****************************************/

pcrlf()

begin
  ptext ("\r\n");
end

/****************************************/

pcircl (radius)
    int radius;

/* Draw a circle on the plotter.
The center of the circle is at absolute
position (x-radius,y) where x and y are
the current values for the pen's position.
Circle is drawn clockwise. */

begin
 psend("CI %d;",radius);
end


/****************************************/

ppen (pen)
  int pen;
/* Select a pen on the plotter.  Pen ranges
from 1-4; if any other value is sent, the 
current pen is returned to its stall and no
new pen is loaded. */

begin
 if      (pen < 0) pen = 0;
 else if (pen > 6) pen = (pen-1) % 6 + 1;
 psend ("PU");
 psend ("SP %d;",pen);
end

/****************************************/

pdown()

/* Routine to put pen down. */

begin
  psend ("PD;");
end

/****************************************/

pup()

begin
 psend ("PU;");
end

/****************************************/

pfont (charset)
    int charset;

/* Routine to set the standard font on the plotter.

	Range of fonts:  0 - 39

	Allowable fonts: 0-4, 6-9, 30-39

*/

begin
  psend ("CS %d;",charset);
end

/****************************************/

ptext (str)
   TEXT *str;

/* Send text string to plotter.  String must
end in zero.  Plotter prints all chars until
it receives a 03. */

begin
 psend (str);
 psend ("\03");
end

/****************************************/

plimit (xmin,xmax,ymin,ymax)
    int xmin,xmax,ymin,ymax;

/* Set the graphic limits on the plotter. */

begin
end

/****************************************/

plsize (x,y)
   int  x,y;

/* Set the label size on the plotter. */

begin
 psend ("SI %d %d;",x,y);
end

/****************************************/

prdot (x,y)
  int x,y;

/* Routine to position the pen at relative coords. */

begin
end

/****************************************/

protat (angle)
   float angle;

/* Rotate all drawing on plotter by ANGLE. */

begin
end

/****************************************/
