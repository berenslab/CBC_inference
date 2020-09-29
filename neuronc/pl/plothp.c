/* Segment PLOTHP */

/* Subroutines for HP 7221A plotter,
to be linked with an outside main(). */

/* Based on Neil Friedman's Fortran
drivers written in Nov, 1978 */

/*  Latest Mod 	26-Mar-82		R.G.Smith	*/

#include <stdio.h>
#include "../h/mdef.h"


/****************************************/

extern psend ();

/* This routine must be in the main segment.
This way, we can direct output wherever we want. */

/* Routine to put char. in output buffer
for plotter.  When buffer is full, it is
dumped to the standard output.  Use this
routine for all graphics commands to plotter.
*/

/* begin
  putc (ch,stdout);
end */

/****************************************/

psend2 (ch1,ch2)
    char ch1,ch2;

/* Send two chars to buffer. */

begin
  psend (ch1);
  psend (ch2);
end

/****************************************/

preset ()

/* Find mechanical limit switch on plotter;
reset all plotter graphics variables. */

begin
 psend2 ('~','_');
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
 plimit(3100,300,13275,10475);
 pgrid (1023,1023);
 plsize (15,30);
 pfont (0,5);
 ppen (1);
 papnt (512.,512.);
end



/****************************************/

pset2()

begin
 pinit();
 plimit (520,380,15720,10380);
 pgrid  (15200,10000);


end

/****************************************/

pveloc (speed)

/* Set max velocity of the plotter pen. */

begin
 if (speed < 1) speed = 1;
 if (speed > 36) speed = 36;
 psend2 ('~','V');
 psbn (speed);
end

/****************************************/

papnt (x,y)
    int x,y;

/* Position the plotter pen at abs postion x,y.  */

begin
  psend ('p');
  pmbp (x,y);
end

/****************************************/

pdraw (x,y)
   int x,y;

/* Move to abs loc (x,y) with the pen down. */

begin
  psend ('q');
  pmbp (x,y);
end

/****************************************/

ppflush()

/* command to tell hplot to pflush(). */

begin
  psend ('z');
end

/****************************************/

pvect (x,y)
  int x,y;

/* Draw a line from current pos. to new
pos. using x,y as delta values. */

begin
 psend ('s');
 ppmb (x,y);
end

/****************************************/

pcrlf()

begin
  ptext (1,"\r\n");
end

/****************************************/

pdash (type,vect,len,units)
      int type,vect[],len,units;

/* Routine to set up dashed lines on the plotter.

	type 	Determines whether a solid line,
		a fixed dashed or variable dashed
		line is desired (see pp. 152-59 of HP man.)

	vect	Array that contains the in order dash and
		space lengths that comprise the dashed line.

	len	Length of vect

	units	Number of plotter units that the pattern fits into.

*/

begin
  int i,j;
  char b;
  BOOL sflag;

 psend ('~');
 if (type == 2) psend ('R');		/* Send dash command header */
 else           psend ('Q');
 if (type < 1 or type > 2) return;
 
 sflag = F;
 for (i=0; i<len; i++)			/* Write vect in msbn format */
  begin
   j = vect [i];
   if (j < 0) j = 0;
   if (j > 31) j = 31;
   if (sflag) b = j + 64;
   else        b = j + 32;
   psbn (b);
   sflag = not sflag;
  end
 
 j = units;				/* Send length to plotter. */
 if (j < 0) j = 1;
 pmbn (j);

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
  psend ('t');
  pmbn (radius);
end


/****************************************/

ppen (pen)
  int pen;
/* Select a pen on the plotter.  Pen ranges
from 1-4; if any other value is sent, the 
current pen is returned to its stall and no
new pen is loaded. */

begin
 psend ('v');
 if      (pen < 0) pen = 0;
 else if (pen > 4) pen = (pen-1) % 4 + 1;
 psbn (pen);
end

/****************************************/

pdown()

/* Routine to put pen down. */

begin
  psend ('q');
end

/****************************************/

pup()

begin
 psend ('p');
end

/****************************************/

pfont (stand,alt)
    int stand,alt;

/* Routine to set the standard and alternate char
fonts on the plotter.

	Range of fonts:  0 - 5

	Alternate font is used when ^N is sent to plotter.
	Standard font resumes with a ^O.
*/

begin
  psend2 ('~','P');
  pmbp (stand,alt);
end

/****************************************/

ptext (type,str)
   int type;
   TEXT *str;

/* Send text string to plotter.  String must
end in zero.  Plotter prints all chars until
it receives a 03. */

begin
 if ((type < 2) or (type > 6)) type = 1;
 type--;
 pfont (type,type);
 psend2 ('~','\047');
 while (*str) psend (*str++);
 psend ('\03');
end

/****************************************/

pnmbr (value,form)
    float value;
    char form;

/* Routine to write a number on the plotter.
FORM is = f for floating, i for integer. */

begin
  TEXT tbuf [10];

 if (form == 'i') sprintf (tbuf,"%d",value);
 if (form == 'f') sprintf (tbuf,"%.2f",value);
 ptext (1,tbuf);
end


/****************************************/

pgrid (x,y)
   int x,y;

/* Routine to set the grid size on the plotter. */

begin
  psend2 ('~','S');
  pmbp (x,y);
end


/****************************************/

plimit (llx,lly,urx,ury)
   int  llx,lly,urx,ury;

/* Set the graphic limits on the plotter. */

begin
 psend2 ('~','W');
 pmbp (llx,lly);
 pmbp (urx,ury);
end

/****************************************/

plsize (x,y)
   int  x,y;

/* Set the label size on the plotter. */

begin
 psend2 ('~','%');
 pmbp  (x,y);
end

/****************************************/

prdot (x,y)
  int x,y;

/* Routine to position the pen at relative coords. */

begin
  psend ('r');
  ppmb (x,y);
end

/****************************************/

protat (angle)
   float angle;

/* Rotate all drawing on plotter by ANGLE. */

begin
 psend ('w');
 pmba (angle);
end

/****************************************/

pmba (angle)
  float angle;

/* Routine for multiple byte angle. */

begin
   FAST char b1,b2,b3;
   FAST int iangle;
   int i;

  b1 = 0x60;				/* Code bits for first char */
  b2 = b3 = 0;
  while (angle > 360.) angle -= 360.;
  while (angle < 0.0)  angle += 360.;
  if (angle >= 180.0)
   begin
    b1 |= 0x08;				/* Or in bit15 of bin angle */
    angle -= 180.0;
   end
  iangle = angle * 182.0444;		/* Make rem into integer (16384/90) */
  b3 = iangle & 0x3f;
  iangle >>= 6;				/* bite off into 6 bit chunks */
  b2 = iangle & 0x3f;
  b1 |= (iangle >> 6) & 0x07;

  psend (b1);
  i = 0;
  if (b2) i=1;
  if (b3) i=2;				/* Send 2 or 3 bytes only if needed */
  if (!(b2 & 0x20)) b2 |= 0x40;		/* Make complement bits */
  if (!(b3 & 0x20)) b3 |= 0x40;
  if (i >= 1) psend (b2);
  if (i == 2)  psend (b3);

end

/****************************************/

pmbn(word)
   int word;

/* Routine for multiple byte number. HP man, p.65. */

begin
  FAST int w1,w2,w3;
  int i;

 i=0;
 word = abs (word);
 if (word >= 1024)
  begin
   w3 = word & 0x3F;
   if (!(w3 & 0x20)) w3 |= 0x40;
   word >>= 6;
   i++;
  end
 if (word >= 16)
  begin
   w2 = word & 0x3f;
   if (!(w2 & 0x20)) w3 |= 0x40;
   word >>= 6;
   i++;
  end
 w1 = word | 0x60;
 psend (w1);
 if (i >= 1) psend (w2);
 if (i == 2) psend (w3);
 
end

/****************************************/

pmbp(word1,word2)
   int word1,word2;

begin
  int i,wmax;
  FAST int w1,w2,w3,w4,w5;
 
  word1 = abs (word1);
  word2 = abs (word2);
  if (word1 > 16383) word1 = 16383;
  if (word2 > 16383) word2 = 16383;

  i = w1 = w2 = w3 = w4 = w5 = 0;
  wmax = max (word1,word2);
  if (wmax >= 2048)
   begin
    w1 = (word1 >> 10) | 0x60;
    w2 = (word1 & 0x3ff) >> 4;
    w3 = ((word1 & 0x0f) << 2) + (word2 >> 12);
    w4 = (word2 & 0xfc0) >> 6;
    w5 = word2 & 0x3f;
    i = 5;
   end
  else if (wmax >= 256)
   begin
    w1 = (word1 >> 7) | 0x60;
    w2 = (word1 & 0x7e) >> 1;
    w3 = ((word1 & 1) << 5) + (word2 >> 6);
    w4 = word2 & 0x3f;
    i = 4;
   end
  else if (wmax >= 32)
   begin
    w1 = (word1 >> 4) | 0x60;
    w2 = ((word1 & 0x0f) << 2) + (word2 >> 6);
    w3 = word2 & 0x3f;
    i = 3;
   end
  else if (wmax >= 4)
   begin
    w1 = (word1 >> 1) | 0x60;
    w2 = ((word1 & 1) << 5) + word2;
    i = 2;
   end
  else 
   begin
    w1 = (word1 << 2) | word2 | 0x60;
    i = 1;
   end
  if (!(w2 & 0x20)) w2 |= 0x40;
  if (!(w3 & 0x20)) w3 |= 0x40;
  if (!(w4 & 0x20)) w4 |= 0x40;
  if (!(w5 & 0x20)) w5 |= 0x40;

  psend (w1);
  if (i > 1) psend (w2);
  if (i > 2) psend (w3);
  if (i > 3) psend (w4);
  if (i > 4) psend (w5);
 
end

/****************************************/

ppmb (word1,word2)
  int word1,word2;

begin
  int i,j;
  int word,mask;
  int w1,w2,w3;
 
 for (j=0; j<2; j++)		/* Run this twice:  for x, then y */
  begin
   if (j == 0)
    begin
     word = word1;
     mask = 0x40;
    end
   else
    begin
     word = word2;
     mask = 0x20;
    end

   i = 0;
   if ((word > 511) or (word < -512))
    begin
     w3 = (word & 0x1f) | mask;
     i++;
     word >>= 5;
    end
   if ((word > 15) or (word < -16))
    begin
     w2 = (word & 0x1f) | mask;
     i++;
     word >>= 5;
    end
   w1 = (word & 0x1f) | mask;
   psend (w1);
   if (i > 0) psend (w2);
   if (i > 1) psend (w3);
  end
end

/****************************************/

psbn (byt)
   int byt;

begin
 if (!(byt & 0x20)) byt = (byt & 0x3f) | 0x40;
 psend (byt);
end

/****************************************/
