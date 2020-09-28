/* Program PLABEL */

/* Program to write labels on pictures. */

/* Latest mod   3-Mar-82			R.G.Smith	*/


/* 	Input statement formats:

	t x y "text"		print text at loc x,y
	l x1 y1 x2 y2		draw line from x1,y1 to x2,y2
	s x y			make following text be size x y
	p x			get pen color x (0 = put back)
	v x			make velocity for draw = x
	g x y			make grid = x,y
	m x y			make plotting limits be x,y
*/ 	
	

#include <stdio.h>
#include "../h/mdef.h"

main()

begin
   int x,y,x2,y2;
   TEXT label [80], str [80];
   char ch;

  pset2();				/* Find mechanical limit */
  					/* set limits, grid prop to limits */

  while ((ch = getc (stdin)) > 0)
   begin
    x = y = x2 = y2 = 0;
    *label = NULL; *str = NULL;
    switch (ch)
     begin

    case 'T':					/* Text at pos. x,y */
    case 't': scanf ("%d %d %*c%[^$/]", &x,&y,label);
    	      papnt (x,y);
	      ptext (1,label);
	     break;

    case 'C':					/* Continue text */
    case 'c': scanf ("%*c%[^$/]",label);
	     break;

    case 'L':					/* Line */
    case 'l': scanf ("%d %d %d %d", &x,&y, &x2, &y2);
 	      papnt (x,y);
	      pvect (x2-x,y2-y);
	     break;

    case 'S':					/* Size of chars */
    case 's': scanf ("%d %d", &x, &y);
	      plsize (x,y);
	     break;

    case 'P':					/* Pen number */
    case 'p': scanf ("%d", &x);
	      ppen (x);
	     break;

    case 'V':					/* draw velocity */
    case 'v': scanf ("%d",&x);
	      pveloc (x);
	     break;

    case 'G':					/* Grid size */
    case 'g': scanf ("%d %d", &x,&y);
	      pgrid (x,y);
	     break;

    case 'M':					/* plot limits */
    case 'm': scanf ("%d %d", &x,&y);
	      plimit (x,y);
	     break;

    default: break;

     end

    gets (str);

/*    fprintf (stderr,"%c %d %d '%s'\r\n",ch,x,y,label);
*/
   end

end

/****************************************/

psend (ch)
   char ch;

/* routine to send char to the plotter */

begin
 putc (ch,stdout);
end

