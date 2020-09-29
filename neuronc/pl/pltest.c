
#include <stdplt.h>

double cwidth();

main()
{
   static double x,y;

	move (0.5,0.5);
        draw (0.9,0.8);
	x = 0.5;
	y = cwidth (x,".");
fprintf (stderr,"x %6.3f y %6.3f\n",x,y);
}
