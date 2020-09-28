extern "C" {
#include <stdio.h>
}
#include "ncio.h"

void default_gtext(char *s)
{
	ncfprintf(stdout,"#gt %s\n",s);	
	fflush(stdout);
}

void default_gpage(void)
{
	ncfprintf (stdout,"#gpage\n");
	fflush (stdout);   
}

void default_gmove(double x, double y)
{
	ncfprintf(stdout,"#gm %g %g\n",x,y);	
   	fflush(stdout);	
}

void default_gdraw(double x, double y)
{
	ncfprintf (stdout,"#gd %g %g\n",x,y);
	fflush (stdout);	
}

void default_grmove(double x, double y)
{
	ncfprintf (stdout,"#grm %g %g\n",x,y);
	fflush (stdout);		
}

void default_grdraw(double x, double y)
{
	ncfprintf (stdout,"#grd %g %g\n", x, y);
	fflush (stdout);	
}

void default_gpen(int x)
{
	ncfprintf (stdout,"#gpen %d\n",x);	
	fflush (stdout);
}

void default_gcrotate(double x)
{
	ncfprintf (stdout,"#gcro %g\n",x);
	fflush (stdout);	
}

void default_grotate(double x)
{
	ncfprintf (stdout,"#gro %g\n",x);
	fflush (stdout);
}

void default_gorigin(double x, double y)
{
	ncfprintf (stdout,"#gor %g %g\n",x,y);
	fflush (stdout);
}

void default_gframe(char *s)
{
	ncfprintf (stdout,"#gfr %s\n",s);
	fflush (stdout);
}

void default_grmframe(char *s)
{
	ncfprintf (stdout,"#grfr %s\n",s);
	fflush (stdout);
}

void default_gcwidth(double x)
{
	ncfprintf (stdout,"#gcw %g\n",x);
	fflush (stdout);
}

void default_gsize(double x)
{
	ncfprintf (stdout,"#gsiz %g\n",x);
	fflush (stdout);
}

void default_gdash(int x)
{
	ncfprintf (stdout,"#gda %d\n",x);
	fflush (stdout);
}

void default_gcirc(double rad, int fill)
{
	ncfprintf (stdout,"#gcirc %g %d\n",rad,fill);
	fflush (stdout);
}

void default_grect(double x1, double y1, double x2, double y2, double x3, double y3,
		double x4, double y4, int fill)
{
	ncfprintf (stdout,"#grect %g %g %g %g %g %g %g %g %d\n",
 			x1,y1,x2,y2,x3,y3,x4,y4,fill);
	fflush (stdout);	
}
		
void default_gtri(double x1, double y1, double x2, double y2, double x3, double y3,
		int fill)
{
	ncfprintf (stdout,"#gtri %g %g %g %g %g %g %d\n",
 		       x1,y1,x2,y2,x3,y3,fill);
	fflush (stdout);	
}

void default_gwindow(double x1, double x2, double y1, double y2)
{
	ncfprintf (stdout,"#gwind %g %g %g %g\n",x1,y1,x2,y2);
	fflush (stdout);	
}

void default_gpurge(void)
{
	ncfprintf (stdout,"#gpurg\n");
	fflush (stdout);
}

void default_gctext(void)
{
	ncfprintf (stdout,"#gctext\n");
	fflush (stdout);	
}

void default_gcgraphics(void)
{
	ncfprintf (stdout,"#gcgraphics\n");
	fflush (stdout);
}

void default_ghinit(void)
{
	ncfprintf (stdout,"#ghinit\n");
	fflush (stdout);
}

void default_ghstart(void)
{
	ncfprintf (stdout,"#ghstart\n");
	fflush (stdout);	
}

void default_ghend(void)
{
	ncfprintf (stdout,"#ghend\n");
	fflush (stdout);
}
