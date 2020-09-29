/* module gprim.c in program nc */

/* graphics primitives for video and text modes */

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "stdplt.h"
#include "gr.h"

#ifdef __cplusplus
}
#endif

#include "ncio.h"
#include "extdraw.h"

extern struct _frame *_dotp;
extern int vidmode;

void gmove (double x, double y)
{
   if (vidmode) move (x,y); 			/* in libP.a */
   else ext_gmove(x, y);
}

void gdraw (double x, double y)
{
   if (vidmode) draw (x,y);			/* in libP.a */
   else ext_gdraw(x, y);
}

void grmove (double x, double y)
{
   if (vidmode) rmove (x,y); 			/* in libP.a */
   else ext_grmove(x, y);
}

void grdraw (double x, double y)
{
   if (vidmode) rdraw (x,y); 			/* in libP.a */
   else ext_grdraw(x, y);
}

void gpen (int x)
{
   if (vidmode) cpen (x); 			/* in libP.a */
   else ext_gpen(x);
}

void gtext(const char *s)
{
  if (vidmode) textf ("%s",s);			/* in libP.a */
  else ext_gtext((char *)s);
}

void glabel(const char *s, double x, double y)
{
  gmove (x, y);
  if (vidmode) textf ("%s",s);			/* in libP.a */
  else ext_gtext((char *)s);
}

void glabel(const char *s, int color, double x, double y)
{
  gmove (x, y);
  gpen (color);
  if (vidmode) textf ("%s",s);			/* in libP.a */
  else ext_gtext((char *)s);
}

void gcrotate (double x)
{
   if (vidmode) crotate (x*DEG,(char *)"."); 		/* in libP.a */
   else ext_gcrotate(x);
}

void grotate (double x)
{
 if (vidmode) rotate (x*DEG);
 else ext_grotate(x);
}

void gorigin(double x, double y)
{
 if (vidmode) origin (x,y);
 else ext_gorigin(x, y);
}

void gframe (const char *s)
{
 frame ((char *)s);
 if (!vidmode) ext_gframe((char *)s);
}

void grmframe (const char *s)
{
 if (vidmode) rmframe ((char *)s);
 else ext_grmframe((char *)s);
}

void gcwidth (double x)
{
 if (vidmode) cwidth (x,(char *)".");
 else ext_gcwidth(x);
}

void gsize (double x)
{
 if (vidmode) scale (x,x);
 else ext_gsize(x);
}

void gdash (int x)
{
 if (vidmode) dash (x);
 else ext_gdash(x);
}

void gwindow (double x1, double y1, double x2, double y2)
{
 if (vidmode) window (x1,y1,x2,y2);
 else ext_gwindow(x1, y1, x2, y2);
}

void gcirc (double rad, int fill)
{
 if (vidmode) _drcirc (rad,fill);
 else ext_gcirc(rad, fill);
}

void grect (double x1, double y1, double x2, double y2,
	    double x3, double y3, double x4, double y4, int fill)
{
 if (vidmode) _drrect (x1,y1,x2,y2,x3,y3,x4,y4,fill);
 else ext_grect(x1, y1, x2, y2, x3, y3, x4, y4, fill);
}

void gtri (double x1, double y1, double x2, double y2,
	    double x3, double y3, int fill)
{
 if (vidmode) _drtri (x1,y1,x2,y2,x3,y3,fill);
 else ext_gtri(x1, y1, x2, y2, x3, y3, fill);
}

void gpurge (void)
{
 if (vidmode) purge ();
 else ext_gpurge();
}

void gpage (void)
{
 if (vidmode) cpage ();
 else ext_gpage();
}

void ghinit (void)
{
 if (vidmode) hidinit ();
 else ext_ghinit();
}

void ghstart (void)
/* start computing hidden line table (in "hide.c") */
{
 if (vidmode) hidstart ();
 else ext_ghstart();
}

void ghend (void)
/* stop computing hidden line table (in "hide.c") */
{
 if (vidmode) hidstop ();
 else ext_ghend();
}

void gctext (void)
/* change to text mode (for mprint{h,e,v}.c drivers) */
{
 if (vidmode) ctext ();
 else ext_gctext();
}

void gcgraphics (void)
/* change to graphics mode (for mprint{h,e,v}.c drivers) */
{
 if (vidmode) cgraphics();
 else ext_gcgraphics();
}

