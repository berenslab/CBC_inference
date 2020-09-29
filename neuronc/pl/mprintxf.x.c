/* Module mprintxf in program vid */

/* XFIG 3.2 driver.  

This module defines device-independent subroutines for "mvid.c".
Nothing real fancy here.  All this program will do for now is
dump out on standard output.  We can redirect at a later stage.
Let's just see it work for now.

*/

#include	<stdio.h>
#include	<time.h>
#include	<sys/types.h>
/* #include	<sys/timeb.h> */
#include	"colors.h"

/* Define image to be 7.5 x 7.5; User can scale
   image to any size if need */

#define STDSIZ 16383.

#define XSCALE 1.820333    /* STDSIZ/(7.5*1200); LW is 1200 dpi */
#define YSCALE 1.820333 

#define GRKFONT 'g'	/* greek font for 's' command */
#define NCOLORS 256

static int col_tran[16] = {0,1,2,3,4,5,6,7,
		           32,10,13,16,19,22,25,7};
                     
static double xpos, ypos;
static double xscale = XSCALE;
static double yscale = YSCALE;
static int xoffset = 600;
static int yoffset = -1200;

#define XHI 9000
#define YHI 9000
static int xhi = XHI;
static int yhi = YHI;
static int scrcolor;
static int fillcolor;
static int fillval=0;
static int depth = 10;

#define DBUFLEN 20000
typedef struct PNTDATA {
	int x;
	int y;
	}pntdata;

static pntdata dbuf[DBUFLEN];
static int npts=0;
static int ncolor=0;

extern char tfont;
extern unsigned int tsize,tangle;
extern char *comment;
extern int scrrev;
extern int backgr;
extern int rotflag;
extern int pargc;
extern char **pargv;

txfinit(void)
{
    printf("#FIG 3.2\n");
    printf("Landscape\n");
    printf("Center\n");
    printf("Inches\n");
    printf("Letter\n");
    printf("100.0\n");
    printf("Single\n");
    printf("-2\n");
    printf("1200 2\n");
    printf("0 32 #808080\n");
    printf("6 %d %d %d %d\n",0+xoffset,-yoffset,xhi+xoffset,yhi-yoffset);
}

txfexit(void)
{
    xfdumppts();
    printf("-6\n");
}
clrcxf(void)
{}
xyconxf(void)
{}

xyvidxf(void)
{}

xfdot(void)
{}

int setxfcol(int pen)
{
   double grayval;
   static int oldcolor;

  oldcolor = scrcolor;
  if (backgr<0 && pen == WHITE) scrcolor = 0;
  else {
     pen = pen % 16; 
     scrcolor = col_tran[pen];
  }
  return (oldcolor);
}

xfdumppts()

{
 int thickness;
 int i,tfill;
 int tdepth;

 if (npts>1) {

   thickness = 2;
   if (fillval) tfill = 20;
   else         tfill = -1;
   printf("2 1 0 %d %d %d %d 0 %d 0.000 0 0 0 0 0 %d\n",
	thickness,ncolor,fillcolor,tdepth=100+ncolor,tfill,npts);
   printf("\t ");
   for (i=0; i<npts; i++) {
	printf("%d %d ", (int)dbuf[i].x,(int)dbuf[i].y);
   }
   printf("\n");
   npts = 0;
 }
}

movexf(int ix, int iy)
{
 extern int scrrev;

 xfdumppts();
 xpos = ix/xscale + xoffset;
 ypos = iy/yscale + yoffset;
 if  (!scrrev) ypos=yhi-ypos;

 npts = 0;
 dbuf[npts].x = xpos;
 dbuf[npts].y = ypos;
 if (++npts >= DBUFLEN) {
	fprintf (stderr,"movexf: too many points in one polygon\n");
	npts= DBUFLEN; 
 }
}

drawxf(int ix, int iy)
{
 extern int scrrev;
 double x,y;

 x = ix/xscale + xoffset;
 y = iy/yscale + yoffset;
 if  (!scrrev) y=yhi-y;
 xpos=x;
 ypos=y;

 ncolor = scrcolor;
 fillcolor = scrcolor;
 dbuf[npts].x = x;
 dbuf[npts].y = y;
 if (++npts >= DBUFLEN) {
	fprintf (stderr,"drawxf: too many points in one polygon\n");
	npts= DBUFLEN; 
 }
}

txfdrcirc(int ix, int iy, int irad, int fill)
{
 double rad,x,y;
 int thickness,tfill,tdepth;

 xfdumppts();
 x = ix/yscale + xoffset;
 y = iy/yscale + yoffset;
 rad = irad/yscale;

 tdepth = 200;
 if (!fill)
   tfill = -1;
 else {
   tfill = 20;
   fillcolor = scrcolor;
 }
 if  (!scrrev) y=yhi-y;
 thickness = 1;
 printf("1 3 0 %d %d %d %d 0 %d 0.000 1 0.000 %d %d %d %d %d %d %d %d\n",
	thickness,ncolor,fillcolor,tdepth+ncolor,tfill,
	(int)x,(int)y,(int)rad,(int)rad,(int)x,(int)y,(int)x+rad,(int)y+rad);

 xpos = x;
 ypos = y;
}

txfdrrect(int ix1, int iy1, int ix2, int iy2, int ix3, int iy3, int ix4, int iy4, int fill)
{
 extern int scrrev;
 int x1,y1,x2,y2,x3,y3,x4,y4;
 int color,tfill,tdepth;
 int thickness,npoints;

 xfdumppts();
 x1 = ix1/xscale + xoffset;
 y1 = iy1/yscale + yoffset;
 x2 = ix2/xscale + xoffset;
 y2 = iy2/yscale + yoffset;
 x3 = ix3/xscale + xoffset;
 y3 = iy3/yscale + yoffset;
 x4 = ix4/xscale + xoffset;
 y4 = iy4/yscale + yoffset;

 tdepth=300;
 if (!fill)
    tfill = -1;
 else {
    tfill = 20;
    fillcolor = scrcolor;
 }
 if  (!scrrev) {
   y1=yhi-y1;
   y2=yhi-y2;
   y3=yhi-y3;
   y4=yhi-y4;
 }

 npoints = 5;
 thickness = 1;
 printf("2 3 0 %d %d %d %d 0 %d 0.000 0 0 0 0 0 %d\n",
	thickness,scrcolor,fillcolor,tdepth,tfill,npoints);
 printf("\t %d %d %d %d %d %d %d %d %d %d\n", 
	    x1,y1,x2,y2,x3,y3,x4,y4,x1,y1);

 xpos = x1;
 ypos = y1;
}

txfdrtri(int ix1, int iy1, int ix2, int iy2, int ix3, int iy3, int fill)
{
 extern int scrrev;
 int x1,y1,x2,y2,x3,y3;
 int color,tfill,tdepth;
 int thickness,npoints;

 xfdumppts();
 x1 = ix1/xscale + xoffset;
 y1 = iy1/yscale + yoffset;
 x2 = ix2/xscale + xoffset;
 y2 = iy2/yscale + yoffset;
 x3 = ix3/xscale + xoffset;
 y3 = iy3/yscale + yoffset;

 tdepth = 400;
 if (!fill)
   tfill = -1;
 else {
   tfill = 20;
   fillcolor = scrcolor;
 }
 if  (!scrrev) {
   y1=yhi-y1;
   y2=yhi-y2;
   y3=yhi-y3;
 }

 npoints = 4;
 thickness = 1;
 printf("2 3 0 %d %d %d 0 %d 0.000 0 0 0 0 0 %d\n",
	thickness,scrcolor,tdepth,tfill,npoints);
 printf("\t %d %d %d %d %d %d %d %d\n", 
	    x1,y1,x2,y2,x3,y3,x1,y1);

 xpos = x1;
 ypos = y1;
}

txffill(int fill)

{
fillval = fill;
}

double xscalexf (double size)
{ 
return ((double)xscale);
}

double yscalexf (double size)
{ 
return ((double)yscale);
}

txfxsize(void)
{ return STDSIZ; }

txfysize(void)
{ return STDSIZ; }

putstrxf(char *str)
{
   static int oldtsize=0,oldtfont=0;
   int height, length, tdepth;
   int x,y,n;
   char *strcat(char *dest, const char *src);

#define MPI 3.14159265358979323846264

  xfdumppts();
  x = xpos;
  y = ypos;
  height = 15;
  tdepth = 500;
  n = strlen(str);
  str = strcat(str,"\\001");
  length = n * height * 0.66;
  printf ("4 0 %d %d 0 0 %g %g 4 %d %d %d %d %s\n", 
	scrcolor,tdepth,(tsize*13500./STDSIZ), 
	tangle*MPI/STDSIZ,
	height, length, x, y, str);
}

