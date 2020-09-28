/* Module mprintxf in program vid */

/* XFIG 3.2 driver.  

This module defines device-independent subroutines for "mvid.c".
Nothing real fancy here.  All this program will do for now is
dump out on standard output.  We can redirect at a later stage.
Let's just see it work for now.

*/

#include	<stdio.h>
#include	<string.h>
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
#define NCOLS 16

static int col_tran[NCOLS] = {0,1,2,3,4,5,6,7,
		           8,10,13,16,19,22,25,7};
                    
#define COLORBASE 100
#define NGCOLS 296
static int gcols_used[NGCOLS] = {0}; /* Check to see if color is used */
				  /*  If not, must print out pseudo-object */
				  /*  for xfig */

static float gcols[NGCOLS][3] = {/* from mprintc.c, taken from ccols.txt */

				/* bright colors */

        0.024,    0,   0.2,     /* blue - green - yellow - red */
        0.082,    0,   0.57,
        0.082,    0,    1,
        0,      0.22,    1,
        0,      0.64,   0.9,
        0,      0.9,    0.5,
        0,        1,     0,
        0.35,     1,     0,
        0.67,     1,     0,
        0.8,    0.8,    0,
        0.9,     0.65,    0,
        0.9,     0.53,    0,
        1,       0.4,    0,
        1,      0.32,    0,
        1,      0.18,    0,
        1,        0,     0,

        0.024,    0,   0.2,     /* blue - green - yellow - red, shows PSPs and spikes */
        0.082,    0,   0.57,
        0.082,    0,    1,
        0,      0.22,    1,
        0,      0.64,   0.9,
        0,      0.93,   0.47,
        0,        1,     0,
        0.35,     1,     0,
        0.67,     1,     0,
        0.98,   0.98,    0,
        1,      0.75,    0,
        1,      0.53,    0,
        1,       0.3,    0,
        1,      0.12,    0,
        1,      0.016,  0.0078,
        1,        0,     0.43,

        0,        0,  0.63,                /* blue - purple - red */
        0.031,    0,  0.86,
        0.047,    0,  0.98,
        0.078,    0,  0.94,
        0.16,     0,  0.86,
        0.24,     0,  0.78,
        0.31,     0,  0.71,
        0.39,     0,  0.63,
        0.47,     0,  0.55,
        0.55,     0,  0.47,
        0.63,     0,  0.39,
        0.71,     0,  0.31,
        0.78,     0,  0.24,
        0.86,     0,  0.16,
        0.94,     0,  0.078,
        1,     0,     0,

        0,     0,   0,          /* black - red */
        0.031, 0,   0,
        0.047, 0,   0,
        0.078, 0,   0,
        0.16,  0,   0,
        0.24,  0,   0,
        0.31,  0,   0,
        0.39,  0,   0,
        0.47,  0,   0,
        0.55,  0,   0,
        0.63,  0,   0,
        0.71,  0,   0,
        0.78,  0,   0,
        0.86,  0,   0,
        0.94,  0,   0,
        1,     0,   0,

        0,   0,     0,          /* black - green */
        0,   0.031, 0,
        0,   0.047, 0,
        0,   0.078, 0,
        0,   0.16,  0,
        0,   0.24,  0,
        0,   0.31,  0,
        0,   0.39,  0,
        0,   0.47,  0,
        0,   0.55,  0,
        0,   0.63,  0,
        0,   0.71,  0,
        0,   0.78,  0,
        0,   0.86,  0,
        0,   0.94,  0,
        0,   1,     0,

        0,   0,   0,            /* black - blue */
        0,   0,   0.031,
        0,   0,   0.047,
        0,   0,   0.078,
        0,   0,   0.16,
        0,   0,   0.24,
        0,   0,   0.31,
        0,   0,   0.39,
        0,   0,   0.47,
        0,   0,   0.55,
        0,   0,   0.63,
        0,   0,   0.71,
        0,   0,   0.78,
        0,   0,   0.86,
        0,   0,   0.94,
        0,   0,   1,

	0.00, 0.00, 0.00,	/* grayscale values */
	0.01, 0.01, 0.01,
	0.02, 0.02, 0.02,
	0.03, 0.03, 0.03,
	0.04, 0.04, 0.04,
	0.05, 0.05, 0.05,
	0.06, 0.06, 0.06,
	0.07, 0.07, 0.07,
	0.08, 0.08, 0.08,
	0.09, 0.09, 0.09,
	0.10, 0.10, 0.10,
	0.11, 0.11, 0.11,
	0.12, 0.12, 0.12,
	0.13, 0.13, 0.13,
	0.14, 0.14, 0.14,
	0.15, 0.15, 0.15,
	0.16, 0.16, 0.16,
	0.17, 0.17, 0.17,
	0.18, 0.18, 0.18,
	0.19, 0.19, 0.19,
	0.20, 0.20, 0.20,
	0.21, 0.21, 0.21,
	0.22, 0.22, 0.22,
	0.23, 0.23, 0.23,
	0.24, 0.24, 0.24,
	0.25, 0.25, 0.25,
	0.26, 0.26, 0.26,
	0.27, 0.27, 0.27,
	0.28, 0.28, 0.28,
	0.29, 0.29, 0.29,
	0.30, 0.30, 0.30,
	0.31, 0.31, 0.31,
	0.32, 0.32, 0.32,
	0.33, 0.33, 0.33,
	0.34, 0.34, 0.34,
	0.35, 0.35, 0.35,
	0.36, 0.36, 0.36,
	0.37, 0.37, 0.37,
	0.38, 0.38, 0.38,
	0.39, 0.39, 0.39,
	0.40, 0.40, 0.40,
	0.41, 0.41, 0.41,
	0.42, 0.42, 0.42,
	0.43, 0.43, 0.43,
	0.44, 0.44, 0.44,
	0.45, 0.45, 0.45,
	0.46, 0.46, 0.46,
	0.47, 0.47, 0.47,
	0.48, 0.48, 0.48,
	0.49, 0.49, 0.49,
	0.50, 0.50, 0.50,
	0.51, 0.51, 0.51,
	0.52, 0.52, 0.52,
	0.53, 0.53, 0.53,
	0.54, 0.54, 0.54,
	0.55, 0.55, 0.55,
	0.56, 0.56, 0.56,
	0.57, 0.57, 0.57,
	0.58, 0.58, 0.58,
	0.59, 0.59, 0.59,
	0.60, 0.60, 0.60,
	0.61, 0.61, 0.61,
	0.62, 0.62, 0.62,
	0.63, 0.63, 0.63,
	0.64, 0.64, 0.64,
	0.65, 0.65, 0.65,
	0.66, 0.66, 0.66,
	0.67, 0.67, 0.67,
	0.68, 0.68, 0.68,
	0.69, 0.69, 0.69,
	0.70, 0.70, 0.70,
	0.71, 0.71, 0.71,
	0.72, 0.72, 0.72,
	0.73, 0.73, 0.73,
	0.74, 0.74, 0.74,
	0.75, 0.75, 0.75,
	0.76, 0.76, 0.76,
	0.77, 0.77, 0.77,
	0.78, 0.78, 0.78,
	0.79, 0.79, 0.79,
	0.80, 0.80, 0.80,
	0.81, 0.81, 0.81,
	0.82, 0.82, 0.82,
	0.83, 0.83, 0.83,
	0.84, 0.84, 0.84,
	0.85, 0.85, 0.85,
	0.86, 0.86, 0.86,
	0.87, 0.87, 0.87,
	0.88, 0.88, 0.88,
	0.89, 0.89, 0.89,
	0.90, 0.90, 0.90,
	0.91, 0.91, 0.91,
	0.92, 0.92, 0.92,
	0.93, 0.93, 0.93,
	0.94, 0.94, 0.94,
	0.95, 0.95, 0.95,
	0.96, 0.96, 0.96,
	0.97, 0.97, 0.97,
	0.98, 0.98, 0.98,
	0.99, 0.99, 0.99,

			     /* see "ccols.txt", "convcol" and "ccols.x" */
			     /* taken from /usr/lib/X11/rgb.txt */

	0.980, 0.941, 0.902,		/* linen */ 
	0.980, 0.922, 0.843,		/* AntiqueWhite */ 
	1.000, 0.937, 0.835,		/* PapayaWhip */ 
	1.000, 0.922, 0.804,		/* BlanchedAlmond */ 
	1.000, 0.855, 0.725,		/* PeachPuff */ 
	1.000, 0.894, 0.710,		/* moccasin */ 
	1.000, 0.980, 0.804,		/* LemonChiffon */ 
	1.000, 0.961, 0.933,		/* seashell */ 
	0.941, 1.000, 0.941,		/* honeydew */ 
	0.961, 1.000, 0.980,		/* MintCream */ 
	0.941, 1.000, 1.000,		/* azure */ 
	0.902, 0.902, 0.980,		/* lavender */ 
	1.000, 0.941, 0.961,		/* LavenderBlush */ 
	1.000, 0.894, 0.882,		/* MistyRose */ 
	0.098, 0.098, 0.439,		/* MidnightBlue */ 
	0.000, 0.000, 0.502,		/* NavyBlue */ 
	0.392, 0.584, 0.929,		/* CornflowerBlue */ 
	0.282, 0.239, 0.545,		/* DarkSlateBlue */ 
	0.416, 0.353, 0.804,		/* SlateBlue */ 
	0.518, 0.439, 1.000,		/* LightSlateBlue */ 
	0.000, 0.000, 0.804,		/* MediumBlue */ 
	0.255, 0.412, 0.882,		/* RoyalBlue */ 
	0.118, 0.565, 1.000,		/* DodgerBlue */ 
	0.000, 0.749, 1.000,		/* DeepSkyBlue */ 
	0.529, 0.808, 0.922,		/* SkyBlue */ 
	0.275, 0.510, 0.706,		/* SteelBlue */ 
	0.690, 0.769, 0.871,		/* LightSteelBlue */ 
	0.678, 0.847, 0.902,		/* LightBlue */ 
	0.686, 0.933, 0.933,		/* PaleTurquoise */ 
	0.000, 0.808, 0.820,		/* DarkTurquoise */ 
	0.282, 0.820, 0.800,		/* MediumTurquoise */ 
	0.251, 0.878, 0.816,		/* turquoise */ 
	0.373, 0.620, 0.627,		/* CadetBlue */ 
	0.400, 0.804, 0.667,		/* MediumAquamarine */ 
	0.498, 1.000, 0.831,		/* aquamarine */ 
	0.333, 0.420, 0.184,		/* DarkOliveGreen */ 
	0.561, 0.737, 0.561,		/* DarkSeaGreen */ 
	0.180, 0.545, 0.341,		/* SeaGreen */ 
	0.235, 0.702, 0.443,		/* MediumSeaGreen */ 
	0.125, 0.698, 0.667,		/* LightSeaGreen */ 
	0.596, 0.984, 0.596,		/* PaleGreen */ 
	0.000, 1.000, 0.498,		/* SpringGreen */ 
	0.486, 0.988, 0.000,		/* LawnGreen */ 
	0.498, 1.000, 0.000,		/* chartreuse */ 
	0.000, 0.980, 0.604,		/* MediumSpringGreen */ 
	0.678, 1.000, 0.184,		/* GreenYellow */ 
	0.196, 0.804, 0.196,		/* LimeGreen */ 
	0.604, 0.804, 0.196,		/* YellowGreen */ 
	0.133, 0.545, 0.133,		/* ForestGreen */ 
	0.420, 0.557, 0.137,		/* OliveDrab */ 
	0.741, 0.718, 0.420,		/* DarkKhaki */ 
	0.941, 0.902, 0.549,		/* khaki */ 
	0.933, 0.910, 0.667,		/* PaleGoldenrod */ 
	0.980, 0.980, 0.824,		/* LightGoldenrodYellow */ 
	1.000, 1.000, 0.878,		/* LightYellow */ 
	1.000, 0.843, 0.000,		/* gold */ 
	0.933, 0.867, 0.510,		/* LightGoldenrod */ 
	0.855, 0.647, 0.125,		/* goldenrod */ 
	0.722, 0.525, 0.043,		/* DarkGoldenrod */ 
	0.737, 0.561, 0.561,		/* RosyBrown */ 
	0.804, 0.361, 0.361,		/* IndianRed */ 
	0.545, 0.271, 0.075,		/* SaddleBrown */ 
	0.627, 0.322, 0.176,		/* sienna */ 
	0.804, 0.522, 0.247,		/* peru */ 
	0.871, 0.722, 0.529,		/* burlywood */ 
	0.961, 0.961, 0.863,		/* beige */ 
	0.961, 0.871, 0.702,		/* wheat */ 
	0.957, 0.643, 0.376,		/* SandyBrown */ 
	0.824, 0.706, 0.549,		/* tan */ 
	0.824, 0.412, 0.118,		/* chocolate */ 
	0.698, 0.133, 0.133,		/* firebrick */ 
	0.647, 0.165, 0.165,		/* brown */ 
	0.914, 0.588, 0.478,		/* DarkSalmon */ 
	0.980, 0.502, 0.447,		/* salmon */ 
	1.000, 0.627, 0.478,		/* LightSalmon */ 
	1.000, 0.647, 0.000,		/* orange */ 
	1.000, 0.549, 0.000,		/* DarkOrange */ 
	1.000, 0.498, 0.314,		/* coral */ 
	0.941, 0.502, 0.502,		/* LightCoral */ 
	1.000, 0.388, 0.278,		/* tomato */ 
	1.000, 0.271, 0.000,		/* OrangeRed */ 
	1.000, 0.000, 0.000,		/* red */ 
	1.000, 0.412, 0.706,		/* HotPink */ 
	1.000, 0.078, 0.576,		/* DeepPink */ 
	1.000, 0.753, 0.796,		/* pink */ 
	1.000, 0.714, 0.757,		/* LightPink */ 
	0.859, 0.439, 0.576,		/* PaleVioletRed */ 
	0.690, 0.188, 0.376,		/* maroon */ 
	0.780, 0.082, 0.522,		/* MediumVioletRed */ 
	0.816, 0.125, 0.565,		/* VioletRed */ 
	1.000, 0.000, 1.000,		/* magenta */ 
	0.933, 0.510, 0.933,		/* violet */ 
	0.867, 0.627, 0.867,		/* plum */ 
	0.855, 0.439, 0.839,		/* orchid */ 
	0.729, 0.333, 0.827,		/* MediumOrchid */ 
	0.600, 0.196, 0.800,		/* DarkOrchid */ 
	0.580, 0.000, 0.827,		/* DarkViolet */ 
	0.541, 0.169, 0.886,		/* BlueViolet */ 
	0.627, 0.125, 0.941,		/* purple */ 
	0.576, 0.439, 0.859,		/* MediumPurple */ 
	};


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

#define DBUFLEN 10000
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
void xfdumppts(void);

#define RNDUP 0.001

/*--------------------------------------------*/

void txfinit(void)
{
	int i;

    printf("#FIG 3.2\n");
    printf("Landscape\n");
    printf("Center\n");
    printf("Inches\n");
    printf("Letter\n");
    printf("100.0\n");
    printf("Single\n");
    printf("-2\n");
    printf("1200 2\n");
    for (i=0; i<NGCOLS; i++) {
       printf("0 %d #%2.2x%2.2x%2.2x\n",i+COLORBASE+NCOLS,
		       		     (int)(gcols[i][0]*255+RNDUP),
		       		     (int)(gcols[i][1]*255+RNDUP),
				     (int)(gcols[i][2]*255+RNDUP));

    }
    printf("6 %d %d %d %d\n",0+xoffset,-yoffset,xhi+xoffset,yhi-yoffset);
}

void txfexit(void)
{
    xfdumppts();
    printf("-6\n");
}
void clrcxf(void)
{}
void xyconxf(void)
{}

void xyvidxf(void)
{}

void xfdot(void)
{}

int setxfcol(int pen)
{
   double grayval;
   static int oldcolor;

  oldcolor = scrcolor;
  if ((backgr<0 || backgr==WHITE) && (pen == WHITE)) scrcolor = BLACK;
  else if ((backgr==BLACK) && (pen == BLACK)) scrcolor = WHITE;
  else {
     if (pen<16) scrcolor = col_tran[pen];
     else {
	     if (pen >= NGCOLS) pen = NGCOLS-1;
	     scrcolor = pen + COLORBASE;

     }
  }
  return (oldcolor);
}

void xfdumppts()

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

void movexf(int ix, int iy)
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

void drawxf(int ix, int iy)
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

void txfdrcirc(int ix, int iy, int irad, int fill)
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
	(int)x,(int)y,(int)rad,(int)rad,(int)x,(int)y,(int)(x+rad),(int)(y+rad));

 xpos = x;
 ypos = y;
}

void txfdrrect(int ix1, int iy1, int ix2, int iy2, int ix3, int iy3, int ix4, int iy4, int fill)
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

void txfdrtri(int ix1, int iy1, int ix2, int iy2, int ix3, int iy3, int fill)
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

void txffill(int fill)

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

int txfxsize(void)
{ return STDSIZ; }

int txfysize(void)
{ return STDSIZ; }

void putstrxf(char *str)
{
   static int oldtsize=0,oldtfont=0;
   int height, length, tdepth;
   int x,y,n;
   // char *strcat(char *dest, const char *src);

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

