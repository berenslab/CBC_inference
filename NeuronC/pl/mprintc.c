/* Module mprintc in program vid */

/* Color Post-script device driver, makes separate pages. */

/* This module with  device-independent subroutines in "mvid.c".  */

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<time.h>
#include	<stdio.h>
#include	<sys/types.h>
/* #include	<sys/timeb.h> */

#include	"colors.h"

/* Define image to be 7.5 x 7.5; User can imbed postscript commands to shrink
   image to any size if need */

#define XSCALE 30.339    /* STDSIZ/(7.5*72); LW is 72 dpi */
#define YSCALE 30.339
#define XHHI 550
#define XHI 550
#define YHI 675

#define STDSIZ 16383.
#define GRKFONT 'g'	/* greek font for 's' command */
#define NCOLORS 256

#define NPSCOLS 328
double pscols[NPSCOLS][3] = {
	0, 0, 0,		/* black */
	.1, .4, 1, 		/* blue */
	0, 1, 0, 		/* green */
	0, 1, 1, 		/* cyan */
	1, 0, 0, 		/* red */
	1, 0, 1, 		/* magenta */
	.8, .6, 0, 		/* brown (dark yel) */
	1, 1, 1, 		/* white */
	.6, .6, .6, 		/* gray */
	.4, .8, 1, 		/* light blue */
	.6, 1, .6, 		/* light green */
	.6, 1, 1, 		/* light cyan */
	1, .6, .6, 		/* light red */
	1, .6, 1, 		/* light magenta */
	1, 1, 0, 		/* yellow */
	.9, .9, .9,		/* bright white */

        0.024,    0,   0.2,	/* blue - green - yellow - red */
        0.082,    0,   0.57,
        0.082,    0,    1,
        0,      0.22,    1,
        0,      0.64,   0.9,
        0,      0.9,   0.5,
        0,        1,     0,
        0.35,     1,     0,
        0.67,     1,     0,
        0.8,    0.8,    0,
        0.9,    0.65,    0,
        0.9,    0.53,    0,
        1,       0.4,    0,
        1,      0.32,    0,
        1,      0.018,   0,
        1,        0,     0,

	0.024,    0,   0.2,		/* blue - green - yellow - red, shows PSPs and spikes */
	0.052,    0,   0.4,
	0.062,    0,    0.6,
	0.082,  0.05,   0.8,
	0.04,   0.22,    1,
	0,      0.32,  0.95,
	0,      0.38,   0.92,
	0,      0.45,   0.9,
	0,      0.54,   0.8,
	0,      0.64,   0.75,
	0,      0.74,   0.7,
	0,      0.85,   0.5,
	0,      0.90,   0.3,
	0,      0.95,   0.1,
	0,        1,     0,
	0.05,     1,     0,
	0.10,     1,     0,
	0.15,     1,     0,
	0.20,     1,     0,
	0.25,     1,     0,
	0.30,  0.99,     0,
	0.35,  0.98,     0,
	0.45,  0.97,     0,
	0.55,  0.95,    0,
	0.65,  0.90,    0,
	0.75,  0.85,    0,
	0.80,  0.83,    0,
	0.82,  0.82,    0,
	0.95,  0.75,    0,
	1,     0.53,    0,
	1,      0.3,    0,
	1,        0,     0,

        0,     0,  0.63,		/* blue - purple - red */
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

        0,     0,   0,		/* black - red */
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

        0,   0,     0,		/* black - green */
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

        0,   0,   0,    	/* black - blue */
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
static double yhi = YHI;

extern char tfont;
extern unsigned int tsize,tangle;
extern char *comment;
extern int scrrev;
extern int backgr;
extern int rotflag;
extern int pargc;
extern char **pargv;
extern char *pagefilenam;
char filenam;

static int seqnum = 0;
static int fontfl = 1;
FILE *filout = NULL;
int gmodec(int argc, char **argv);

/*----------------------------------------------*/

int tcinit(void)
{
gmodec(pargc,pargv);
}
void tcexit(void)
{
   fprintf(filout,"stroke showpage\n");
   fprintf(filout,"restore\n");
/*   fprintf(filout,"clear\n"); */
}
int tcpage(void)
{
  gmodec(pargc,pargv);
}

void clrcc(void)
{}
void xyconc(void)
{}

void xyvidc(void)
{}

void cdot(void)
{}

int setccol(int pen)
{
   double grayval;
   static int oldcolor,scrcolor;

  if ((backgr < 0 || backgr==WHITE) && pen==WHITE) pen = BLACK;
  oldcolor = scrcolor;
  scrcolor = pen;

  fprintf(filout,"stroke\n");
  if (pen > NPSCOLS) pen = NPSCOLS; 
  else if (pen < 0) pen = 0;
  fprintf(filout,"%g %g %g setrgbcolor\n",     pscols[pen][0],
					pscols[pen][1],
					pscols[pen][2]);
  return (oldcolor);
}

void movec(int ix, int iy)
{
 extern int scrrev;

 xpos = ix/xscale;
 ypos = iy/yscale;
 if (xpos >XHHI) xpos=XHHI;
 if (xpos <0) xpos=0;
 if (ypos >yhi) ypos=yhi;
 if (ypos <0) ypos=0;
 if (scrrev) ypos = yhi-ypos;
 fprintf(filout,"s\n%-5.4g %-5.4g m\n",xpos,ypos);
}

void drawc(int ix, int iy)
{
 extern int scrrrev;
 double x,y;

 x = ix/xscale;
 y = iy/yscale;
 if  (scrrev) y=yhi-y;
 if (x>XHHI) x=XHHI;
 if (x<0) x=0;
 if (y>yhi) y=yhi;
 if (y<0) y=0;
 fprintf(filout,"%-5.4g %-5.4g l\n",x,y);
 xpos=x;
 ypos=y;
}

void tcdrcirc(int ix, int iy, int irad, int fill)
{
 double rad,x,y;

 x = ix/yscale;
 y = iy/yscale;
 rad = irad/yscale;
 if  (scrrev) y=yhi-y;
 if (x>XHHI) x=XHHI;
 if (x<0) x=0;
 if (y>yhi) y=yhi;
 if (y<0) y=0;
 fprintf(filout,"s %-5.4g %-5.4g %-5.4g 0 360 arc",x,y,rad);
 if (!fill) fprintf(filout," s\n");
 else       fprintf(filout," fill\n");
 movec(ix,iy);
 xpos = x;
 ypos = y;
}

void tcdrrect(int ix1, int iy1, int ix2, int iy2, int ix3, int iy3, int ix4, int iy4, int fill)
{
 extern int scrrev;
 double x1,y1,x2,y2,x3,y3,x4,y4;
 int color;

 x1 = ix1/xscale;
 y1 = iy1/yscale;
 x2 = ix2/xscale;
 y2 = iy2/yscale;
 x3 = ix3/xscale;
 y3 = iy3/yscale;
 x4 = ix4/xscale;
 y4 = iy4/yscale;

 if  (scrrev) {
   y1=yhi-y1;
   y2=yhi-y2;
   y3=yhi-y3;
   y4=yhi-y4;
 }

 if (!fill && backgr>=0) {
   color = setccol(backgr);
   fprintf(filout,"%-5.4g %-5.4g m\n",x1,y1);
   fprintf(filout,"%-5.4g %-5.4g l\n",x2,y2);
   fprintf(filout,"%-5.4g %-5.4g l\n",x3,y3);
   fprintf(filout,"%-5.4g %-5.4g l\n",x4,y4);
   fprintf(filout,"%-5.4g %-5.4g l\n",x1,y1);
   fprintf(filout," fill\n");
   setccol(color); 
 }
 fprintf(filout,"%-5.4g %-5.4g m\n",x1,y1);
 fprintf(filout,"%-5.4g %-5.4g l\n",x2,y2);
 fprintf(filout,"%-5.4g %-5.4g l\n",x3,y3);
 fprintf(filout,"%-5.4g %-5.4g l\n",x4,y4);
 fprintf(filout,"%-5.4g %-5.4g l\n",x1,y1);
 if (!fill) fprintf(filout,"s\n");
 else       fprintf(filout,"fill\n");
 movec(ix1,iy1);
 xpos = x1;
 ypos = y1;
}

void tcdrtri(int ix1, int iy1, int ix2, int iy2, int ix3, int iy3, int fill)
{
 extern int scrrev;
 double x1,y1,x2,y2,x3,y3;
 int color;

 x1 = ix1/xscale;
 y1 = iy1/yscale;
 x2 = ix2/xscale;
 y2 = iy2/yscale;
 x3 = ix3/xscale;
 y3 = iy3/yscale;

 if  (scrrev) {
   y1=yhi-y1;
   y2=yhi-y2;
   y3=yhi-y3;
 }

 if (!fill && backgr>=0) {
   color = setccol(backgr);
   fprintf(filout,"%-5.4g %-5.4g m\n",x1,y1);
   fprintf(filout,"%-5.4g %-5.4g l\n",x2,y2);
   fprintf(filout,"%-5.4g %-5.4g l\n",x3,y3);
   fprintf(filout,"%-5.4g %-5.4g l\n",x1,y1);
   fprintf(filout," fill\n");
   setccol (color);
 }
 fprintf(filout,"%-5.4g %-5.4g m\n",x1,y1);
 fprintf(filout,"%-5.4g %-5.4g l\n",x2,y2);
 fprintf(filout,"%-5.4g %-5.4g l\n",x3,y3);
 fprintf(filout,"%-5.4g %-5.4g l\n",x1,y1);
 if (!fill) fprintf(filout," s\n");
 else      fprintf(filout," fill\n");
 movec(ix1,iy1);
 xpos = x1;
 ypos = y1;
}

void tcfill(int fill)

{
  if (fill) fprintf(filout,"fill\n");
}

double xscalec (double size)
{ 
return ((double)xscale);
}

double yscalec (double size)
{ 
return ((double)yscale);
}

int tcxsize(void)
{ return STDSIZ; }

int tcysize(void)
{ return STDSIZ; }

void putstrc(char *str)
{
   static int oldtsize=0,oldtfont=0;

  if (tangle) fprintf(filout,"%-5.3g rotate\n",(tangle*360./32768.));
  if (tsize!=oldtsize || oldtfont!=tfont || fontfl) {
    if (tsize==0) tsize=10;
    if (tfont == GRKFONT) fprintf(filout,"/Symbol findfont\n");
    else                  fprintf(filout,"/Times-Roman findfont\n");
    fprintf(filout,"%g scalefont setfont\n",(tsize*13500./STDSIZ));
    fontfl = 0;
  }
  fprintf(filout,"(%s) show s\n",str);
  if (tangle) fprintf(filout,"-%-5.3g rotate\n",(tangle*360./32768.));
  oldtsize=tsize;
  oldtfont=tfont;
}


char *makefilenam(char *filnam)

/* make the next file name from a template */

{
  static char sbuf[80];

  sprintf (sbuf,"%s%04d.ps",filnam,++seqnum);
  return sbuf;
};

void openfilout (char *filnam)

{
   if (filout && filout!=stdout) {
      tcexit();
      fflush(filout);
      fclose(filout);
   }
   fontfl = 1;
   filnam = makefilenam(filnam);
   if (!(filout = fopen(filnam,"w"))) {
      fprintf (stderr,"vid, mprintp can't open file %s\n",filnam);
      exit(0);
   }
}

int gmodec(int argc, char **argv)
{
/* Variables declared here.... */
   int	x,y;
   long	clock;
   char	date[30];

    time(&clock);
    strcpy(date, ctime(&clock));
    if (pagefilenam) openfilout(pagefilenam);
    else filout = stdout;
 
    fprintf(filout,"%%!PS-Adobe1.1\n");
    fprintf(filout,"%%%%Creator: Vid - Display routine by R.G. Smith and");
    fprintf(filout," A. Hsu.\n");
    fprintf(filout,"%%%%Title: %s\n",argv[0]);
    fprintf(filout,"%%%%CreationDate: %s",&date[0]);
    fprintf(filout,"%%%%Pages: 1\n");
    fprintf(filout,"%%%%BoundingBox: 50 108 600 680\n");
    fprintf(filout,"%%%%EndComments\n");
    fprintf(filout,"%%%%BeginProlog\n");
    fprintf(filout,"/m {moveto} def\n");
    fprintf(filout,"/l {lineto} def\n");
    fprintf(filout,"/s {stroke} def\n");
    fprintf(filout,"%%%%EndProlog\n");
    fprintf(filout,"%%%%Page: 0 1\n");
    fprintf(filout,"\nnewpath\n");
    fprintf(filout,"save\n");
    fprintf(filout,"%% This is for text at the top of page..\n");
    fprintf(filout,"/CommentFont\n  /Helvetica findfont 10 scalefont def\n");
    fprintf(filout,"/Date (%s) def\n", &date[0]);
/*    if (comment)  fprintf(filout,"/Comment (%s) def\n",comment);
    else          fprintf(filout,"/Comment ([         ]) def\n");
*/
    /* Let's try drawing an 8 x 8 box for a frame...*/
    fprintf(filout,"%% Border Box Begins here\n");
    fprintf(filout,"%% 0 0 0 setrgbcolor\n");
    fprintf(filout,"%% 0 0 moveto\n");
    fprintf(filout,"%% 650 0 rlineto\n");
    fprintf(filout,"%% 0 800 rlineto\n");
    fprintf(filout,"%% -650 0 rlineto\n");
    fprintf(filout,"%% closepath\n");
    fprintf(filout,"%% fill\n");
    fprintf(filout,"%% 1 setlinewidth\n");

if (backgr>=0) {
      int pen;

    /* the equivalent of "setccol(backgr)" -- but can't use it here. */

    pen = backgr;
    fprintf(filout,"stroke\n");
    if (pen > NPSCOLS) pen = NPSCOLS; 
    else if (pen < 0) pen = 0;
    fprintf(filout,"%g %g %g setrgbcolor\n",	pscols[pen][0],
					pscols[pen][1],
					pscols[pen][2]);
    fprintf(filout,"\nnewpath\n");
    fprintf(filout,"0 0 moveto\n");
    fprintf(filout,"650 0 rlineto\n");
    fprintf(filout,"0 800 rlineto\n");
    fprintf(filout,"-650 0 rlineto\n");
    fprintf(filout,"closepath\n");
    fprintf(filout,"fill\n\n");
}

    fprintf(filout,"%% This is where you can scale, rotate and translate the image\n");
    fprintf(filout,"%%     if you feel the need to manipulate the image further...\n");
    if (rotflag) {
      fprintf(filout,"%6.3f rotate\n",rotflag*90.0);
      fprintf(filout,"30 -580  translate\n");
    }
    else {
      fprintf(filout,"0 rotate\n");
      fprintf(filout,"50 108 translate\n");
    }
    fprintf(filout,"1 1 scale\n");
    fprintf(filout,"0 560 moveto\n");
    fprintf(filout,"CommentFont setfont\n");
    fprintf(filout,"Date show\n");
/*    fprintf(filout,"0 -72 moveto\n");
    fprintf(filout,"Comment show\n");
*/

    fprintf(filout,"%% Vid data begins here...\n");
    fprintf(filout,".5 setlinewidth\n");
}
