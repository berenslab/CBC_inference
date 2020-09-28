
/***************************************************/
/*                                                 */
/*                     MPRINTX                     */
/*                                                 */
/*                                                 */ 
/*                     Wei  Han                    */
/*                                                 */
/*                     July, 90                    */
/*                                                 */ 
/***************************************************/

/* Modifications by Andrew the Hsu 11/1990 */

#include <X11/Xlib.h>
#include <X11/Xutil.h>
/* #include <X11/Xos.h>  */
/* #include <X11/Xatom.h> */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "icon_bitmap"
#include "colors.h"
#define BITMAPDEPTH 1

#define XSCALE 18
#define YSCALE 18
#define XHHI 4000
#define XHI 1000
#define YHI 1000 
#define STDSIZ 16383.

static Display *display;
static unsigned long foreground_pixel, background_pixel, border_pixel;
static char *progname;
static int screen_num;
static Window win;
static GC gc;
static XFontStruct *font_info;
static unsigned int width, height;
static XSizeHints size_hints;
static int window_size = 0;
static Screen *screen_ptr;
static Pixmap pix;
static Colormap default_cmap;
static char colornam[40] = {0}; 
static char *color_name[] = {
 /*       "black","magenta","green","cyan","red","orange","brown","white",
        "grey","LightBlue","aquamarine","CadetBlue","firebrick","blue",
        "yellow","SpringGreen"};
*/
        "black","DodgerBlue","green","cyan","red","magenta","peru","white",
        "gray60","LightBlue","aquamarine","CadetBlue","firebrick",
	"VioletRed","yellow", "snow"};

#define NVCOLS 16
#define NNCOLS 32
#define NPCOLS 16
#define NRCOLS 16
#define NDCOLS 16
#define NBCOLS 16

/*
double vcols[NVCOLS][3] = { .02, .02,  1.0,	/* red-yel-grn-blue, pastel */
/*			    .1,  .2,   .9,
			    .1,  .35,  .75,
			    .2,  .4,   .7,
			    .2,  .5,   .6,
			    .2,  .7,   .5,
			    .2,  .8,   .4,
			    .2,  .9,   .3,
			    .3,  .9,   .2,
			    .5,  .9,   .2,
			    .6,  .95,  .2,
			     1,  .7,   .1,
			     1,  .6,   .1,
			     1,  .5,   .1,
			     1,  .4,   .2,
			     1,  .4,    1};

*/

double vcols[NVCOLS][3] = { 				/* red - yel - grn - blue */
				0.024,    0,   0.2,
				0.082,    0,   0.57,
				0.082,    0,    1,
   				0,      0.22,    1,
   				0,      0.64,   0.9,
   				0,      0.9,    0.5,
   				0,        1,     0,
				0.35,     1,     0,
				0.67,     1,     0,
				0.83,   0.83,    0,
   				0.9,    0.65,    0,
   				0.9,    0.53,    0,
   				1,       0.4,    0,
   				1,      0.32,    0,
   				1,      0.18,   0,
   				1,        0,     0,
			  };

double ncols[NNCOLS][3] = { 				/* red - yel - grn - blue, shows PSPs and spikes */
				0.024,    0,   0.2,
				0.052,    0,   0.4,
				0.062,    0,   0.6,
				0.082,  0.05,  0.8,
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
			  };

double pcols[NPCOLS][3] = {				/* red - purple - blue */ 
				   0,     0,  0.63,
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
			  };

double rcols[NRCOLS][3] = {				/* black - red */ 
				   0,     0,  0,
				0.031,    0,  0,
				0.047,    0,  0,
				0.078,    0,  0,
				0.16,     0,  0,
				0.24,     0,  0,
				0.31,     0,  0,
				0.39,     0,  0,
				0.47,     0,  0,
				0.55,     0,  0,
				0.63,     0,  0,
				0.71,     0,  0,
				0.78,     0,  0,
				0.86,     0,  0,
				0.94,     0,  0,
   				1,     .15,   .15,
			  };

double dcols[NDCOLS][3] = {				/* black - green */ 
				0,   0,      0,
				0,   0.031,  0,
				0,   0.047,  0,
				0,   0.078,  0,
				0,   0.16,   0,
				0,   0.24,   0,
				0,   0.31,   0,
				0,   0.39,   0,
				0,   0.47,   0,
				0,   0.55,   0,
				0,   0.63,   0,
				0,   0.71,   0,
				0,   0.78,   0,
				0,   0.86,   0,
				0,   0.94,   0,
   				0,   1,      0
			  };

double bcols[NBCOLS][3] = {				/* black - blue */ 
				0,   0,   0,    
				0,   0,   0.031,
				0,   0,   0.047,
				0,   0,   0.078,
				0,   0,   0.16, 
				0,   0,   0.24, 
				0,   0,   0.31, 
				0,   0,   0.39, 
				0,   0,   0.47, 
				0,   0,   0.60, 
				0,   0,   0.67, 
				0,   0,   0.80, 
				0,   0,   0.92, 
				0.05, 0.05, 0.95, 
				0.1,  0.1, 0.98, 
   				0.15, 0.15, 1,   
			  };

#define NGCOLS 100

#define NCCOLS 100
double ccols[NCCOLS][3] = { 
			     /* see "ccols.txt", "convcol", and "ccols.x" */
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

static XColor exact_def;

#define NCOLORS 1024
static unsigned int col_val[NCOLORS];

static int xpos, ypos;
static double xscale=XSCALE;
static double yscale=XSCALE;
static int yhi=YHI;
static int xoffset = 5;
static int yoffset = 5;

#define DBUFLEN 10000
static XPoint dbuf[DBUFLEN] = {0};
static int npts=0;
static int fillval=0;
static int ncolor=0;

#define CWIDTH  (.017 * STDSIZ)       /* char matrix size for enh graph */
#define CHEIGHT (.04 * STDSIZ)
static int cwidth=CWIDTH;
static int cheight=CHEIGHT;

extern char *display_name;
extern int scrcolor;
extern int backgr;
static int backgrcol=BLACK;
static int foregrcol=WHITE;
extern int scrrev;
extern double scrsiz;
extern int pargc;
extern char **pargv;
extern int multfl;
static int gmodefl=0;
static char rname[]={"vid"};
static char *gargc[2]={rname, NULL};
int gmodex(int argc, char **argv);
void xdumppts(void);
void tmodex(void);
void xyvidx(int x, int y);
void movex(int x, int y);
void drawx(int x, int y);
void setorg(int x,int y);
void linsegx(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
void load_font(XFontStruct **font_info);
void get_GC(void);

/*-----------------------------------------*/

void txinit(void)
{

  if (!gmodefl || multfl) {
    gmodefl = 1;
    gmodex(1,gargc);
  }
}

void txexit(void)
{
  xdumppts();
  gmodefl = 0;
  tmodex();
}

void xyconx(int x, int y)
{
  int xloc,yloc;

#define XMAX 80
#define YMAX 27
 x = (x<XMAX? x : XMAX);
 y = (y<YMAX? y : YMAX);
 x = (x<0? 0 : x);
 y = (y<0? 0 : y);
 xloc = x * cwidth;
 yloc = y * cheight;
 xyvidx(xloc,yloc);
}

void xyvidx(int x, int y)
{
movex (x,y);
setorg(x,y);
}

void putstrx (char *str)
{
  int width,len;

   len = strlen(str);
   width=XTextWidth(font_info,str,len);
   XDrawString(display,pix,gc,xpos+xoffset,ypos+yoffset,str,len);
   XDrawString(display,win,gc,xpos+xoffset,ypos+yoffset,str,len);
   XFlush(display);
}

void xdot(int xd, int yd)
{
 extern int scrrev;

 xdumppts();
 xd/=xscale;
 yd/=yscale;
 if  (!scrrev) yd=yhi-yd;
 linsegx(xd,yd,xd,yd);
 xpos = xd;
 ypos = yd;
}
 
int setxcol(int pen)
{
  static int oldcol=0, temp;

  temp = oldcol;
  XSetForeground(display,gc,col_val[pen]);
  if (backgr == WHITE && pen==WHITE) pen = BLACK;
  scrcolor = pen;
  oldcol = pen;
  return temp;
}

void xdumppts(void)

{
 int i;

 
 XSetForeground(display,gc,col_val[ncolor]);
 if (npts>1) {

  for (i=0; i<(npts-1); i++) {
      XDrawLine(display, pix, gc, 
		(int)dbuf[i].x,(int)dbuf[i].y,
		(int)dbuf[i+1].x,(int)dbuf[i+1].y);
      XDrawLine(display, win, gc, 
		(int)dbuf[i].x,(int)dbuf[i].y,
                (int)dbuf[i+1].x,(int)dbuf[i+1].y);
   }

   if (fillval) {
     XFillPolygon(display, pix, gc, dbuf, npts, Nonconvex, CoordModeOrigin);
     XFillPolygon(display, win, gc, dbuf, npts, Nonconvex, CoordModeOrigin);
    }
   npts = 0;
  }
 XSetForeground(display,gc,col_val[scrcolor]);
}
 
void movex(int x, int y)
{
 extern int scrrev;

 xdumppts();
 x/=xscale;
 y/=yscale;
 xpos = x;
 ypos = y;
 if (!scrrev) ypos = yhi-ypos;

 npts = 0;
 dbuf[npts].x = xpos+xoffset;
 dbuf[npts].y = ypos+yoffset;
 if (++npts >= DBUFLEN) {
        fprintf (stderr,"movex: too many points in one polygon\n");
        npts= DBUFLEN;
 }

}
 
void drawx(int x, int y)
{
 extern int scrrrev;

 x/=xscale;
 y/=yscale;
 if  (!scrrev) y=yhi-y;

 if (x>XHHI) x=XHHI;
 if (x<0) x=0;
 if (y>yhi) y=yhi;
 if (y<0) y=0;
 if (xpos >XHHI) xpos=XHHI;
 if (xpos <0) xpos=0;
 if (ypos >yhi) ypos=yhi;
 if (ypos <0) ypos=0;
 xpos=x;
 ypos=y;

 ncolor = scrcolor;
 dbuf[npts].x = x+xoffset;
 dbuf[npts].y = y+yoffset;
 if (++npts >= DBUFLEN) {
        fprintf (stderr,"drawx: too many points in one polygon\n");
        npts= DBUFLEN;
 }
}

double xscalex (double size)
{
   double temp;

   temp = xscale;
   xscale = STDSIZ / size;
   return(temp);
}

double yscalex (double size)
{
   double temp;

   temp = yscale;
   yscale = STDSIZ / size;
   yhi = size;
   return(temp);
}

/**********************************************/

int txxsize(void)
{ return STDSIZ; }

/**********************************************/

int txysize(void)
{ return STDSIZ; }

/**********************************************/

#define SMALL 1
#define OK 0
#define COLMAX 255.0

int gmodex(int argc, char **argv)
{
     int x = 0, y = 0;            /* window position */
     unsigned int border_width = 4; /* border fourn pixels wide */
     unsigned int display_width,display_height; 
     char *window_name = "XVid";
     char *icon_name = "XVid";
     Pixmap icon_pixmap;
     XEvent report;
     int i, totcolors = 0;
     int count = 1;
     double winscal = .999;
     double g;

     /* connect to X server */

     progname = argv[0];
     if ((display=XOpenDisplay(display_name))==NULL)
     {
	(void) fprintf(stderr,
		  "basicwin:cannot connect to X server %s\n",
		  XDisplayName(display_name));
        exit(-1);
     }
     /* get screen size form display structure macro */
     screen_num=DefaultScreen(display);
     screen_ptr=DefaultScreenOfDisplay(display);

     display_width=DisplayWidth(display,screen_num);
     display_height=DisplayHeight(display,screen_num);
     
     /*size window with enough room for text */
     height= display_height * 0.5 * scrsiz;
     width=  height;
     xscalex(((double) width)*winscal);
     yscalex(((double) height)*winscal);

 /* fprintf (stderr, "width == %d, height == %d\n", width, height); /* */

/**************************/
/* Allocates colors for the colormap */

        for(i = 0; i < NCOLORS; i++) {
                col_val[i] = i % 128;		/* set colors roughly first */
        }
        default_cmap = DefaultColormap(display,screen_num);
        for(i = 0; i < NBCOLS; i++) {
/*              fprintf(stderr,"allocating color %s\n",color_name[i]); */
                if(!XParseColor(display,default_cmap,color_name[i],&exact_def)) {
                     fprintf(stderr,"mprintx: can't parse %s\n",color_name[i]);
                     /* exit(0); */
                }

                if (!XAllocColor(display, default_cmap,&exact_def)) {
                  fprintf(stderr,"mprintx: can't allocate '%s'\n",color_name[i]);
                        /* exit(0); */
                }

                col_val[i+totcolors] = exact_def.pixel;
        }
        totcolors += NBCOLS;
        for(i = 0; i < NVCOLS; i++) {
		sprintf (colornam,"rgbi:%g/%g/%g", vcols[i][0], 
						   vcols[i][1], 
						   vcols[i][2]);
                if(!XParseColor(display,default_cmap,colornam,&exact_def)) {
                     fprintf(stderr,"mprintx: can't parse %s\n",colornam);
                     /* exit(0); */
                }

                if (!XAllocColor(display, default_cmap,&exact_def)) {
                  fprintf(stderr,"mprintx: can't allocate '%s'\n",colornam);
                        /* exit(0); */
                }

                col_val[i+totcolors] = exact_def.pixel;
        }
        totcolors += NVCOLS;
        for(i = 0; i < NNCOLS; i++) {
		sprintf (colornam,"rgbi:%g/%g/%g", ncols[i][0], 
						   ncols[i][1], 
						   ncols[i][2]);
                if(!XParseColor(display,default_cmap,colornam,&exact_def)) {
                     fprintf(stderr,"mprintx: can't parse %s\n",colornam);
                     /* exit(0); */
                }

                if (!XAllocColor(display, default_cmap,&exact_def)) {
                  fprintf(stderr,"mprintx: can't allocate '%s'\n",colornam);
                        /* exit(0); */
                }

                col_val[i+totcolors] = exact_def.pixel;
        }
        totcolors += NNCOLS;
        for(i = 0; i < NPCOLS; i++) {
		sprintf (colornam,"rgbi:%g/%g/%g", pcols[i][0], 
						   pcols[i][1], 
						   pcols[i][2]);
                if(!XParseColor(display,default_cmap,colornam,&exact_def)) {
                     fprintf(stderr,"mprintx: can't parse %s\n",colornam);
                     /* exit(0); */
                }

                if (!XAllocColor(display, default_cmap,&exact_def)) {
                  fprintf(stderr,"mprintx: can't allocate '%s'\n",colornam);
                        /* exit(0); */
                }

                col_val[i+totcolors] = exact_def.pixel;
        }
        totcolors += NPCOLS;
        for(i = 0; i < NRCOLS; i++) {
		sprintf (colornam,"rgbi:%g/%g/%g", rcols[i][0], 
						   rcols[i][1], 
						   rcols[i][2]);
                if(!XParseColor(display,default_cmap,colornam,&exact_def)) {
                     fprintf(stderr,"mprintx: can't parse %s\n",colornam);
                     /* exit(0); */
                }

                if (!XAllocColor(display, default_cmap,&exact_def)) {
                  fprintf(stderr,"mprintx: can't allocate '%s'\n",colornam);
                        /* exit(0); */
                }

                col_val[i+totcolors] = exact_def.pixel;
        }
        totcolors += NRCOLS;
        for(i = 0; i < NDCOLS; i++) {
		sprintf (colornam,"rgbi:%g/%g/%g", dcols[i][0], 
						   dcols[i][1], 
						   dcols[i][2]);
                if(!XParseColor(display,default_cmap,colornam,&exact_def)) {
                     fprintf(stderr,"mprintx: can't parse %s\n",colornam);
                     /* exit(0); */
                }

                if (!XAllocColor(display, default_cmap,&exact_def)) {
                  fprintf(stderr,"mprintx: can't allocate '%s'\n",colornam);
                        /* exit(0); */
                }

                col_val[i+totcolors] = exact_def.pixel;
        }
        totcolors += NDCOLS;
        for(i = 0; i < NBCOLS; i++) {
		sprintf (colornam,"rgbi:%g/%g/%g", bcols[i][0], 
						   bcols[i][1], 
						   bcols[i][2]);
                if(!XParseColor(display,default_cmap,colornam,&exact_def)) {
                     fprintf(stderr,"mprintx: can't parse %s\n",colornam);
                     /* exit(0); */
                }

                if (!XAllocColor(display, default_cmap,&exact_def)) {
                  fprintf(stderr,"mprintx: can't allocate '%s'\n",colornam);
                        /* exit(0); */
                }

                col_val[i+totcolors] = exact_def.pixel;
        }
        totcolors += NBCOLS;
        for(i = 0; i < NGCOLS; i++) {	/* grayscale */
		g = ((double)i) / NGCOLS;
		g = pow(g,2.0) + 0.03;	/* make simple gamma curve */
		g *= 0.97;
		sprintf (colornam,"rgbi:%g/%g/%g",g,g,g); 
                if(!XParseColor(display,default_cmap,colornam,&exact_def)) {
                     fprintf(stderr,"mprintx: can't parse %s\n",colornam);
                     /* exit(0); */
                }

                if (!XAllocColor(display, default_cmap,&exact_def)) {
                  fprintf(stderr,"mprintx: can't allocate '%s'\n",colornam);
                        /* exit(0); */
                }
                col_val[i+totcolors] = exact_def.pixel;
        }
        totcolors += NGCOLS;
        for(i = 0; i < NCCOLS; i++) {
		sprintf (colornam,"rgbi:%g/%g/%g", ccols[i][0], 
						   ccols[i][1], 
						   ccols[i][2]);
                if(!XParseColor(display,default_cmap,colornam,&exact_def)) {
                     fprintf(stderr,"mprintx: can't parse %s\n",colornam);
                     /* exit(0); */
                }

                if (!XAllocColor(display, default_cmap,&exact_def)) {
                  fprintf(stderr,"mprintx: can't allocate '%s'\n",colornam);
                        /* exit(0); */
                }
                col_val[i+totcolors] = exact_def.pixel;
        }
/* fprintf(stderr,"Colors allocated successfully\n"); */

     if (backgr >= 0) {
       backgrcol=backgr;
     }

     /* create opaque window */
     win=XCreateSimpleWindow(display,RootWindow(display,screen_num),
	     x,y,width,height,border_width,
	     col_val[foregrcol],
	     col_val[backgrcol]);
	
     /* Create Pixmap to display graphs */
     pix = XCreatePixmap(display,
			RootWindow(display,screen_num),
			width,
			height,
			DefaultDepth(display,screen_num));

     /* create pixmap of depth 1 (bitmap) for icon */
     icon_pixmap = XCreateBitmapFromData(display, win, icon_bitmap_bits,
		icon_bitmap_width, icon_bitmap_height);

     /* initialize size hint property for window manager */
/*
     size_hints.flags=PPosition |PSize|PMinSize;
     size_hints.x=x;
     size_hints.y=y;
     size_hints.width=2*width;
     size_hints.height=2*height;
     size_hints.min_width=1000;
     size_hints.min_height=1000;
*/

/* Improved size hints */
     size_hints.flags=PPosition |PSize|PMinSize;
     size_hints.x=200;
     size_hints.y=200;
     size_hints.width=2*width;
     size_hints.height=2*height;
     size_hints.min_width=400*scrsiz;
     size_hints.min_height=400*scrsiz;

     /* set properties for window manager (always before mapping) */
     XSetStandardProperties(display,win,window_name,icon_name,
	  icon_pixmap, argv, argc, &size_hints);

     /* select event types wanted */
     /* XSelectInput(display, win, ExposureMask|KeyPressMask |
	      ButtonPressMask | StructureNotifyMask); */
      XSelectInput(display, win, ExposureMask| 
  	      ButtonPressMask | StructureNotifyMask); 

     load_font(&font_info);

     /* create GC for text and drawing */
     get_GC();

    /* Now fill the Pixmaps with Black since some servers default
                to a White background */
    setxcol(backgrcol);
    XFillRectangle(display, pix, gc, 0, 0, width, height);


     /* Display window */
     XMapWindow(display, win);

/********************************/
        report.type = 0;
    while (count-- > 0) {

 	XNextEvent(display,&report);
	switch(report.type) {
	case Expose:
		XClearWindow(display,win);
/*
     		draw_text(width,height);
     		draw_graphics(width, height);
*/
		XCopyArea(display,pix,win,gc,0,0,width,height,0,0);
           break;
        case ConfigureNotify:
            width = report.xconfigure.width;
	    height = report.xconfigure.height;
	    xscalex((double) width*winscal);
	    yscalex((double) height*winscal);
	    if ((width< size_hints.min_width) ||
		 ( height< size_hints.min_height))
               window_size=SMALL;
            else
	       window_size=OK;
            break;
         default:
	    break;
         }
    }
}
/*******************************/
void dealwithwindow(int parm)
{
XEvent	report;
int count = 1;
int newevent;

/*  The trick here is to make sure that the process that called this
    subroutine sleeps when the display is done but is beeing
    viewed by the user. The window continues to need display
    events because it may be moved and covered by other windows.
    If we don't sleep, then this process takes up cpu time,
    polling for new events.
*/

  
do {
     static int last=0;

   if (parm) {	/* Finished displaying window. wait here for user to look */

      if (!last) { 
        xdumppts();		/* finish drawing last polygon */
        last = 1;
      }
      XMaskEvent(display, ExposureMask| 
	ButtonPressMask | StructureNotifyMask, &report); 
      newevent = 1;

    }
    else {   /* parm == 0, when displaying moves and draws */
      count--;
      newevent = XCheckMaskEvent(display, ExposureMask | 
		ButtonPressMask | StructureNotifyMask, &report); 
   }

   if (newevent) {
        switch(report.type) {
        case ButtonPress:
        case Expose:
		XClearWindow(display,win);
		XCopyArea(display,pix,win,gc,0,0,width,height,0,0);
		XFlush(display);
           break;
        case ConfigureNotify:
            width = report.xconfigure.width;
            height = report.xconfigure.height;
            if ((width< size_hints.min_width) ||
                 ( height< size_hints.min_height))
               window_size=SMALL;
            else
               window_size=OK;
            break;
         default:
            break;
         }
	}
     } while (count > 0); 
}


/*******************************/
void tmodex(void)
{
	     getchar();
             XUnloadFont(display, font_info->fid); 
	     XFreeGC(display, gc);
             XCloseDisplay(display); 
}

 
/********************************/
void get_GC(void)
{

   unsigned long valuemask = 0; /* ignore XGCvalues and use defaults */
   XGCValues values;
   unsigned int line_width = 1;
   int line_style = LineSolid;
   int cap_style = CapRound;
   int join_style = JoinRound;
   int dash_offset =0;
   static char dash_list[] = {
	12,24   };
   int list_length = 2;

   /* Create default graphics contect */
   gc = XCreateGC (display, win, valuemask, &values);

   /* specify font */
   XSetFont(display, gc, font_info->fid);

   XSetForeground(display,gc, col_val[foregrcol]); 

   /* set line attributes */
   XSetLineAttributes(display,gc,line_width,line_style,cap_style,join_style);

   /* set dashes to be line_width in length */
   XSetDashes(display, gc,dash_offset,dash_list,list_length);

}

/***************************/
void load_font(XFontStruct **font_info)
{
   char *fontname= "fixed";

   /* access font */
   if ((*font_info=XLoadQueryFont(display,fontname)) == NULL)
   {
      (void) fprintf(stderr, "Basic: Cannot open fixed font\n");
      exit(-1);
   }

} 

/***************************/
void draw_text(unsigned int win_width, unsigned int win_height)
{
   int y=20;   /* offset form corner of window */
   char *string1 = "Vid";
   char *string2 = "X11 version";
   char *string3 = "Feb 1994";
   int len1, len2, len3;
   int width1, width2, width3;
   char cd_height[50], cd_width[50], cd_depth[50];

   /* need length for both XTextWidth and XDrawString */
   len1 = strlen(string1);
   len2 = strlen(string2);
   len3 = strlen(string3);
   
   /* get string widths for centering */
   width1=XTextWidth(font_info,string1,len1);
   width2=XTextWidth(font_info,string2,len2);
   width3=XTextWidth(font_info,string3,len3);

   /* output text, centered on each line */
   XDrawString(display,pix,gc,(win_width - width1)/2,y,string1,len1);
   XDrawString(display,pix,gc,(win_width - width2)/2,
	     (int)(win_height -30), string2,len2);
   XDrawString(display, pix, gc,(win_width-width3)/2,
             (int)(win_height -15), string3,len3);
}
/***********************/

void draw_graphics(unsigned int window_width, unsigned int window_height)
{
   int x,y;
   unsigned int width, height;

   height =window_height;
   width = window_width;
   x = window_width/2 - width/2;
   y = window_height/2 - height/2;

   XDrawRectangle(display, pix, gc, x, y, width, height);

   while (x <= window_width/2+width/2-10) {
   x+=10;
   linsegx(x, window_height/2+height/2-10, x, window_height/2+height/2);
   }
   x = window_width/2 - width/2;
   while (x <= window_width/2+width/2-100) {
   x+=100;
   linsegx(x, window_height/2+height/2-20, x, window_height/2+height/2);
   }
   while (y <= window_height/2+height/2-10) {
   y+=10;
   linsegx(window_width/2-width/2, y, window_width/2-width/2+10, y);
   }
   y = window_height/2 - height/2;
   while (y <= window_height/2+height/2-100) {
   y+=100;
   linsegx(window_width/2-width/2, y, window_width/2-width/2+20, y);

   }
}
/**************************/
void linsegx(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
   XDrawLine(display, pix, gc, x1+xoffset, y1+yoffset, x2, y2);
   XDrawLine(display, win, gc, x1+xoffset, y1+yoffset, x2, y2);
   XFlush(display);
}


/*************************/
void polygon(int x, int y, int radius, int vertex)
            /*center*/ 
                   
{
    int index=0;
    unsigned int x1=0, y1=0, x2=0, y2=0; 
    double angle;

    angle=3.1415926*2/vertex;
    x1=x+(int)(radius*sin(angle/2));
    y1=y-(int)(radius*cos(angle/2));
    while (index++ !=vertex) {
       x2=x1+(int)(radius*sin(angle/2)*2*cos(index*angle));
       y2=y1+(int)(radius*sin(angle/2)*2*sin(index*angle));
       linsegx (x1, y1, x2, y2);
       x1=x2; y1=y2;
    }
}
/*************************/
void xerase(void)
{
   int oldcol;

  xdumppts();
  oldcol = setxcol(backgrcol);			/* set to background for erase */
  XFillRectangle(display, pix, gc, 0, 0, width, height); 
  XClearWindow(display, win);
/* 
  draw_text(width, height);
  draw_graphics(width, height);
*/
  XFlush(display); 
  oldcol = setxcol(oldcol);			/* set the color back again */
  xyconx (0,24);				/* home up cursor */
}
/*************************/

void xcircle(int x, int y, int radius, int fill)
                          
/* Draw a circle at location (x,y) with "radius".
*/

{
 xdumppts();
 x/=xscale;
 y/=yscale;
 radius/=yscale;
 if  (!scrrev) y=yhi-y;

 if (x>XHHI) x=XHHI;
 if (x<0) x=0;
 if (y>yhi) y=yhi;
 if (y<0) y=0;

 x += xoffset - radius;
 y += yoffset - radius;  

  if (radius < 1) {
    XDrawPoint(display, pix, gc, x, y);
    XDrawPoint(display, win, gc, x, y);
  }
  else if (fill) {
    XFillArc(display, pix, gc, x, y, 2*radius, 2*radius, 0, 64*360);
    XFillArc(display, win, gc, x, y, 2*radius, 2*radius, 0, 64*360);
  }  /*  */
  else {
    XDrawArc(display, pix, gc, x, y, 2*radius, 2*radius, 0, 64*360); 
    XDrawArc(display, win, gc, x, y, 2*radius, 2*radius, 0, 64*360); /* */
  }
   XFlush(display);
}

/*************************/

void xellipse(int x, int y, int radius1, int radius2, int fill)
                          
/* Draw an ellipse at location (x,y) with major and minor axes.
*/

{
 xdumppts();
 x/=xscale;
 y/=yscale;
 radius1/=yscale;
 radius2/=yscale;
 if  (!scrrev) y=yhi-y;

 if (x>XHHI) x=XHHI;
 if (x<0) x=0;
 if (y>yhi) y=yhi;
 if (y<0) y=0;

 x += xoffset - radius1;
 y += yoffset - radius1;  

  if (radius1 < 1) {
    XDrawPoint(display, pix, gc, x, y);
    XDrawPoint(display, win, gc, x, y);
  }
  else if (fill) {
    XFillArc(display, pix, gc, x, y, 2*radius1, 2*radius2, 0, 64*360);
    XFillArc(display, win, gc, x, y, 2*radius1, 2*radius2, 0, 64*360);
  }  /*  */
  else {
    XDrawArc(display, pix, gc, x, y, 2*radius1, 2*radius2, 0, 64*360); 
    XDrawArc(display, win, gc, x, y, 2*radius1, 2*radius2, 0, 64*360); /* */
  }
   XFlush(display);
}
/*************************/

void xrect (int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, int fill)
{
 xdumppts();
 x1 /= xscale;
 y1 /= yscale;
 x2 /= xscale;
 y2 /= yscale;
 x3 /= xscale;
 y3 /= yscale;
 x4 /= xscale;
 y4 /= yscale;
 if  (!scrrev) y1=yhi-y1;
 if  (!scrrev) y2=yhi-y2;
 if  (!scrrev) y3=yhi-y3;
 if  (!scrrev) y4=yhi-y4;

 x1 += xoffset;
 y1 += yoffset;
 x2 += xoffset;
 y2 += yoffset;
 x3 += xoffset;
 y3 += yoffset;
 x4 += xoffset;
 y4 += yoffset;

  XDrawLine(display, pix, gc, x1, y1, x2, y2);          /* draw rectangle */
  XDrawLine(display, pix, gc, x2, y2, x3, y3);
  XDrawLine(display, pix, gc, x3, y3, x4, y4);
  XDrawLine(display, pix, gc, x4, y4, x1, y1);

  XDrawLine(display, win, gc, x1, y1, x2, y2);          /* draw rectangle */
  XDrawLine(display, win, gc, x2, y2, x3, y3);
  XDrawLine(display, win, gc, x3, y3, x4, y4);
  XDrawLine(display, win, gc, x4, y4, x1, y1);

 if (fill) {

#define NUMXPT 4

    static XPoint xpnt[NUMXPT] = {0};

   xpnt[0].x = x1;
   xpnt[0].y = y1;
   xpnt[1].x = x2;
   xpnt[1].y = y2;
   xpnt[2].x = x3;
   xpnt[2].y = y3;
   xpnt[3].x = x4;
   xpnt[3].y = y4;

   XFillPolygon(display, pix, gc, xpnt, NUMXPT, Convex, CoordModeOrigin);
   XFillPolygon(display, win, gc, xpnt, NUMXPT, Convex, CoordModeOrigin);


 }
 XFlush(display);
}

/*************************/

void xtri (int x1, int y1, int x2, int y2, int x3, int y3, int fill)
{
 xdumppts();
 x1 /= xscale;
 y1 /= yscale;
 x2 /= xscale;
 y2 /= yscale;
 x3 /= xscale;
 y3 /= yscale;
 if  (!scrrev) y1=yhi-y1;
 if  (!scrrev) y2=yhi-y2;
 if  (!scrrev) y3=yhi-y3;

 x1 += xoffset;
 y1 += yoffset;
 x2 += xoffset;
 y2 += yoffset;
 x3 += xoffset;
 y3 += yoffset;

  XDrawLine(display, pix, gc, x1, y1, x2, y2);          /* draw triangle */
  XDrawLine(display, pix, gc, x2, y2, x3, y3);
  XDrawLine(display, pix, gc, x3, y3, x1, y1);

  XDrawLine(display, win, gc, x1, y1, x2, y2);          /* draw triangle */
  XDrawLine(display, win, gc, x2, y2, x3, y3);
  XDrawLine(display, win, gc, x3, y3, x1, y1);

 if (fill) {

#define NUMXPTT 3

    static XPoint xpntx[NUMXPTT] = {0};

   xpntx[0].x = x1;
   xpntx[0].y = y1;
   xpntx[1].x = x2;
   xpntx[1].y = y2;
   xpntx[2].x = x3;
   xpntx[2].y = y3;

//fprintf (stderr,"xtri %d\n",fill);

   XFillPolygon(display, pix, gc, xpntx, NUMXPTT, Convex, CoordModeOrigin);
   XFillPolygon(display, win, gc, xpntx, NUMXPTT, Convex, CoordModeOrigin);


 }
 XFlush(display);
}

void txfill(int val)
{
 fillval = val;
}

/*************************/
/*EOF */

