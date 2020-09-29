/* Segment ncplot in Program nc */

/* Plots neuronal circuits */

/*      May 88                  R.G. Smith */

extern "C" {

#include <stdio.h>
#include <string.h>
#include <math.h>
char *fgets(char *s, int size, FILE *stream);
#ifdef CPML
#include <cpml.h>
#endif

}


#include "nc.h"
#include "ndef.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncelem.h"
#include "ncomp.h"
#include "ncplot.h"
#include "control.h"
#include "colors.h"
#include "gprim.h"
#include "ncio.h"

#define DEBUG

#ifdef DEBUG
#include "ncdebug.h"
#endif

extern plotfr plotnod[PLOTNODSIZ];      /* node numbers, zero sample voltage */

extern int plotnum;
int listplots[PLOTNODSIZ] = {0};	/* list of command line forced plots */
int leftaxes[PLOTNODSIZ] = {0};		/* list of command line left axes */
int rightaxes[PLOTNODSIZ] = {0};	/* list of command line right axes */
int aplots[PLOTNODSIZ] = {0};		/* assigned (logical) plots */
int nplots[PLOTNODSIZ] = {0};		/* number of traces so far in plot */
int nnplots[PLOTNODSIZ] = {0};		/* total named traces in plot */
int plotdone[PLOTNODSIZ] = {0};		/* plots axis has been done */
int plotlabl[PLOTNODSIZ] = {0};		/* plot has label */
int naplots = 0;			/* number of assigned plots */

double splots[PLOTNODSIZ] = {0};	/* list of plot sizes */
double scplots[PLOTNODSIZ] = {0};	/* cum list of plot sizes */
double yaplot[PLOTNODSIZ] = {0};	/* y axis, indexed by aplot */
double ybplot[PLOTNODSIZ] = {0};	/* y base, indexed by aplot */
double yplcharsiz[PLOTNODSIZ] = {0};	/* y trace label size */

int numlaxis=0;
int numraxis=0;

static int dttplots = 0;		/* total plots read from input file */
static double xplotmargl;		/* X position of left y axis */
static double xplotmargr;		/* margin width to right of X axis */
static int numyaxis;			/* number of Y axes */
extern double charsiz;			/* char size for plot labels */
const char *rootframe = "/root";	/* name of root plot frame */

double plotxmax;			/* xmax for plot, endexp or setxmax */
double plotxmin;
extern double ssetxmax;			/* command-line overrides */
extern double ssetxmin;

extern int nonodes;                    	/* =1 -> no node numbers in display */
extern int rowmode;                    	/* =1 -> numbers for plot on row */
extern int nofilfl;			/* =1 ->  means no filename on plot */
extern int finfl;			/* =1 ->  means final plot, no labels */
extern int setnox;			/* =1 ->  means no x values */
extern int yprecision;			/* precision of y values in printout */
int yfieldwidth;			/* total field width of y values */

extern int setvid;                      /* vidmode set on */
extern int unsetvid;                    /* vidmode set off */
extern int ncerror;                     /* error flag set by "warning()" */
extern int tabmode;                     /* print tabs on output */
extern int setplsep;                    /* set separate plots on output */
extern int disp_ray;			/* make ray-traced output for povray */

extern char *infile;                    /* input filename pointer */
extern char *srcfile;                   /* source filename pointer */
extern char *plotlabel;                 /* label for plot */
extern char *stfile;                    /* stimulus filename pointer */
static char stfilnam[60]={0};		/* stimulus file name */
static char infilnam[60]={0};		/* input file name */
static char srcfilnam[60]={0};		/* source file name */
static char plotlabelnam[100]={0};	/* label name */
static char complin[400]={0};		/* "#c xx comps" line from infile*/

int pinitfl = 0;                 	/* == 1 means plotinit() run */

#ifdef __cplusplus
extern "C" {
#endif

#include "stdplt.h"
#include "gr.h"
double exp(double), log(double);
void exit(int code);
double strtod(const char *, char **);
char *strncpy(char *dest, const char *src, size_t n);
char *strncat(char *dest, const char *src, size_t n);

#ifdef __cplusplus
}
#endif

int notinit(double var);
int notinit(double var);
int notinit(int   var);

double ncabs(double x);
int looklist(int i);
int lookaxes(int i, int list[]);
void prtext (const char *adj, double csize, const char *str);
void prnum (const char *adj,double csize, const char *fmt,
		double val1,double val2, double val3);
void drsymb (double csize, int ch);
void plotinit (int totplots);
void autograph(FILE *instream, char *inbuf, int numpl);
void drawframe (int nplot, int totplots, int useframe);
char *emalloc(unsigned int n);

double recallx (int ymax, int ypos, int xpos, double *arr);
double recally (int ymax, int ypos, int xpos, double *arr);
void storeplot(double yval, double xval, 
		int ymax, int ypos, int xpos, double *arr);
double synfilt(lpfilt *lpnt, double transrate);
lpfilt *maksynfilt(int stype, int nfilt,double offset,double tfall,double val,
		double *timec, double timestep);
lpfilt *makbessfilt(unsigned int options, int order, double alpha);
void rprintf (char *fmt,double x, double y, double z, double size, int color);

void plotpen(int pen, int nplot);
char *gfrname();
int getaplot (int i);

void drchar (double csize, int ch); 			/* char */
void drcirc (double csize, int fill); 			/* circle */
void drsqr  (double csize, int fill);			/* square */
void drrect (double csize, double aspect, int fill);	/* rect */
void drtri  (double csize, int fill);			/* tri */
void drtriu (double csize, int fill);			/* ups tri */


char *print_version(double version);
double scanversion(char *vbuf);

/*------------------------------------*/


static double yplotmarg=0.0;

static char colortab[8] ={BLUE,GREEN,CYAN,RED,MAGENTA,BROWN,WHITE,YELLOW};
static char color2  [8] ={LTBLUE,LTGREEN,LTCYAN,LTRED,LTMAG,BROWN,WHITE,YELLOW};
static char color3  [8] ={BLUE,GREEN,CYAN,RED,MAGENTA,BROWN,WHITE,YELLOW};
static char color4  [8] ={LTBLUE,LTGREEN,LTCYAN,LTRED,LTMAG,BROWN,WHITE,YELLOW};
static char color5  [8] ={BLUE,GREEN,CYAN,RED,MAGENTA,BROWN,WHITE,YELLOW};
static char color6  [8] ={LTBLUE,LTGREEN,LTCYAN,LTRED,LTMAG,BROWN,WHITE,YELLOW};
static char color7  [8] ={BLUE,GREEN,CYAN,RED,MAGENTA,BROWN,WHITE,YELLOW};
static char color8  [8] ={LTBLUE,LTGREEN,LTCYAN,LTRED,LTMAG,BROWN,WHITE,YELLOW};
static char color9  [8] ={BLUE,GREEN,CYAN,RED,MAGENTA,BROWN,WHITE,YELLOW};
static char color10 [8] ={LTBLUE,LTGREEN,LTCYAN,LTRED,LTMAG,BROWN,WHITE,YELLOW};

/*------------------------------------*/

int findplfr (int p, int tot)

{
   int i,n,found;

   for (found=i=n=0; i<=tot; i++) {
     if (strcmp(gfrname(),plotnod[i].plframe)==0) {
        if (n++ == p) {
           found=1;
           break;
        }
     }
  }
  if (found==1) return i;
  else          return 0; 
}

/*------------------------------------*/

void mplot(double yval, double xval, int mtotplots, int plotpos, int pflag)

/* plot a value from one or more nodes on graph */

{
   double x,y,xaxis,yrange,xrange,yaxis,ybase,csize;
   double *arr;
   int n,nn,pen,plpos,plposn,aplot;
   int ch;
   char *plframe;
   lpfilt *filt;
   Symbol *varpen;
#define NUMBUFSIZ 50
	     char buf[NUMBUFSIZ];
#undef NUMBUFSIZ

#ifdef DEBUG
  if (debug & NCPLOT) ncfprintf (stderr,"mplot  ");
#endif

#define XABSMARGL  .03
#define XPLOTMARGL (XABSMARGL + .05)
#define XPLOTMARGR .06
#define XPLOTLABEL .02
#define YPLOTMARG  .03
#define YAXISSPC  .028
#define YAXISWIDTH .1
#define YNODEWIDTH .02
#define YLABELPOS   .035
#define YLABELPOSN  .05
#define LARGETIC   .015
#define SMALLTIC   .007
#define LABELTIC   .016
#define LABELOFFS  .3
#define LABELAX    .02

 if (!pinitfl) return;                   /* return if graph not initialized */

 plpos = looklist(plotpos); 		   /* possibly modify plot order */
 if (plpos > dttplots-1)
      plpos = 0;			   /* use default ranges and maxmins */

 if ((filt=plotnod[plpos].filt) != (lpfilt*)NULL) {
    yval = synfilt(filt,yval);
 }

 if ((arr=plotnod[plpos].arr) != NULL) {	/* plot into predefined array */
	  int arrindex;
     if ((arrindex=plotnod[plpos].arrindex) < plotnod[plpos].maxindex*2) {
          arr[arrindex++] = xval;
          arr[arrindex++] = yval;
          plotnod[plpos].arrindex = arrindex;
     }
#ifdef DEBUG 
 if (debug & NCPLOT) ncfprintf (stderr,"mplotend0\n");
#endif
  return;
 }

 if (!vidmode) {                         /* if no graphics, print numbers */

  if (ncabs(xval) < 1e-12) xval = 0;
  //fprintf (stderr, "setxmin %g\n", setxmin);
  //fprintf (stderr, "setxmax %g\n", setxmax);
   if (!((!notinit(setxmin) && (xval < setxmin)) || 
     (!notinit(setxmax) && (xval > setxmax)))) {
     if (tabmode) {
      if (plotpos==0) ncfprintf (stdout,"%-8.8g\t",xval);
      sprintf (buf,"%-*2$.*3$g",yval,yfieldwidth,yprecision);
      ncfprintf (stdout,"%s\t",buf);
    } else {
      if (plotpos==0) ncfprintf (stdout,"%-8.8g ",xval);
      sprintf (buf,"%-*2$.*3$g",yval,yfieldwidth,yprecision);
      ncfprintf (stdout,"%s  ",buf);
    }
    if (plotpos == (mtotplots-1)) ncfprintf (stdout,"\n");
   }
  fflush (stdout);
#ifdef DEBUG 
 if (debug & NCPLOT) ncfprintf (stderr,"mplotend1\n");
#endif
  return;
}

 if (plsep) mtotplots = naplots;

 if (mtotplots == 0) mtotplots = 1;
 if (!plsep) mtotplots = 1;                /* if separate plots */

 if (strcmp(gfrname(),rootframe) != 0)     /* if diff than root fr */
    plposn = findplfr(plpos,plotnum);
 else plposn = plpos;

 //plposn = plpos;

 // nn = plpos+plposn;
 nn = plpos;

 if (nn >=PLOTNODSIZ) nn= PLOTNODSIZ-1;

 yrange = plotnod[plposn].yrange;
 xrange = plotnod[plposn].xrange;
 if (yrange == 0) yrange = 1e-20;
 if (xrange == 0) xrange = 1e-20;
 y  = yval - plotnod[plposn].pymin; 
 y /= yrange;                     		/* normalize y to max of 1 */
 xaxis = 1 - (xplotmargl + xplotmargr);   	/* calc size of plot x axis */
 x = xval  - plotnod[plposn].pxmin;
 x *= xaxis / xrange;
 x += xplotmargl;

 if (plsep) {
          int pln;

        aplot = getaplot(plotpos);
        yaxis = yaplot[aplot];
        ybase = ybplot[aplot];
        y *= yaxis;                             /* make y fit into screen */
        y += ybase;                             /* y offset */
        pen = plotnod[plposn].ppen;
	if (pen==NULLVAL) pen = colortab[plposn];
        if (pen > 0) cpen (pen);
 }
 else {                         /* draw colors, not separate plots */
        yaxis = (1.0 - yplotmarg) - YPLOTMARG * numyaxis;
        ybase = yplotmarg;
        y *= yaxis;                       /* make y fit into screen */
        y += ybase;                       /* no offset (too confusing!) */
        n = plposn; 
        pen = plotnod[n].ppen;
        if (pen == NULLVAL) pen = colortab[n];
        if (pen > 0) cpen (pen);
/*      if (dashfl) dash ((long)(pen & 7)); /* */
 }
 x = max(-.5,x); 
 x = min(2.5,x); 
 y = max(-.2,y); 
 y = min(1.2,y); 
 if (plotnod[plposn].charfl==' ')	/* reset if turned off */
     plotnod[plposn].charmode = PLINES;

 if (pen <= 0) return;  		/* do nothing for negative pen */
 csize = plotnod[plposn].csize;
 if (csize < 0) {
     ncfprintf (stderr, "mplot: negative size %g\n", csize);
     return;
 }
 ch = plotnod[plposn].charfl;
 plframe=plotnod[plposn].plframe;
 if (*plframe) frame (plframe);
 if (!((!notinit(plotnod[plposn].pxmin) && (xval < plotnod[plposn].pxmin)) || 
     (!notinit(plotnod[plposn].pxmax) && (xval > plotnod[plposn].pxmax+RNDUP)))) {
   switch (plotnod[plposn].charmode) {

    case NOLINES:				/* plot chars with no lines */
       move (x,y);
       drsymb (csize,ch);
       break;

    case   LINES:				/* plot chars with lines */
       if (!notinit(plotnod[nn].oldx)) {
         move (plotnod[nn].oldx,plotnod[nn].oldy);  /* plot lines */
         if ((plotnod[nn].oldx > x) && !allowrev) move (x,y);
         else                                     draw (x,y);
         move (plotnod[nn].oldx,plotnod[nn].oldy);  /* prev point */
         drsymb (csize,ch);		/* draw symbol at previous point */
         move (x,y);
         drsymb (csize,ch);		/* draw symbol at current point */
       }
       break;

    case   PLINES:				/* plot just lines */
       if (!notinit(plotnod[nn].oldx)) {
          /* if (!dashfl) move (plotnod[nn].oldx,plotnod[nn].oldy);*/
          move (plotnod[nn].oldx,plotnod[nn].oldy);
       if ((plotnod[nn].oldx > x) && !allowrev) move (x,y);
       else                                          draw (x,y);
       }
       else move (x,y);				/* move to first point */
       break;

    default:
       break;

  } /* switch */

  plotnod[nn].olddx = plotnod[nn].oldx; /* use plot order because */
  plotnod[nn].olddy = plotnod[nn].oldy; /* plots may be duplicated */
  plotnod[nn].oldx = x;
  plotnod[nn].oldy = y;
  if (pflag) gpurge();
 } /* if (x > ) */
 if (*plframe) frame ("..");

#ifdef DEBUG 
 if (debug & NCPLOT) ncfprintf (stderr,"mplotend2\n");
#endif
}

/*------------------------------------*/

void drsymb (double csize, int ch)

/* draw a symbol at the data point */

{
  if (ch > '@') {
      drchar (csize, ch);
  }
  else 
    switch (ch) {
      case '!':   drcirc (csize, 0); break;		/* open   circle */
      case '@':   drcirc (csize, 1); break;		/* filled circle */
      case '#':   drsqr (csize, 0); break;		/* open   square */
      case '&':   drsqr (csize, 1); break;		/* filled square */
      case '^':   drrect (csize, 3.0, 0);break;		/* open   rect */
      case '~':   drrect (csize, 3.0, 1);break;		/* filled rect */
      case '_':   drrect (csize, 0.33, 0);break;	/* open   rect */
      case '=':   drrect (csize, 0.33, 1);break;	/* filled rect */
      case '(':   drtri (csize, 0);break;		/* open       tri */
      case ')':   drtri (csize, 1);break;		/* filled     tri */
      case '{':   drtriu (csize, 0);break;		/* open   ups tri */
      case '}':   drtriu (csize, 1);break;		/* filled ups tri */
      default:    drchar (csize, ch);break;
  }
}

/*------------------------------------*/

void disperr(void)
{
    static int displayed=0;

 if (vidmode && !displayed) {
   displayed = 1;
   cpen (WHITE);
   move (0.9, 0.9);
   cwidth (.025,".");
   prtext ("l",0.025,"Error\n");
   cwidth (charsiz,".");
 }
}

/*------------------------------------*/

double geticinc(double range)
{
   int itic;
   double dtic, ticinc, srange;

#define TFACTOR 3.0

  range = ncabs(range);
  srange = range / TFACTOR; 
  dtic = log(srange) / LN10;		/* take log10 */
  itic = (int)floor(dtic);
  ticinc = exp ((double)itic * LN10);
  if (range / ticinc < 5) ticinc *= .5;
  return (ticinc); 
}

/*------------------------------------*/

void initpl(int charind)
{
   int i;

 for (i=0; i<PLOTNODSIZ; i++) {
  plotnod[i].cnod1 = 0;
  plotnod[i].cnod2 = 0;
  plotnod[i].cnod3 = 0;
  plotnod[i].cnod4 = 0;
  plotnod[i].pxmax  = LARGENUM;
  plotnod[i].pxmin  = LARGENUM;
  plotnod[i].pymax  = LARGENUM;
  plotnod[i].pymin  = LARGENUM;
  plotnod[i].olddx  = LARGENUM;
  plotnod[i].olddy  = LARGENUM;
  plotnod[i].oldx   = LARGENUM;
  plotnod[i].oldy   = LARGENUM;
  plotnod[i].yval   = 0;
  plotnod[i].xrange = 0;
  plotnod[i].yrange = 0;
  plotnod[i].pmod   = VREC;
  plotnod[i].pmod2  = 0;
  plotnod[i].pval   = 0;
  plotnod[i].ppen   = NULLVAL;
  plotnod[i].plotn  = NULLVAL;
  plotnod[i].plotsiz= 1;
  plotnod[i].plotval= LARGENUM;
  plotnod[i].plotval2=LARGENUM;
  plotnod[i].plotval3=LARGENUM;
  plotnod[i].plotval4=LARGENUM;
  plotnod[i].plotval5=LARGENUM;
  plotnod[i].plotval6=LARGENUM;
  plotnod[i].plotval7=LARGENUM;
  plotnod[i].vpen   = (int *)NULL;
  plotnod[i].vpenc  = (int *)NULL;
  plotnod[i].func   = (int *)NULL;
  plotnod[i].funcc  = (int *)NULL;
  plotnod[i].var    = (double *)NULL;
  plotnod[i].arr    = (double *)NULL;
  plotnod[i].arrindex = 0;
  plotnod[i].maxindex = 0;
  plotnod[i].filt   = (lpfilt*)NULL;
  plotnod[i].plframe[0]  = 0;
  plotnod[i].plname[0]  = 0;
  if (i>=charind) plotnod[i].csize  = 0;
  plotnod[i].automax = 0;
  if (plotnod[i].charmode < LINES || plotnod[i].charmode > PLINES) {
    plotnod[i].charfl = 0;
    plotnod[i].charmode  = PLINES;
  }
 }
 plotnum = -1;
 frame (rootframe); 
}

/*------------------------------------*/

void plotpen(int pen, int nplot)
{
 if (nplot>=PLOTNODSIZ) nplot = PLOTNODSIZ-1;
 if (nplot < 0) nplot = 0;
 plotnod[nplot].ppen = pen;

 if (! vidmode) {
    ncfprintf (stdout,"#dcolor %d plot %d\n",pen,nplot);
 }

}

/*------------------------------------*/

void plotvpen(Symbol *vpen, int nplot)
{
 if (nplot>=PLOTNODSIZ) nplot = PLOTNODSIZ-1;
 if (nplot < 0) nplot = 0;
 plotnod[nplot].vpen = (int *)vpen;
}

/*------------------------------------*/

void plotvpenc(double (*vpen)(int,double,double), int nplot)
{
 if (nplot>=PLOTNODSIZ) nplot = PLOTNODSIZ-1;
 if (nplot < 0) nplot = 0;
 plotnod[nplot].vpenc = (int *)vpen;  // C++
}

/*------------------------------------*/

void plotvpenc(double (*vpen)(int,double,double))
{
 if (plotnum>=PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
 if (plotnum < 0) plotnum = 0;
 plotnod[plotnum].vpenc = (int *)vpen;  // C++
}

/*------------------------------------*/

void plotfunc(Symbol *func, int nplot)
{
 if (nplot>=PLOTNODSIZ) nplot = PLOTNODSIZ-1;
 if (nplot < 0) nplot = 0;
 plotnod[nplot].func = (int *)func;
}

/*------------------------------------*/

void plotfunc(double (*func)(double, double))
{
 if (plotnum>=PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
 if (plotnum < 0) plotnum = 0;
 plotnod[plotnum].funcc = (int *)func;  // C++
}

/*------------------------------------*/

void plotfunc(double (*func)(double, double, double))
{
 if (plotnum>=PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
 if (plotnum < 0) plotnum = 0;
 plotnod[plotnum].funcc = (int *)func;  // C++
}

/*------------------------------------*/

void plotfunc(double (*func)(double, double, double, double))
{
 if (plotnum>=PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
 if (plotnum < 0) plotnum = 0;
 plotnod[plotnum].funcc = (int *)func;  // C++
}

/*------------------------------------*/

void plotfunc(double (*func)(double, double, double, double, double))
{
 if (plotnum>=PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
 if (plotnum < 0) plotnum = 0;
 plotnod[plotnum].funcc = (int *)func;  // C++
}

/*------------------------------------*/

void plotfunc(double (*func)(double, double, double, double, double, double))
{
 if (plotnum>=PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
 if (plotnum < 0) plotnum = 0;
 plotnod[plotnum].funcc = (int *)func;  // C++
}

/*------------------------------------*/

void plotfunc(double (*func)(double, double, double, double, double, double, double))
{
 if (plotnum>=PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
 if (plotnum < 0) plotnum = 0;
 plotnod[plotnum].funcc = (int *)func;  // C++
}

/*------------------------------------*/

void plotfunc(double (*func)(double, double, double, double, double, double, double, double))
{
 if (plotnum>=PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
 if (plotnum < 0) plotnum = 0;
 plotnod[plotnum].funcc = (int *)func;  // C++
}

/*------------------------------------*/

void plotvar(double *var)
{
 if (plotnum>=PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
 if (plotnum < 0) plotnum = 0;
 plotnod[plotnum].var = var;  // C++
}

/*------------------------------------*/

void plotarr(double *arr, int maxindex)
{
 if (plotnum>=PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
 if (plotnum < 0) plotnum = 0;
 plotnod[plotnum].arr = arr;  // C++
 plotnod[plotnum].maxindex = maxindex;
}

/*------------------------------------*/

void plotcsiz(double size, int nplot)
{
 if (nplot>=PLOTNODSIZ) nplot = PLOTNODSIZ-1;
 plotnod[nplot].csize = size;

				/* don't allow changing char during plot: */
/* if (! vidmode) {
    ncfprintf (stdout,"#zcsize %g plot %d\n",size,nplot);
 }			/* */

}

/*------------------------------------*/

void plotcsiz(double size)
{
 if (plotnum>=PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
 plotnod[plotnum].csize = size;

				/* don't allow changing char during plot: */
/* if (! vidmode) {
    ncfprintf (stdout,"#zcsize %g plot %d\n",size,nplot);
 }			/* */

}

/*------------------------------------*/

void plotchar(int lchar, int mode, int nplot)

{
 if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
 if (lchar==' ') 				/* space erases plot char */
   plotnod[nplot].charmode = PLINES; 
 else {						/* set on plot char */
   plotnod[nplot].charfl = (char)lchar;
   plotnod[nplot].charmode = mode; 
 }

				/* don't allow changing char during plot: */
/*
 if (! vidmode) {
    ncfprintf (stdout,"#lchar '%c' mode %c plot %d\n",lchar,mode,nplot);
 }			/*  */

}

/*------------------------------------*/

void plotfilt(int nfilt, double *timec,  int nplot)

{
   if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
   plotnod[nplot].filt = maksynfilt(LP, nfilt, 0.0,0.0,0.0,timec,ploti); 
}

/*------------------------------------*/

void plotfilt(int nfilt, double *timec)

{
   if (plotnum>PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
   plotnod[plotnum].filt = maksynfilt(LP, nfilt, 0.0,0.0,0.0,timec,ploti); 
}

/*------------------------------------*/

void plotfilt(int nfilt, double *timec, double initval)

{
   if (plotnum>PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
   plotnod[plotnum].filt = maksynfilt(LP, nfilt, 0.0,0.0,initval,timec,ploti); 
}

/*------------------------------------*/

void plotfilt(int nfilt, double *timec, double initval, int nplot)

{
   if (plotnum>PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
   plotnod[nplot].filt = maksynfilt(LP, nfilt, 0.0,0.0,initval,timec,ploti); 
}

/*------------------------------------*/

void plotbessfilt(double cutoff)

/* Make a fourth-order lowpass Bessel filter using matched * z-transform */
/*  To make filter with bilinear transform method, leave out opt_z */
/* see: http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html */
/* code from: http://www-users.cs.york.ac.uk/~fisher/software/mkfilter/current  */

{
   unsigned int options;
   int order;
   double alpha;

   if (plotnum>PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
   options= (opt_be | opt_lp | opt_z);
   // options= (opt_be | opt_lp );
   plotnod[plotnum].filt = makbessfilt(options,order=4,alpha=ploti*cutoff); 
}

/*------------------------------------*/

void plotname(const char *plname,  int nplot)

{
   if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
   strncpy(plotnod[nplot].plname,plname,PNAMSIZ-1); 
}

/*------------------------------------*/

void plotname(const char *plname)

{
   if (plotnum>PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
   strncpy(plotnod[plotnum].plname,plname,PNAMSIZ-1); 
}

/*------------------------------------*/

void plotn(int plotn,  int nplot)

{
   if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
   plotnod[nplot].plotn = plotn; 
}

/*------------------------------------*/

void plotn(int plotn)

{
   if (plotnum>PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
   plotnod[plotnum].plotn = plotn; 
}

/*------------------------------------*/

void plotsize(double plotsiz,  int nplot)

{
   if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
   plotnod[nplot].plotsiz = plotsiz; 
}

/*------------------------------------*/

void plotsize(double plotsiz)

{
   if (plotnum>PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
   plotnod[plotnum].plotsiz = plotsiz; 
}

/*------------------------------------*/

void plotmode(int pmod,  int nplot)

{
   if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
   plotnod[nplot].pmod = pmod; 
}

/*------------------------------------*/

void plotmode(int pmod)

{
   if (plotnum>PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
   plotnod[plotnum].pmod = pmod; 
}

/*------------------------------------*/

void plotval(double plotvl)

{
   if (plotnum>PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
   plotnod[plotnum].plotval = plotvl; 
   plotnod[plotnum].plotval2 = LARGENUM; 
   plotnod[plotnum].plotval3 = LARGENUM; 
   plotnod[plotnum].plotval4 = LARGENUM; 
   plotnod[plotnum].plotval5 = LARGENUM; 
   plotnod[plotnum].plotval6 = LARGENUM; 
   plotnod[plotnum].plotval7 = LARGENUM; 
}

/*------------------------------------*/

void plotval(double plotvl,  int nplot)

{
   if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
   plotnod[nplot].plotval = plotvl; 
   plotnod[nplot].plotval2 = LARGENUM; 
   plotnod[nplot].plotval3 = LARGENUM; 
   plotnod[nplot].plotval4 = LARGENUM; 
   plotnod[nplot].plotval5 = LARGENUM; 
   plotnod[nplot].plotval6 = LARGENUM; 
   plotnod[nplot].plotval7 = LARGENUM; 
}

/*------------------------------------*/

void plotval(double plotvl, double plotvl2, int nplot)

{
   if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
   plotnod[nplot].plotval  = plotvl; 
   plotnod[nplot].plotval2 = plotvl2; 
   plotnod[nplot].plotval3 = LARGENUM; 
   plotnod[nplot].plotval4 = LARGENUM; 
   plotnod[nplot].plotval5 = LARGENUM; 
   plotnod[nplot].plotval6 = LARGENUM; 
   plotnod[nplot].plotval7 = LARGENUM; 
}

/*------------------------------------*/

void plotval(double plotvl, double plotvl2, double plotvl3, int nplot)

{
   if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
   plotnod[nplot].plotval  = plotvl; 
   plotnod[nplot].plotval2 = plotvl2; 
   plotnod[nplot].plotval3 = plotvl3; 
   plotnod[nplot].plotval4 = LARGENUM; 
   plotnod[nplot].plotval5 = LARGENUM; 
   plotnod[nplot].plotval6 = LARGENUM; 
   plotnod[nplot].plotval7 = LARGENUM; 
}

/*------------------------------------*/

void plotval(double plotvl, double plotvl2, double plotvl3, double plotvl4, int nplot)

{
   if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
   plotnod[nplot].plotval  = plotvl; 
   plotnod[nplot].plotval2 = plotvl2; 
   plotnod[nplot].plotval3 = plotvl3; 
   plotnod[nplot].plotval4 = plotvl4; 
   plotnod[nplot].plotval5 = LARGENUM; 
   plotnod[nplot].plotval6 = LARGENUM; 
   plotnod[nplot].plotval7 = LARGENUM; 
}

/*------------------------------------*/

void plotval(double plotvl, double plotvl2, double plotvl3, double plotvl4, double plotvl5, int nplot)

{
   if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
   plotnod[nplot].plotval  = plotvl; 
   plotnod[nplot].plotval2 = plotvl2; 
   plotnod[nplot].plotval3 = plotvl3; 
   plotnod[nplot].plotval4 = plotvl4; 
   plotnod[nplot].plotval5 = plotvl5; 
   plotnod[nplot].plotval6 = LARGENUM; 
   plotnod[nplot].plotval7 = LARGENUM; 
}

/*------------------------------------*/

void plotval(double plotvl, double plotvl2, double plotvl3, double plotvl4, double plotvl5, 
			    double plotvl6, int nplot)

{
   if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
   plotnod[nplot].plotval  = plotvl; 
   plotnod[nplot].plotval2 = plotvl2; 
   plotnod[nplot].plotval3 = plotvl3; 
   plotnod[nplot].plotval4 = plotvl4; 
   plotnod[nplot].plotval5 = plotvl5; 
   plotnod[nplot].plotval6 = plotvl6; 
   plotnod[nplot].plotval7 = LARGENUM; 
}

/*------------------------------------*/

void plotval(double plotvl, double plotvl2, double plotvl3, double plotvl4, double plotvl5, 
			    double plotvl6, double plotvl7, int nplot)

{
   if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
   plotnod[nplot].plotval  = plotvl; 
   plotnod[nplot].plotval2 = plotvl2; 
   plotnod[nplot].plotval3 = plotvl3; 
   plotnod[nplot].plotval4 = plotvl4; 
   plotnod[nplot].plotval5 = plotvl5; 
   plotnod[nplot].plotval6 = plotvl6; 
   plotnod[nplot].plotval7 = plotvl7; 
}

/*------------------------------------*/

void plotrst(int totplots)
{
   int i;

 for (i=0; i<totplots; i++)                     /* restart all graphs */
   plotnod[i].oldy = plotnod[i].oldx = LARGENUM;

 if (! vidmode) {
    ncfprintf (stdout,"#restart plot...\n");
 }

}

/*------------------------------------*/

void plotinit (int totplots)
{
    int i,n, pmod2;
    int nnodpos,charfl,charmode,plotn,ppen,useframe; 
    const char *modech,*plframe;
    char fstr[40];
    static char modey[5];

#ifdef DEBUG 
 if (debug & NCPLOT) ncfprintf (stderr,"plotinit %d\n",totplots);
#endif

pinitfl = 1;
yfieldwidth = yprecision + 5;

if (dttplots < totplots) dttplots = totplots;

if (!notinit(ssetxmax))  setxmax = ssetxmax;  /* command-line override */
if (!notinit(setxmax))  plotxmax = setxmax;  /* control var override */
else		        plotxmax = endexp;

if (!notinit(ssetxmin)) setxmin  = ssetxmin;
if (!notinit(setxmin))  plotxmin =  setxmin;
else		        plotxmin =  simtime;

if (!notinit(setymax)) plmax = setymax;
if (!notinit(setymin)) plmin = setymin;


if (!vidmode) {                                /* init non-graphics display */
   if (complin[0]=='#')
	ncfprintf (stdout,"%.200s",complin);	/* print "#c xx comps..." */
   if (stfile)
        ncfprintf (stdout,"#f stim file %s\n",stfile);
   if (infile && srcfile) { 
        if (strcmp (infile,srcfile)!=0)
        ncfprintf (stdout,"#f input file %s\n",infile);
    }
    else if (infile)
        ncfprintf (stdout,"#f input file %s\n",infile);
    if (srcfile)
        ncfprintf (stdout,"#f source file %s\n",srcfile);
    if (plotlabel)
        ncfprintf (stdout,"#f plotlabel %s\n",plotlabel);
    ncfprintf (stdout,"#p %d plots plsep %d plmax %g plmin %g\n",
					totplots,plsep,plmax,plmin);
    ncfprintf (stdout,"#e begin %g  endexp %g  rseed  %d version %s\n",
					plotxmin,plotxmax,rseed,
					print_version(ncversion));
    for (i=0; i<totplots; i++) {
      n = looklist(i);      	/* command line overrides order */


      if (plotnod[n].plname[0])
          ncfprintf (stdout,"#n \"%s\"\n", plotnod[n].plname);

      switch(plotnod[n].pmod) {
	 case NRECGV:
         case VREC:     modech = "V"; break;
         case WREC:     modech = "V"; break;
	 case NRECGI:
         case IREC:     modech = "I"; break;
         case MREC:     modech = "M"; break;
         case LREC:     modech = "L"; break;
         case GREC:     modech = "X"; break;
         case SREC:     modech = "S"; break;
         case FREC:     
			if ((pmod2 = plotnod[n].pmod2) > 0) {
			   switch (pmod2) {
         		      case IREC: case NRECGI: sprintf (modey,"UI"); break; 
         		      case NRECG0: case NRECG1: case NRECG2: case NRECG3:
			         sprintf (modey,"UG%d",pmod2 - NRECG0); break;
         		      case NRECG: sprintf (modey,"UG%d",plotnod[n].pval); break;
			      default: modey[0] = 0;
			    }
			    modech = modey;
			} else { modech = "U";}
			break;
	 case CREC_GLU:  case CREC_AMPA: case CREC_KAINATE: case CREC_NMDA: case CREC_CNQX:
	 case CREC_GABA: case CREC_BIC:  case CREC_PTX:  case CREC_GLY: case CREC_CHRC:
	 case CREC_STRY:  case CREC_CAMP: case CREC_CGMP: case CREC_PH: case CREC_ATP:
			 sprintf (modey,"N%d",plotnod[n].pmod - CREC_GLU);
			 modech = modey;
			 break;
	 case CREC_CABUF:
	 case CREC_CABUFB:
	 case CREC_CACONC:
	 case CREC_CAS:
			 sprintf (modey,"Y%d",plotnod[n].pval);
			 modech = modey;
			 break;
	 case CREC_CAV: case CREC_CAIT: case CREC_CAIP: 
	 case CREC_CAIE: case CREC_CAIPE: case CREC_CICR:
			 sprintf (modey,"C%d",plotnod[n].pmod - CREC_CAV);
			 modech = modey;
			 break;
	 case NRECA0: case NRECA1: case NRECA2: case NRECA3: case NRECA4:
	 case NRECH0: case NRECH1: case NRECH2: case NRECH3: case NRECH4:
	 case NRECB0: case NRECB1: case NRECB2: case NRECB3: case NRECB4:
	 case NRECC0: case NRECC1: case NRECC2: case NRECC3: case NRECC4:
	 case NRECC9:	modech = "F"; break; 
         case NRECA8:
         case NRECA9:   modech = "R"; break;
         case NRECG0: case NRECG1: case NRECG2: case NRECG3:
         		sprintf (modey,"G%d",plotnod[n].pmod - NRECG0);
			modech = modey;
		      break;

         case NRECG:    sprintf (modey,"G%d",plotnod[n].pval);
			modech = modey;
		      break;
 
	 case NRECGC: case NRECGM: case NRECGH:
	   		modech = "F"; break;
         case CABLE:    modech = "V"; break;
         default:       modech = "0"; break;
      }

    nnodpos = 1;
    if (plotnod[n].cnod2 > NULLVAL) nnodpos++;
    if (plotnod[n].cnod3 > NULLVAL) nnodpos++;
    if (plotnod[n].cnod4 > NULLVAL) nnodpos++;

  if (notinit(plotnod[n].pxmax))            /* default values for xmax,xmin? */
  switch (nnodpos) {

  case 4:
   ncfprintf(stdout,"#x node %-4d %-2d %-2d %-2d mode %s\n",
                plotnod[n].cnod1, plotnod[n].cnod2, 
                plotnod[n].cnod3, plotnod[n].cnod4, 
		modech);
   break;
  case 3:
   ncfprintf(stdout,"#x node %-4d %-2d %-2d %-2d mode %s\n",
                plotnod[n].cnod1, plotnod[n].cnod2,
		plotnod[n].cnod3, NULND, modech);
   break;
  case 2:
   ncfprintf(stdout,"#x node %-4d %-2d %-2d %-2d mode %s\n",
                plotnod[n].cnod1, plotnod[n].cnod2, NULND, NULND, modech);
   break;
  case 1:
   ncfprintf(stdout,"#x node %-4d %-2d %-2d %-2d mode %s\n",
                plotnod[n].cnod1, NULND,  NULND,  NULND, modech);
   break;

    }  /* switch (nnodpos) */

  else
  switch (nnodpos) {		/* not default values for xmax, xmin */
   case 4:
    ncfprintf(stdout,
	"#x node %-4d %-2d %-2d %-2d mode %s xmax %-8.3g xmin %-8.3g\n",
                plotnod[n].cnod1, plotnod[n].cnod2, 
                plotnod[n].cnod3, plotnod[n].cnod4, 
		modech, 
                plotnod[n].pxmax, plotnod[n].pxmin);
   break;
   case 3:
    ncfprintf(stdout,
	"#x node %-4d %-2d %-2d %-2d mode %s xmax %-8.3g xmin %-8.3g\n",
                plotnod[n].cnod1, plotnod[n].cnod2, 
		plotnod[n].cnod3, NULND, modech, 
                plotnod[n].pxmax, plotnod[n].pxmin);
   break;
  case 2:
   ncfprintf(stdout,
	"#x node %-4d %-2d %-2d %-2d mode %s xmax %-8.3g xmin %-8.3g\n",
                plotnod[n].cnod1, plotnod[n].cnod2, NULND, NULND, modech, 
                plotnod[n].pxmax, plotnod[n].pxmin);
   break;
  case 1:
   ncfprintf(stdout,
	"#x node %-4d %-2d %-2d %-2d mode %s xmax %-8.3g xmin %-8.3g\n",
                plotnod[n].cnod1, NULND, NULND, NULND, modech, 
                plotnod[n].pxmax, plotnod[n].pxmin);
   break;
}

    charfl = plotnod[n].charfl;
    if (charfl==0) charfl = ' '; 
    charmode = plotnod[n].charmode;
    plframe = plotnod[n].plframe; 
    if ((plotn=plotnod[n].plotn)==NULLVAL) plotn = -1;
    if ((ppen=plotnod[n].ppen)==NULLVAL)   ppen = -1;
    if (*plframe) sprintf (fstr,"fr %s",plframe);
     else strcpy (fstr, "       ");
    if (charmode==PLINES) {
     if (notinit(plotnod[n].pymax))          /* default values for ymax,ymin? */
       ncfprintf(stdout,"#y %s   pl %d siz %g pen %d styl %c\n", 
		fstr,plotn,plotnod[n].plotsiz,
		ppen,charmode);
     else
       ncfprintf(stdout,
	"#y %s   pl %d siz %g pen %d styl %c ymax %-8.3g ymin %-8.3g\n", 
		fstr,plotn,plotnod[n].plotsiz,
		ppen,charmode, plotnod[n].pymax, plotnod[n].pymin);
    } /* if (charmode==) */
    else {
       double csize;
     csize = plotnod[n].csize;
     if (notinit(plotnod[n].pymax))          /* default values for ymax,ymin? */
       ncfprintf(stdout, 
	"#y %s   pl %d siz %g pen %d styl %c char '%c' siz %g\n", 
		fstr,plotn,plotnod[n].plotsiz,
		ppen,charmode,charfl,csize);
     else
       ncfprintf (stdout, 
	"#y %s  pl %d siz %g pen %d styl %c char '%c' siz %g ymax %-8.3g ymin %-8.3g\n",
                fstr,plotn,plotnod[n].plotsiz,
		ppen,charmode,charfl,csize,
		plotnod[n].pymax, plotnod[n].pymin);
    } /* if (charmode==) */

   }    /* for (i=0; i<totplots; )  */

      ncfprintf (stdout,"# \n");
      ncfprintf (stdout,"%-8.8s ","#node--> ");
      for (i=0; i<totplots; i++) {
#define NUMBUFSIZ 50
	     char buf[NUMBUFSIZ];
          n = looklist(i);      	/* command line overrides order */
#undef NUMBUFSIZ

  switch (plotnod[n].pmod) {
	case VREC: case IREC: case MREC: case SREC:  case FREC: case WREC:
	case CREC_GLU:  case CREC_AMPA: case CREC_KAINATE: case CREC_NMDA: case CREC_CNQX:
	case CREC_GABA: case CREC_BIC:  case CREC_PTX:  case CREC_GLY: case CREC_CHRC:
	case CREC_STRY:  case CREC_CAMP: case CREC_CGMP: case CREC_PH: case CREC_ATP:

	case CREC_CABUF:
	case CREC_CABUFB:
	case CREC_CACONC:
	case CREC_CAS: 
	case CREC_CAV: case CREC_CAIT: case CREC_CAIP: 
	case CREC_CAIE: case CREC_CAIPE: case CREC_CICR:

    	case  LREC: case  CABLE:
  		if (plotnod[n].cnod3 > NULLVAL) 
		  sprintf (buf,"[%-d %-d %-d]",
			plotnod[n].cnod1,plotnod[n].cnod2,plotnod[n].cnod3);
  		else if (plotnod[n].cnod2 > NULLVAL) 
		  sprintf (buf,"[%-d %-d]",plotnod[n].cnod1,plotnod[n].cnod2);
  		else
		  sprintf (buf,"[%-d]",plotnod[n].cnod1);
		break;
	case  GREC:
		  sprintf (buf,"(%-d)",n);
		break;
    	case NRECA0: case NRECA1: case NRECA2: case NRECA3: case NRECA4:
    	case NRECA8: case NRECA9:
	case NRECH0: case NRECH1: case NRECH2: case NRECH3: case NRECH4:
    	case NRECB0: case NRECB1: case NRECB2: case NRECB3: case NRECB4:
    	case NRECC0: case NRECC1: case NRECC2: case NRECC3: case NRECC4:
    	case NRECC9:

	case NRECG0: case NRECG1: case NRECG2: case NRECG3:
	case NRECG: case NRECGV: case NRECGI: case NRECGC: case NRECGM: 
	case NRECGH:
  		sprintf (buf,"(%-d)",plotnod[n].cnod1);
		break;
	}   /* switch */
        ncfprintf (stdout,"%-*2$.*3$s  ",buf,yfieldwidth,yprecision+6);
       }

        ncfprintf (stdout,"\n");
        fflush (stdout);
        return;			/* return if not vidmode */

   }  /* if (!vidmode) */

  if (finfl) { nonodes=1; nofilfl=1; }

  for (useframe=i=0; i<totplots; i++) {     /* check for frames used */ 
    n = looklist(i);      		    /* command line overrides order */
    if (plotnod[n].plframe[0]) useframe=1; 
  }

  /* The idea here is that "drawframe" draws one X-axis and maybe
     multiple Y-axes.  If one wants to use frames, it might make
     sense to have a new X-axis in each frame. However this has
     not been fully developed, so for now we allow frames without
     making a new X-axis for each */

  //useframe = 0;		/* only call "drawframe()" once -- only 1 X-axis */

  if (useframe) {
      char *frames[PLOTNODSIZ];
      int k;

    for (i=0; i<PLOTNODSIZ; i++) frames[i]=(char *)NULL;

    for (i=0; i<totplots; i++) {            /* check for first frame */
      n = looklist(i);      		    /* command line overrides order */
      if (plotnod[n].plframe[0]) {	    /* diff than root frame?  */
         for (k=0; k<PLOTNODSIZ; k++) {	    /* check against existing frames */ 
           if (frames[k]) {
              if (strcmp (frames[k],plotnod[n].plframe)==0) {
                 plotnod[n].newframe = 0;   /* frame already started */
                 break;
              }
           }
           else {
            frames[k] = plotnod[n].plframe;
            plotnod[n].newframe = 1;
            break;
           }  
        }
      }   /* if plotnod[n].plframe[0] */
      if (plotnod[n].newframe)
         drawframe (n,totplots,useframe);	/* draw new frame with X-axis */
         //drawframe (n,totplots,0);
    } 
  }
  else drawframe (0,totplots,useframe);

#ifdef DEBUG 
  if (debug & NCPLOT) ncfprintf (stderr,"plotinitend\n");
#endif

}

/*------------------------------------*/

int getlengthi (int n)

/* get length of a number that will be displayed */

{
  char tbuf[20];

  sprintf (tbuf,"%g ",(double)n);
  return strlen (tbuf); 
}

/*------------------------------------*/

int getlengthd (double n)

/* get length of a number that will be displayed */

{
  char tbuf[20];

  sprintf (tbuf,"%g", n);
  return strlen (tbuf); 
}

/*------------------------------------*/

int getaplot (int i)

/* Find logical plot. */
/* "aplots" must already contain list of logical plots */

{
     int aplot,n,pln;

 n = looklist(i);      	/* command line overrides order */

 if (plsep) {
     if ((pln=plotnod[n].plotn) < 0) pln = i;
     aplot = aplots[pln];
  }
  else aplot = i;			/* logical plot */
 return aplot;
}

/*------------------------------------*/

void draw_plotlabel(int dispsize)
{
     int i,sm,u,lg,labelen;
     double lablx, lably;
     char *pllabel=NULL;

#define islgupper(c)      ((c) == 'M' || (c) == 'W' || (c) == '_')
#define isupper(c)        ('A' <= (c) && (c) <= 'Z')
#define islgdigit(c)      ((c)=='0' || '2' <= (c) && (c) <= '9')
#define issmall(c)        ((c)=='-' || (c)=='-' || (c)=='i' || (c)=='I')

   // if (disp_ray) return;

   lablx = 0.97;
   lably = 0.98;
   if (plotlabel) {
	pllabel = plotlabel;
   }
   else if (infile && srcfile) {		/* if both infile and srcfil defined */
	pllabel = infile;			/*  print infile only. */
   }
   else if (srcfile) {
	pllabel = srcfile;
   }
   else if (infile) { 
	pllabel = infile;
   }
   // if (pllabel==NULL) pllabel = "stdin";
   if (pllabel!=NULL) {
     labelen = strlen(pllabel);
     for (i=sm=u=lg=0; i<labelen; i++) {
          if (isupper(pllabel[i]) || islgdigit(pllabel[i])) u++;
          if (islgupper(pllabel[i])) lg++;
          if (issmall(pllabel[i])) sm++;
     }
     // lablx -= charsiz * ((labelen-u)*0.72 + u*1.1 + lg*0.3) ;
     lablx -= charsiz * ((labelen)*0.72 - sm*0.2 + u*0.1 + lg*0.3) ;	// compromise between vid and ps chars
     if (disp_ray) {
       rprintf(pllabel,lablx*dispsize*0.5,lably*dispsize*0.5,1,charsiz*dispsize,BLACK);
     } else {
       gcwidth (charsiz);
       gpen (WHITE);
       gmove (lablx, lably);
       prtext ("l",charsiz,pllabel);
     }
   }
}

#undef islgupper
#undef isupper
#undef isdigit

/*------------------------------------*/

void drawframe (int nplot, int totplots, int useframe)

/* Draw one frame on the screen using parameters
    stored in "plotnod[]".
*/

{
    double axmin,axmax,range,naxmin,naxmax,nrange,srange,yaxisX;
    double yaxis,ybase,ticmin,ticinc;
    double plaxis, plticmin, pltic, calibtic;
    double ppmin, ppmax, pxnum, pynum, floormul;
    double ticypos,ticymin,ticymax,ycharsiz,totsiz;
    double ylablwidth, ylablwidtht, ylablwidthmax, ylwidth, ylwidthmax;
    int ticyminfl,ticymaxfl,lablint,nolabltic; 
    int i,j,k,m,n,ttotplots,itic,itici,ntics,hitics,calibmv,calibms;
    int nnodpos,plotype,pen,charfl,charmode,ticdir;
    int numlaxisx, setaxis, rightYax, axpen, aplot;
    int yticchar, yticcharmax;
    int plmaxlen, plminlen, pllenmax;
    const char *calstr, *xcalstr; 
    const char *plframe;
    char fstr[40];

#ifdef DEBUG 
  if (debug & NCPLOT) ncfprintf (stderr,"drawframe\n");
#endif


				/* fix max, min for x values */
 for (i=0; i<totplots; i++) {
   n = looklist(i);      	/* command line overrides order */
   plotype = plotnod[n].pmod;
   if (notinit(plotnod[n].pxmax))                /* find default values for */
        switch (plotype) {
	  default:
                plotnod[n].pxmax  = plotxmax;     /*  x axis max and min */
                plotnod[n].pxmin  = plotxmin;
                break;
        case GREC:
        case FREC:
        case SREC:
		if (n==0) {
                  plotnod[n].pxmax  = plotxmax; /*  x axis max and min */
                  plotnod[n].pxmin  = plotxmin;
		}
		else {
                  plotnod[n].pxmax  = plotnod[0].pxmax;
                  plotnod[n].pxmin  = plotnod[0].pxmin;
   		  if (notinit(plotnod[n].pxmax)) {    /* find default values for */
                    plotnod[n].pxmax  = plotxmax;     /*  x axis max and min */
                    plotnod[n].pxmin  = plotxmin;
		  }
		}
                break;
        }  /* switch */

   if (!notinit(setxmax)) plotnod[n].pxmax = setxmax; 
   if (!notinit(setxmin)) plotnod[n].pxmin = setxmin; 

   plotnod[n].xrange = plotnod[n].pxmax - plotnod[n].pxmin;
 }

				/* fix max, min for y values */

 for (i=0; i<totplots; i++) {
   n = looklist(i);      		/* command line overrides order */
   if (notinit(plotnod[n].pymax))        /* find default values for */
        switch (plotnod[n].pmod) {
	case NRECGV:
	case NRECGI:
	case CREC_CAV: 
        case VREC:
        case WREC:
        case MREC:
        case IREC:
                plotnod[n].pymax  = plmax;      /*  y axis max and min */
                plotnod[n].pymin  = plmin;
                break;
        case LREC:
                plotnod[n].pymax  = 1000; 	/*  y axis max and min */
                plotnod[n].pymin  = 0;
                break;
        case GREC:
        case FREC:
        case SREC:
		if (n==0) {
                  plotnod[n].pymax  = plmax; /* y axis max and min */
                  plotnod[n].pymin  = plmin;
		}
		else {
                  plotnod[n].pymax  = plotnod[0].pymax; 
                  plotnod[n].pymin  = plotnod[0].pymin;
   		  if (notinit(plotnod[n].pymax)) {        /* find default values for */
                    plotnod[n].pymax  = plmax; 		  /* y axis max and min */
                    plotnod[n].pymin  = plmin;
		  }
		}
                break;
	case NRECA0: case NRECA1: case NRECA2: case NRECA3:
	case NRECA4: 
                plotnod[n].pymax  = 0.05; 	/*  y axis max and min */
                plotnod[n].pymin  = 0.0;
                break;
	case NRECH0: case NRECH1: case NRECH2: case NRECH3: case NRECH4:
                plotnod[n].pymax  =  0.05; 	/*  y axis max and min */
                plotnod[n].pymin  = -0.05;
                break;
	case NRECA8:
	case NRECA9:
                plotnod[n].pymax  = 500; 	/*  y axis max and min */
                plotnod[n].pymin  = 0;
                break;

	case CREC_GLU:  case CREC_AMPA: case CREC_KAINATE: case CREC_NMDA: case CREC_CNQX:
	case CREC_GABA: case CREC_BIC:  case CREC_PTX:  case CREC_GLY: case CREC_CHRC:
	case CREC_STRY:  case CREC_CAMP: case CREC_CGMP:  case CREC_ATP:
                plotnod[n].pymax  = 100e-6;	/*  y axis max and min */
                plotnod[n].pymin  = 0;
                break;

	case CREC_PH:
                plotnod[n].pymax  = 14;	/*  y axis max and min */
                plotnod[n].pymin  = 0;
                break;

	case CREC_CABUF:
	case CREC_CABUFB:
	case CREC_CACONC:
	case CREC_CAS: 
		if (plotnod[n].pval==0) {
                 plotnod[n].pymax  = 1e-3;	/*  y axis max and min */
                 plotnod[n].pymin  = 0;
		}
		else {
                 plotnod[n].pymax  = 100e-6;	/*  y axis max and min */
                 plotnod[n].pymin  = 0;
		}
                break;

   case CREC_CAIT: case CREC_CAIP: case CREC_CAIE: case CREC_CAIPE: case CREC_CICR:

                plotnod[n].pymax  = 1e-9;	/*  y axis max and min */
                plotnod[n].pymin  = 0;
                break;


	case NRECB0: case NRECB1: case NRECB2: case NRECB3: case NRECB4:
	case NRECC0: case NRECC1: case NRECC2: case NRECC3: case NRECC4:
	case NRECC9:

	case NRECG1: case NRECG2: case NRECG3:
	case NRECGC: case NRECGM: case NRECGH:
                plotnod[n].pymax  = 2; 		/*  y axis max and min */
                plotnod[n].pymin  = 0;
                break;
	case NRECG: if (plotnod[n].pval>0) {
                     plotnod[n].pymax  = 2; 	/*  y axis max and min */
                     plotnod[n].pymin  = 0;
		    }
			/* no break, if =0, go to next case */	
	case NRECG0: 				/* conductance of chan */
                plotnod[n].pymax  = 1e-8;	/*  y axis max and min */
                plotnod[n].pymin  = 0;
		break;
	case CABLE:
                plotnod[n].pymax  = plmax;      /*  y axis max and min */
                plotnod[n].pymin  = plmin;
               	break; 
        }

	/* set actual plot size from overrides */

   //if (i==0) {		/* if first plot in list */
      if (!notinit(setymax)) plotnod[n].pymax = setymax; 
      if (!notinit(setymin)) plotnod[n].pymin = setymin; 
   //}
   plotnod[n].yrange = plotnod[n].pymax - plotnod[n].pymin;
 }

 for (i=0; i<totplots; i++) {                     /* restart all graphs */
   plotnod[i].oldy = plotnod[i].oldx = LARGENUM;    /* use plot order because */
 }						/* plots may be duplicated */

 /* find number and location of logical plots, set by "plotn = " */
 
 for (i=0; i<PLOTNODSIZ; i++) {
   aplots[i] = 0;			/* histogram to find plots used */
   splots[i] = 0;			/* size of first trace in plot */
 }
 for (i=0; i<totplots; i++) {		/* make histogram of all plots */
        double plotsiz;

   n = looklist(i);
   if ((k=plotnod[n].plotn) < 0) k = i;
   if (k<0) k = 0;
   if (k>=PLOTNODSIZ) k = PLOTNODSIZ-1;
   if (aplots[k]==0) {
	plotsiz = plotnod[n].plotsiz;
        if (plotsiz==0) plotsiz = 1;
	splots[n] = plotsiz;			/* set first one's size */
   }
   aplots[k]++;
 }
 for (naplots=i=0; i<PLOTNODSIZ; i++) {
   if (aplots[i] > 0) {
	aplots[i] = ++naplots;	/* set logical plot, find num of plots */ 
   }				/*  value is 1 too much ( > 0 ) */
 }
 for (totsiz=i=0; i<PLOTNODSIZ; i++) {	/* normalize total plot size to 1 */
    totsiz += splots[i];
 }
 if (totsiz < 1) totsiz = 1;		/* allow total size to be small */
 for (i=0; i<PLOTNODSIZ; i++) {		/* normalize total plot size */
    if (totsiz>0) splots[i] *= naplots/totsiz;
 }
 for (totsiz=k=i=0; i<PLOTNODSIZ; i++) {  /* For all logical plots, */
	int first, plotn;		  /*  make cumulative yaxis size */
        double plotsiz;

    if (aplots[i]>0) {
        for (first=1,k=0; k<PLOTNODSIZ; k++) {  /* For each plot assigned */
           n = looklist(k);      		/* command line overrides order */

	  if ((plotn=plotnod[n].plotn)==NULLVAL) plotn = n;
	  if (plotn==i) {	      		/*  to this logical plot */
	    scplots[n] = totsiz;		/* save its ybase */
	    if (first) {
		first = 0;
		plotsiz = splots[n];
	    }
	    else {
	       splots[n] = plotsiz;	/* override size of other traces */
	    }				/*  so all are the same */

	    /* ncfprintf (stderr,"i %-2d n %-2d t %g sc %g\n",
						i,n,totsiz,scplots[n]); /* */
	  }
	}
	totsiz += plotsiz;
        aplots[i]--;			/* restore aplot value */
    }
 }

 
 setaxis = 0;
 if (useframe) {
   ttotplots = 1;
   numyaxis = 1;
 }
 else {
   if (totplots == 0) totplots = 1;
   if (plsep) {
      if (naplots)
	ttotplots = naplots;
      else
	ttotplots = totplots;
   }
   else       ttotplots = 1;
   numyaxis = numlaxis+numraxis;
   setaxis = (numyaxis>0); 
   if (numyaxis<=0) {  			/* if no y axes set on command line */
       leftaxes[0]=1;
       numlaxis = 1;
       numyaxis = ttotplots;
   }
   if (setaxis) plsep = 0; 
 }
  
 /* setmode (fileno(stdout), O_BINARY); 		/* for microsoft C */
 /* frame (rootframe);  */

 cwidth (charsiz,".");

 if (ncerror) {
    cpen (WHITE);
    move (0.9, 0.9);
    cwidth (.025,".");
    prtext ("l",0.025,"Error\n");
    cwidth (charsiz,".");
 }
 cpen (BLUE);

 yplotmarg = YPLOTMARG + 2*XPLOTLABEL+charsiz*2; /* calculate y axis length */
 if (plsep) {
     ycharsiz = ((ttotplots > 4) ? charsiz * 3.0 / pow(ttotplots,0.7): charsiz);
     ybase = yplotmarg;
  }
  else {		/* make y axis smaller if x axis will be smaller */
     ycharsiz = charsiz;
     ybase = yplotmarg;
  }

		/* find max width of left margin from labels */

 for (i=0; i<PLOTNODSIZ; i++) {
   plotlabl[i] = 0;
   nnplots[i] = 0;
   yplcharsiz[i] = charsiz * 1.4;
 }
 for (i=0; i<totplots; i++) {		/* count labeled traces in plot */
     int n, aplot;
    
   n = looklist(i);

   aplot = getaplot(i);

   if (plotnod[n].plname[0]) {
      plotlabl[aplot] = 1; 
      nnplots[aplot]++;			/* count named traces in plot */
   }
   else if (!plotlabl[aplot]) {
      nnplots[aplot]++;			/* count numbered traces in plot */
   }
 }

 yticcharmax = 0;
 ylablwidthmax = ylwidthmax = 0;
 for (i=0; i<totplots; i++) {		/* check left margin size */
     int l1, l2;
     int yticchar;
     double t1,t2;
     double yclimit;

   n = looklist(i);
   aplot = getaplot(i);

   switch (plotnod[n].pmod) {		/* Find calib multiplier. */
     case NRECGV: 			/*  This is also done below, */
     case CREC_CAV: 			/*  should optimize. */
     case  WREC:
     case  VREC:
     case  CABLE:
          if (ncabs(plotnod[n].pymax) < 1) {      /* y axis calibration units */
             calibmv = 1; 
          }
	  break;
     default: calibmv = 0;
	     break;
   }
   t1 = plotnod[n].pymax;
   t2 = plotnod[n].pymin;
   if (calibmv) {
     t1 *= 1000.0;
     t2 *= 1000.0;
     t1 = int(t1);
     t2 = int(t2);
   }
   l1 = getlengthd(t1);  /* size of tic labels */ 
   l2 = getlengthd(t2); 
   yticchar = max(l1,l2);
   yticchar = min(yticchar,6);

   if (yticcharmax < yticchar) yticcharmax = yticchar;

   nnodpos = 0;
   if (plotnod[n].plname[0])    nnodpos = strlen(plotnod[n].plname);
   else {
     if (plotnod[n].cnod1 > NULLVAL) nnodpos+=getlengthi (plotnod[n].cnod1);
     if (plotnod[n].cnod2 > NULLVAL) nnodpos+=getlengthi (plotnod[n].cnod2);
     if (plotnod[n].cnod3 > NULLVAL) nnodpos+=getlengthi (plotnod[n].cnod3);
   }
   if (nonodes) ylablwidth = 3;
   else         ylablwidth = nnodpos+1; 

   /* ncfprintf (stderr,"%s ylablwidth %g yticchar %d ycharsiz %g \n",
   		plotnod[n].plname,ylablwidth,yticchar,ycharsiz); /* */

		/* Set char size for plot labels, for use below. */
		/* Problem is that when plot is assigned smaller size */
		/*  char size must also be set smaller to fit. */

   yaxis = splots[n] / ttotplots * (1-yplotmarg*1.5) -
					YAXISSPC*10.0/(ttotplots+10.0);

   yclimit = yaxis * 0.8 / nnplots[aplot];   /* set maximum char size */
   if (yplcharsiz[aplot] > yclimit) yplcharsiz[aplot] = yclimit;
 
   ylablwidtht = ylablwidth * yplcharsiz[aplot];
 
   if (ylablwidtht > 0.2) {
	   double t;
        /* ncfprintf (stderr,"1 ylablwidtht %g aplot %d yplcharsiz %g\n", ylablwidtht,aplot,yplcharsiz[aplot]); /* */
	t = 0.2 / ylablwidth * pow(0.028/yplcharsiz[aplot],1.4);
	if (yplcharsiz[aplot] > t) yplcharsiz[aplot] = t;
        /* ncfprintf (stderr,"1 ylablwidtht %g aplot %d yplcharsiz %g\n", ylablwidtht,aplot,yplcharsiz[aplot]); /* */
   }
   
   ylablwidtht = ylablwidth * yplcharsiz[aplot];
   if (ylablwidthmax < ylablwidtht) ylablwidthmax = ylablwidtht; /* total label width when not plsep */

   ylwidth = 0.9*ylablwidth*yplcharsiz[aplot] + (yticchar-0.5)*ycharsiz; /* total */
   if (ylwidthmax < ylwidth) ylwidthmax = ylwidth; 		      /* total label width for plsep */
   /* ncfprintf (stderr,"ylwidth %g ylablwidth %g yplcharsiz %g pow %g yticchar %d ycharsiz %g\n", 
		   ylwidth,ylablwidth,yplcharsiz[aplot],pow(yplcharsiz[aplot]/0.029,1),yticchar,ycharsiz); /* */
 }
    /* ncfprintf (stderr,"ylwidthmax %g\n", ylwidthmax); /* */

 n = looklist(nplot);

 numlaxisx = numlaxis - 1;
 if (numlaxisx<0) numlaxisx = 0;

 if (plsep) {
   xplotmargl = ylwidthmax; 
   xplotmargr = XPLOTMARGR;
 }
 else {
   xplotmargl = ylablwidthmax+(numlaxis>0)*(yticcharmax+.1)*ycharsiz + (numlaxisx) * YAXISWIDTH;
   xplotmargr = XPLOTMARGR + (numraxis) * YAXISWIDTH;
 }

 plframe = plotnod[looklist(nplot)].plframe;
 if (*plframe) frame (plframe);
 cwidth (charsiz,".");

 move (xplotmargl,  ybase);                 /* draw x axis */
 draw (1-xplotmargr,ybase);

 plaxis = 1 - (xplotmargl + xplotmargr);       	/* calc size of plot x axis */

 axmax = plotnod[n].pxmax;
 axmin = plotnod[n].pxmin;
 range = plotnod[n].xrange;
 if (range == 0) range = 1e-20;
 ticinc = geticinc(range);                      /* find x axis tic incr */

 if (axmin >= 0) floormul = 1.000001;
 else            floormul = 0.999999;
 itic   = (int)floor((axmin * floormul) / ticinc);
 ticmin = itic * ticinc;  
 naxmin = ticmin;
 if (axmax < 0) floormul = 1.000001;
 else           floormul = 0.999999;
 itic   = (int)ceil((axmax * floormul) / ticinc); 
 naxmax = itic * ticinc;  
 nrange = (naxmax - naxmin);

/* plotnod[n].xrange = nrange;			/* rescale x axis */
/* plotnod[n].pxmax = naxmax;
 plotnod[n].pxmin = naxmin;
*/

 if (!useframe) {
  for (i=0; i<totplots; i++) {			/* rescale other x axes */
     int nn;
     double trange;

   nn = looklist(i); 
   trange = plotnod[nn].xrange;
   plotnod[nn].xrange *= nrange/range;
   plotnod[nn].pxmax += (naxmax-axmax)/range*trange;
   plotnod[nn].pxmin += (naxmin-axmin)/range*trange;
  }
 }
 range = nrange;
 axmin = naxmin;
 axmax = naxmax;

 ntics = (int)(range / ticinc + 1.000001); 

 nolabltic = 0;
 lablint = 10;
 if (ntics >= 11) {
    hitics = 1;			/* use big and small tics */
    if (ntics<50) {
      lablint = 5;
      nolabltic = 1;
      if (ntics<30) nolabltic = 0;
    }
 }
 else {
    hitics = 0;
    ticinc/=10.0; 
    ntics = (int)(range / ticinc + 1.000001); 
 }
 floormul = 1.000001;
 itic   = (int)((axmin * floormul) / ticinc);
 ticmin = itic * ticinc;  
 plticmin = ((ticmin - axmin) / range) * plaxis;   /* physical loc of xmin */
 plticmin += xplotmargl;
 pltic = (ticinc / range) * plaxis;

 calibtic = ticinc * 10.0;

/* if (hitics) calibtic = ticinc * 10.0;
 else        calibtic = ticinc; */
 
 calibms = 0;
 xcalstr = (char*)NULL;
 switch (plotnod[n].pmod) {

       default:
        if (calibtic < 1) {                     /* x axis calibration units */
             calibms = 1; 
             xcalstr = "msec";
        }
        else xcalstr = "Seconds";
        break;
 
    case  GREC:
         xcalstr = "";
      break;

   } /* switch (plotnod[n].pmod) */
 
 for (itici=itic,j=0; j<ntics; j++,itici++) {     /* draw x axis tics */
   move (plticmin, ybase);
   if ((itici % lablint) == 0) 
        draw (plticmin, ybase - LARGETIC);    /* always draw larger tics */
   else 
   if (hitics) draw (plticmin, ybase - SMALLTIC); /* smaller tics */

   if ((itici % lablint) == 0){
     if (nolabltic && (abs(itici%10) == 5)) ;
     else {
        cpen (WHITE);
        pxnum = ((int)(axmin/ticinc * 1.000001) + j) * ticinc;
        if (ncabs(pxnum) < (range / 100)) pxnum = 0;
        if (calibms) pxnum *= 1000;
	move(plticmin,ybase-XPLOTLABEL-charsiz*1.2); /* position of label */
	prnum("c",charsiz,"%g",pxnum,0,0);   	 /* label x tics */
        cpen (BLUE);
     }
   }
   plticmin += pltic;
 } 

 cpen (WHITE);					/* label the end of axes */

 pxnum = axmin;
 if (ncabs(pxnum) < (range / 100)) pxnum = 0;
 if (calibms) pxnum *= 1000;
 if (hitics && ntics < 61) {
    if ((itic % lablint) < lablint/2 && (itic%lablint) > 0) {
      move (xplotmargl,ybase-XPLOTLABEL-charsiz*1.2);
      prnum("c",charsiz,"%g",pxnum,0,0);	/* label xmin */
    }
 }
 
 pxnum = axmax; 
 if (ncabs(pxnum) < (range / 100)) pxnum = 0;
 if (calibms) pxnum *= 1000;
 if (hitics && ntics < 61) { 
    if ((--itici % lablint) > lablint/2) {
 	move (1-xplotmargr, ybase-XPLOTLABEL-charsiz*1.2); 
 	prnum("c",charsiz,"%g",pxnum,0,0); 		/* label xmax */
    }
 }

 if (!finfl) {
  if (plotnod[n].pmod != GREC) {
   move ((1+xplotmargl-xplotmargr)*0.5, ybase-XPLOTLABEL-charsiz*3);
   prtext ("c",charsiz,xcalstr);                /* draw x axis calib label */
  }
 }
 if (*plframe) frame ("..");

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 /* draw y axes */

 for (i=0; i<PLOTNODSIZ; i++) {
   plotdone[i] = 0;
   yaplot[i] = 0;
   ybplot[i] = 0;
 }
 for (m=i=0; m<numyaxis; i++,m++) {        /* draw (maybe multiple) y axis */
   rightYax = 0;
   if (useframe) n = looklist(nplot);
   else if (plsep) {
      n = looklist(i);  /* index for axis color, data */
   }
   else if (i>=numlaxis) {
      rightYax = 1;		/* do a right axis */
      n = looklist(lookaxes(i-numlaxis, rightaxes)); /* for axis color, data */
   }
   else {
      rightYax = 0;		/* do a left axis */
      n = looklist(lookaxes(i, leftaxes));
   }

   plframe = plotnod[looklist(n)].plframe;	/* change to correct frame */
   if (*plframe) frame (plframe);
   cwidth (charsiz,".");


   if (plsep || i || setaxis) {
      pen = plotnod[n].ppen;
      if (pen == NULLVAL) pen = colortab[n];
   }
   else {
     pen = plotnod[looklist(nplot)].ppen;
     if (pen == NULLVAL) pen = colortab[0];
   }
   if (plsep || numyaxis>ttotplots) axpen=pen;		/* color of y axis */
   else                             axpen=BLUE;		/* color of y axis */
   cpen(axpen);

   aplot = getaplot(i);

   if (plotdone[aplot]) {			/* only do y-axis once */
      m--;
      continue;
   }
   else plotdone[aplot] = 1;

	/* at this point in plsep mode, we're only drawing active axes */

   if (plsep) { 		
					/* yaxis mod by plotsize */
      yaxis = splots[n] / ttotplots * (1-yplotmarg*1.5) -
					YAXISSPC*10.0/(ttotplots+10.0);
      ybase = scplots[n] / ttotplots * (1-yplotmarg*1.5) + yplotmarg;
      yaplot[aplot] = yaxis;
      ybplot[aplot] = ybase;

/* ncfprintf (stderr,"i %-2d n %-2d m %-2d aplot %-2d yax %-6.3g yb %-6.3g plotn %-2d pen %-2d size %g\n",
		i,n,m,aplot, yaxis, ybase,plotnod[n].plotn,
		plotnod[n].ppen, splots[n]); 			/* */
   }
   else {
      yaxis = (1.0 - yplotmarg) / ttotplots - YPLOTMARG * numyaxis;
      ybase = yplotmarg;
   }

   if (numyaxis>ttotplots) {
     if (rightYax) yaxisX = 1 - xplotmargr + (aplot-numlaxis)*YAXISWIDTH; /* R */ 
     else          yaxisX = xplotmargl - aplot * YAXISWIDTH;	/* on left */
   }
   else	       yaxisX = xplotmargl;

   if (axpen > 0) {
    move (yaxisX, ybase);			/* draw y axis */
    draw (yaxisX, ybase+yaxis);
   }

   plaxis = yaxis;                              /* size of plot y axis */
   axmax = plotnod[n].pymax;
   axmin = plotnod[n].pymin;
   range = axmax - axmin;
   if (range == 0) range = 1e-20;
   plotnod[n].yrange = range;

   /* The problem here is to find the best tick spacing. */
   /* In "plsep" mode we need to take into account both the range */
   /*  and the size of the plot */

   if (plsep) {					/* compute scaled range */
      //srange = range/(naplots*.7*yaplot[aplot]); /* find y axis tic incr */
      srange = splots[n] * range/(naplots * 0.7*yaplot[aplot]); /* find y axis tic incr */
   }
   else
      srange = range;
   ticinc = geticinc(srange);			/* find y axis tic incr */

   if (i==0 && plotnod[n].automax) {		/* rescale y axes */
       int k;

     if (axmin >= 0) floormul = 1.000001;
     else            floormul = 0.999999;
     itic   = (int)floor((axmin * floormul) / ticinc);
     ticmin = itic * ticinc;  
     naxmin = ticmin;
     if (axmax < 0) floormul = 1.000001;
     else           floormul = 0.999999;
     itic   = (int)ceil((axmax * floormul) / ticinc);
     naxmax = itic * ticinc;  
     nrange = (naxmax - naxmin);
     plotnod[n].yrange = nrange;		/* first y axis */
     plotnod[n].pymax = naxmax;
     plotnod[n].pymin = naxmin;

     for (k=1; k<totplots; k++) {		/* other y axes */
       double trange;
       int p;

       p = looklist(k); 
       trange = plotnod[p].yrange;
       plotnod[p].pymax += (naxmax-axmax)/range*trange;
       plotnod[p].pymin += (naxmin-axmin)/range*trange;
       plotnod[p].yrange = plotnod[p].pymax - plotnod[p].pymin;
     }
   range = nrange;
   axmin = naxmin;
   axmax = naxmax;
   }  /* if i */

   if (axpen <= 0) continue;
					/* find tics as if plot fills display */
   ntics = (int)(srange / ticinc + 1.000001);	
   nolabltic = 0;
   lablint = 10;
   if (ntics >= 11) {
      hitics = 1;		/* use big and small tics */
      if (ntics<50) {
        lablint = 5;
        nolabltic = 1;
        if (ntics<30) nolabltic = 0;
      }
   }
   else {
      hitics = 0;
      ticinc/=10.0; 
      ntics = (int)(range / ticinc + 1.000001); 
   }
   ntics = (int)(range / ticinc + 1.000001); 	/* find actual # of tics */
   itic  = (int)((axmin * 1.000001) / ticinc);
   ticmin = itic * ticinc;

   plticmin = ((ticmin - axmin) / range) * plaxis;
   plticmin += ybase;
   pltic = (ticinc / range) * plaxis;

   calibtic = ticinc * 10.0;

   /* if (hitics) calibtic = ticinc * 10.0;
   else        calibtic = ticinc; */

   calibmv = 0;
   calstr = (char *)NULL;
   switch (plotnod[n].pmod) {

    case NRECGV: 
    case CREC_CAV: 
    case  VREC:
    case  WREC:
    case  CABLE:
        if (calibtic < 1) {                     /* y axis calibration units */
             calibmv = 1; 
             if (nonodes) calstr = "mVolts";
             else         calstr = "mV";
        }
        else {
             calibmv = 0;
             if (nonodes) calstr = "Volts";
             else         calstr = "V";
        }
        break;
 
    case NRECGI: 
    case CREC_CAIT: case CREC_CAIP: case CREC_CAIE: case CREC_CAIPE: case CREC_CICR:
    case MREC: 
    case IREC: 
          if (nonodes) calstr = "Amperes";
          else         calstr = "A";
      break;

    case  LREC:
          if (nonodes) calstr = "Photons/sec";
          else         calstr = "Ph";
      break;

    case CREC_GLU:  case CREC_AMPA: case CREC_KAINATE: case CREC_NMDA: case CREC_CNQX:
    case CREC_GABA: case CREC_BIC:  case CREC_PTX:  case CREC_GLY: case CREC_CHRC:
    case CREC_STRY:  case CREC_CAMP: case CREC_CGMP:  case CREC_ATP:

   case CREC_CABUF:
   case CREC_CACONC:
   case CREC_CAS:
        if (nonodes) calstr = "Molar";
        else         calstr = "M";
	break;

    case CREC_PH:
        if (nonodes) calstr = "pH";
        else         calstr = "pH";
	break;

    case NRECA0: case NRECA1: case NRECA2: case NRECA3: case NRECA4:
    case NRECH0: case NRECH1: case NRECH2: case NRECH3: case NRECH4:
       calstr = "Frac";
	break;
    case NRECB0: case NRECB1: case NRECB2: case NRECB3: case NRECB4:
        if (nonodes) calstr = "Molar";
        else         calstr = "M";
	break;
    case NRECC0: case NRECC1: case NRECC2: case NRECC3: case NRECC4:
    case NRECC9:
    case NRECG1: case NRECG2: case NRECG3: 
    case NRECGC: case NRECGM: case NRECGH: 
       calstr = "Frac";
	break;
    case NRECA8:
    case NRECA9:
        if (nonodes) calstr = "Rate";
        else         calstr = "R";
	break;

    case NRECG: if (plotnod[n].pval>0) {
                 calstr = "Frac"; 
		 break;
	      }
			/* no break, if =0, go to next case */	
    case NRECG0:
        if (nonodes) calstr = "Siemens";
        else         calstr = "S";
	break;

    case  FREC: 
	  switch (plotnod[n].pmod2) {
	    case IREC:
	    case NRECGI:
    	    case MREC: 
    	    case CREC_CAIT: case CREC_CAIP: case CREC_CAIE: case CREC_CAIPE: case CREC_CICR:
               if (nonodes) calstr = "Amperes";
	       else         calstr = "A";
	       break;
	    case NRECG0: case NRECG1: case NRECG2: case NRECG3: case NRECG:
		if (nonodes) calstr = "Siemens";
		else         calstr = "S";
		break;
 	  }
	  break;
	  
    case  SREC:
    case  GREC:
         calstr = "";
      break;

   } /* switch (plotnod[n].pmod) */

   ticymin = ybase- ycharsiz/2;
   ticymax = ybase+yaxis- ycharsiz/2;
   ticyminfl=ticymaxfl=0;
   if (rightYax) ticdir = -1;
   else          ticdir = 1;
   for (itici=itic,j=0; j<ntics; j++,itici++) {   /* draw y axis tics */
     move (yaxisX, plticmin);
     if ((itici % lablint) == 0)
          draw (yaxisX-LARGETIC*ticdir, plticmin); /* always draw larger tics */
     else
      if (hitics) draw (yaxisX-SMALLTIC*ticdir, plticmin);  /* smaller tics */

     if ((itici % lablint) == 0) {
       if (ttotplots < 20) {			/* arbitrary limit on label printout */
         if (nolabltic && (abs(itici%10) == 5)) ;
         else {
	   cpen (WHITE);
           pynum = ((int)(axmin/ticinc * 1.000001) + j) * ticinc;
           if (ncabs(pynum) < (range / 100)) pynum = 0;
           if (calibmv) pynum *= 1000.0;
	   ticypos = plticmin- ycharsiz*.4;
	   if (ticypos==ticymin) ticyminfl = 1;
	   if (ticypos==ticymax) ticymaxfl = 1;       /*label y tics*/
           move (yaxisX-(LABELTIC+ycharsiz*LABELOFFS)*ticdir,ticypos);
	   if (ticdir>0) prnum ("r",ycharsiz,"%g",pynum,0,0);
	   else          prnum ("l",ycharsiz,"%g",pynum,0,0);
           cpen (axpen);
         }
       }
     }
     plticmin += pltic;
   } 
   ppmax = plotnod[n].pymax;
   ppmin = plotnod[n].pymin;
   if (calibmv) {
        ppmax = plotnod[n].pymax * 1000;
        ppmax = int(ppmax + 0.5);
        ppmin = plotnod[n].pymin * 1000;
        ppmin = int(ppmin + 0.5);
   }
   plmaxlen = getlengthd (ppmax);
   plminlen = getlengthd (ppmin);
   pllenmax = max(plmaxlen,plminlen);
   pllenmax = min(pllenmax,6);		/* X location of max calib label */
   pllenmax += 1;

   if (!ticyminfl) {
     cpen (WHITE);

     if (hitics && ntics < 61) {
        if ((itic % lablint) < lablint/2 && (itic%lablint) > 0) {
     	  move (yaxisX-(LABELTIC+ycharsiz*LABELOFFS)*ticdir, ticymin);
	  if (ticdir>0) prnum ("r",ycharsiz,"%.2g",ppmin,0,0); /* display ymin*/
	  else          prnum ("l",ycharsiz,"%.2g",ppmin,0,0);
        }
     }

     if (!nonodes && !plsep) {
       move (yaxisX-(LABELAX+ycharsiz*LABELOFFS)*ticdir,ticymin+0.035);
       if (ticdir>0) prtext ("r",ycharsiz,calstr); /* draw y axis calib label */
       else          prtext ("l",ycharsiz,calstr); /*   near bottom of Y axis */
       cpen (axpen);
     }
   }

   if (!ticymaxfl) {

     cpen (WHITE);

     if (hitics && ntics < 61) {
        if ((--itici % lablint) > lablint/2) {
	  move (yaxisX-(LABELTIC+ycharsiz*LABELOFFS)*ticdir, ticymax);
	  if (ticdir>0) prnum ("r",ycharsiz,"%.2g",ppmax,0,0); /* display ymax*/
	  else          prnum ("l",ycharsiz,"%.2g",ppmax,0,0);
        }
     }

    if (!finfl) { 
     if (nonodes) { 
       /*move (yaxisX-LABELAX*ticdir,ticymax+YLABELPOSN); /* At top of Y axis,*/
       move (yaxisX,ticymax+YLABELPOSN); /* At top of Y axis, */
       if (ticdir>0) prtext ("r",ycharsiz,calstr); /* draw y axis calib label */
       else          prtext ("l",ycharsiz,calstr);
     }
     else {
       //move (yaxisX-(+LABELAX*ticdir),ticymax-YLABELPOS); 
       move (yaxisX-(pllenmax*ycharsiz*ticdir),ticymax); 
       if (ticdir>0) prtext ("r",ycharsiz,calstr); /* draw y axis calib label */
       else          prtext ("l",ycharsiz,calstr);
     }
    }
     cpen (axpen);
   }
   if (*plframe) frame ("..");

 }              /* for (i; i<ttotplots; ) */

if (!nonodes) {
 for (i=0; i<PLOTNODSIZ; i++) {
    plotlabl[i] = 0;
    nplots[i] = 0;
 }
 if (plsep)
  for (m=i=0; m<totplots; i++,m++) {		/* draw plot labels */
      int pln;

   n = looklist(i); 
   aplot = getaplot(i);

   plframe = plotnod[looklist(n)].plframe;	/* change to correct frame */
   if (*plframe) frame (plframe);
   cwidth (charsiz,".");

   pen = plotnod[n].ppen;
   if (pen == NULLVAL) pen = colortab[n];
   if (pen <=0) continue;
   cpen (pen);					/* change pen for node num */

   ybase = ybplot[aplot]-charsiz*.5;
   yaxis = yaplot[aplot];

   if (plotnod[n].plname[0]) {			/* plot name label */
      nplots[aplot]++;				/* count of traces on plot */
      plotlabl[aplot] = 1;
						/* display plot name */
      move (XABSMARGL, ybase+yaxis/2+(nnplots[aplot]*.5-nplots[aplot]+.5) * 
						yplcharsiz[aplot]*1.25);
      prtext ("l",yplcharsiz[aplot],plotnod[n].plname);
   }
   else {

    if (plotlabl[aplot]) {		    /* if plot has already been named */
      continue;
    }
    nplots[aplot]++;				/* count traces in plot */
						/* display node numbers */
    move(XABSMARGL, ybase+yaxis/2+(nnplots[aplot]/2-nplots[aplot]+1) *
						yplcharsiz[aplot]*1.16); 

    if (plotnod[n].pmod==GREC) {
      rmove (0.05, 0);                              /* display node number */
      prnum ("r",yplcharsiz[aplot],"%2g",
   		(double)n,0,0); 	         /* just plot number */
    }
    else if (plotnod[n].cnod3 > NULLVAL) {
      prnum ("l",yplcharsiz[aplot],"%g %g %g",
		(double)plotnod[n].cnod1,
		(double)plotnod[n].cnod2,
		(double)plotnod[n].cnod3);	 /* (display cramped cnod3) */
    }
    else if (plotnod[n].cnod2 > NULLVAL) {
      prnum ("l",yplcharsiz[aplot],"%2g %2g",
		(double)plotnod[n].cnod1,
		(double)plotnod[n].cnod2,0);	 /* (no room for cnod3) */
    }
    else
      prnum ("l",yplcharsiz[aplot],"%2g",
		(double)plotnod[n].cnod1,0,0); 	 /* (no room for cnod2) */
    } 
   if (*plframe) frame (",.");
 }
 else {     /* !plsep */

  for (i=0; i<totplots; i++) {                   /* draw node nums */
    n = looklist(i); 
    plframe = plotnod[looklist(n)].plframe;	/* change to correct frame */
    if (*plframe) frame (plframe);
    cwidth (charsiz,".");
    pen = plotnod[n].ppen;
    if (pen == NULLVAL) pen = colortab[n];       /* default pen is plot num */
    cpen (pen);                                  /* change pen for node num */
    ybase = yplotmarg;
    move (XABSMARGL,ybase+yaxis/2+(totplots/2-i)*charsiz*1.8); /* disp node num*/
  if (plotnod[n].plname[0]) {
     prtext ("l",charsiz*1.4,plotnod[n].plname);
  }
  else 
  if (plotnod[n].pmod==GREC) {
   rmove (0.05, 0);
   prnum ("r",charsiz,"%2g",
   		(double)n,0,0); 	         /* just plot number */
  }
  else if (plotnod[n].cnod3 > NULLVAL) {
   prnum ("l",charsiz,"%2g %g %g  ",
		(double)plotnod[n].cnod1,
		(double)plotnod[n].cnod2,
		(double)plotnod[n].cnod3);	 /* (display cramped cnod3) */
  }
  else if (plotnod[n].cnod2 > NULLVAL) {
   prnum ("l",charsiz,"%2g %g  ",
		(double)plotnod[n].cnod1,
		(double)plotnod[n].cnod2,0);	 /* (no room for cnod3) */
  }
  else
   prnum ("l",charsiz,"%2g    ",
   		(double)plotnod[n].cnod1,0,0); 	 /* (no room for cnod3) */
   if (*plframe) frame (",.");
  }      /* for (i;;) */
 }     /* else !plsep */
}    /* if (!nonodes) */

 if (!nofilfl) {
   plframe = plotnod[looklist(nplot)].plframe;	/* change to correct frame */
   if (*plframe) frame (plframe);
   draw_plotlabel(1.0);
   if (*plframe) frame ("..");
 }

 gpurge();

#ifdef DEBUG 
 if (debug & NCPLOT) ncfprintf (stderr,"drawframe_end\n");
#endif

}

/*------------------------------------*/

void prnum (const char *adj,double csize,const char *fmt,
		double val1,double val2,double val3) 

/* print a number, formatted with either
    left, right, or center justification. */

{
   char numbuf[30];

  sprintf (numbuf,fmt,val1,val2,val3);
  prtext (adj,csize,numbuf);
}

/*------------------------------------*/

void prtext (const char *adj, double csize, const char *str)

/* print a string, formatted with either
    left, right, or center justification. */

{
    double len,xoffs,yoffs;
    int large,small,caps;
    const char *p;

  if (!str) return;
  yoffs = 0;
  if ((len=strlen(str))==0) return;   		   /* length of number */
  large=small=caps=0;
  for (p=str;*p; p++) {
     if (*p=='m') large++;
     else if (*p=='w') large++;
     else if (*p=='i') small++;
     else if (*p=='.') small++;
     else if (*p==',') small++;
     else if (*p=='-') small++;
     else if (*p==':') small++;
     else if ((*p>='A')&& (*p<='Z')) caps++;
  }
  len += (large+caps)*0.7 - (small*0.5);   /* proportional spacing */

  switch (*adj) {

    case 'c': xoffs=len/2.0;    	   /* center justification */
   	     break;
    case 'l': xoffs=0;    		   /* left justification */
	     break;
    case 'r': xoffs=len;		   /* right justification */
	     break;
    default: ncfprintf (stderr,"prtext: bad adjust '%c'\n",*adj);
	      return;
	     break;
  }
  xoffs *= csize * 0.78;
  gcwidth (csize);
  grmove (-xoffs, yoffs);
  gtext (str);
  gcwidth (charsiz);
}

/*------------------------------------*/

int looklist(int i)
         

/* lookup plot in command line list and possibly
    return different plot number.
    This rearranges order of displayed plots.
*/

{
    int retval;
    extern int listplots[PLOTNODSIZ];

   if (i<0) i = 0;
   retval = i;
   if (listplots[0]) {
	retval=listplots[i]-1;	 		/* possibly different order */
  	if (retval<0) retval = 0;		/*   from command line */
	if (retval>=dttplots) retval = dttplots - 1;
   }
  if (retval>=PLOTNODSIZ) retval = PLOTNODSIZ-1;
  return retval;
}

/*------------------------------------*/

int lookaxes(int i, int list[])
         

/* lookup plot in command line list and possibly
    return different plot number.
    This rearranges order of displayed plots.
*/

{
   int j;

   if (i<0) i = 0;
   if (list[0]) {
	j=list[i]-1;	 	/* possibly different order */
  	if (j<0) j = 0;		/*   from command line */
	if (j>=dttplots) j = dttplots - 1;
   }
  else j=i;
  return j;
}

/*------------------------------------*/

void getplot(FILE *instream, int numpl, int setauto)

/* reads in data from file for plotting */

{
#define INBUFSIZ 10000

   static char *inbuf,tfil[20],vbuf[40];
   int i,f,k,n,ttplots,pmod,pmod2,pval,pen,plotn,np,npl,nplot,narg;
   int xpos,lchar,charmode,graphmem;
   char *str, *ptr, pmodc[10], plframe[40];
   double val, xval,csize,plotsiz;
   double x,y;
   static int automode = 1;
   double *plotarr,*plotp;

 if (!(inbuf = (char *)emalloc (INBUFSIZ))) {
	ncfprintf (stderr,"Error, getplot: can't allocate line buffer\n");
        return;
 }

  if (setvid) vidmode = 1;
  else if (unsetvid) vidmode = 0;

    if (rowmode) {	/* if rowmode, allocate memory to swap X and Y */

#define GRAPHMEM (PLOTNODSIZ*100000)

    graphmem = GRAPHMEM;

     if (!(plotarr = (double *)emalloc (graphmem * sizeof(double)))) {
       graphmem /= 4;
       if (!(plotarr = (double *)emalloc (graphmem * sizeof(double)))) {
   	   ncfprintf (stderr,"Error, getplot: can't allocate memory\n");
           return;
       }
     }
    xpos = 0;
   }
  frame (rootframe); 
  for (npl= -1; fgets(inbuf,INBUFSIZ,instream); ) {
   if (inbuf[0] == '#') {
     if (!setauto)
      switch (inbuf[1])  {

       case 'c':
                strcpy (complin,inbuf);
		break;
       case 'd':
                sscanf(inbuf,"%*s %d %*s %d\n",&pen,&nplot);
                if (nplot>=PLOTNODSIZ) nplot = PLOTNODSIZ-1;
                if (pen == NULLVAL) pen = colortab[nplot];
		plotpen(pen,nplot);
                break;
       case 'l':
               sscanf(inbuf,"%*s '%c' %*s %c %*s %d\n",(char*)&lchar,(char*)&charmode,&nplot);
                if (nplot>=PLOTNODSIZ) nplot = PLOTNODSIZ;
		plotchar(lchar,charmode,nplot);
                break;
       case 'z':
                sscanf(inbuf,"%*s %lg %*s %d\n",&csize,&nplot);
                if (nplot>PLOTNODSIZ) nplot = PLOTNODSIZ-1;
		plotcsiz(csize,nplot);
                break;
       case 'f':
                narg = sscanf(inbuf,"%*s %s %*s\n",tfil);
		if (strcmp(tfil,"stimfile") == 0) {
                   narg = sscanf(inbuf,"%*s stim file %60s\n",stfilnam);
		   stfile = stfilnam;
		}
		else if (strcmp(tfil,"input") == 0) {
                   narg = sscanf(inbuf,"%*s input file %60s\n",infilnam);
/*		   infile = infilnam;     /* ignore "infile" here, */
		}			  /*  set infile from command line */
		else if (strcmp(tfil,"source") == 0) {
                   narg = sscanf(inbuf,"%*s source file %60s\n",srcfilnam);
		   srcfile = srcfilnam;
		}
		else if (strcmp(tfil,"plotlabel") == 0) {
                   narg = sscanf(inbuf,"%*s plotlabel %100s\n",plotlabelnam);
		   plotlabel = plotlabelnam;
		}
                break;

       case 'g':

		switch (inbuf[2]) {
   
    		case 'c': switch(inbuf[3]) {

		  case 'i': sscanf (inbuf,"%*s %lg%*[ ,]%d",&x,&f); /* circle */
			  gcirc(x,f);
			  break;

		  case 'r': sscanf (inbuf,"%*s  %lg", &x);
			  gcrotate (x);
			  break;

		  case 't': gctext ();			/* text mode */
			  break;			/* antiquated */

		  case 'g': gcgraphics ();		/* graphics mode */
			  break;			/* antiquated */

		   case 'w': sscanf (inbuf,"%*s  %lg", &x);
			  gcwidth(x);	
			  break;
			}
			break;
		
    		case 'd': switch(inbuf[3]) {

    			case ' ':			/* draw */
			 sscanf (inbuf,"%*s %lg%*[ ,]%lg", &x,&y);
			 gdraw (x,y);
			 break;

			case 'a':			/* dash */
			 sscanf (inbuf,"%*s %lg", &x);
   			 gdash ((int)x);
			 break;

			}  /* case "gd?" */
			break;
	
    		case 'f': sscanf (inbuf,"%*s %s", tfil);	/* frame */
			  gframe (tfil);
			  break;

    		case 'h': switch(inbuf[3]) {

    			case 'i': ghinit ();		/* hidinit */
			 break;

			case 's': ghstart();		/* hidstart */
			 break;

			case 'e': ghend();		/* hidstop */
			 break;

			}  /* case "gh?" */
			break;
	
    		case 'm': sscanf (inbuf,"%*s  %lg%*[ ,]%lg", &x,&y);
			  gmove (x,y);
			  break;

    		case 'o': sscanf (inbuf,"%*s  %lg%*[ ,]%lg", &x,&y);
			  gorigin (x,y);
			  break;

    		case 'r': switch(inbuf[3]) {

			  case 'd':
			   sscanf (inbuf,"%*s  %lg%*[ ,]%lg", &x,&y);
			   grdraw (x,y);
			   break;

			  case 'e': {				/* rectangle */
				double x1,y1,x2,y2,x3,y3,x4,y4;
				int fill;
	
			  sscanf (inbuf,
  "%*s  %lg%*[ ,]%lg%*[ ,]%lg%*[ ,]%lg%*[ ,]%lg%*[ ,]%lg%*[ ,]%lg%*[ ,]%lg%*[ ,]%d",
			 &x1,&y1,&x2,&y2,&x3,&y3,&x4,&y4,&fill);
			    grect (x1,y1,x2,y2,x3,y3,x4,y4,fill);
			   }
			   break;

			  case 'f':
    			   sscanf (inbuf,"%*s %s", tfil);	/* rmframe */
			   grmframe (tfil);
			   break;

			  case 'm':
			   sscanf (inbuf,"%*s  %lg%*[ ,]%lg", &x,&y);
			   grmove (x,y);
			   break;
    		
			  case 'o':
			   sscanf (inbuf,"%*s  %lg", &x);
			   grotate (x);
			  break;

			 }
			 break;  /* case "gr?" */

    		case 'p': switch(inbuf[3]) {

    			case 'a': gpage ();
			  	  break;
    			case 'e': sscanf (inbuf,"%*s %lg", &x);
			  	  gpen ((int)x);
			  	  break;
    			case 'u': gpurge ();
			  	  break;

		}   /* case "gp?" */
		break;

    		case 's': sscanf (inbuf,"%*s %lg", &x);
			  gsize (x);
			  break;

    		case 't': switch(inbuf[3]) {

    		  case ' ': i = strlen (&inbuf[4]);
			  inbuf[i+3] = 0;
			  gtext (&inbuf[4]);
			  break;
		  case 'r': {				/* triangle */
				double x1,y1,x2,y2,x3,y3;
				int fill;
	
			  sscanf (inbuf,
  "%*s  %lg%*[ ,]%lg%*[ ,]%lg%*[ ,]%lg%*[ ,]%lg%*[ ,]%lg%*[ ,]%d",
			 &x1,&y1,&x2,&y2,&x3,&y3,&fill);
			    gtri (x1,y1,x2,y2,x3,y3,fill);
			   }
			  break;

		} /* case "gt?" */
		break;

		}  /* switch (inbuf[2] */
		break;		/* case "g?"  (graphics primitives */

      case 'p':
                narg = sscanf(inbuf,"%*s %d %*s %*s %d %*s %lg %*s %lg\n",
				&dttplots,&plsep,&plmax,&plmin);
		if (narg<4) {
                  narg = sscanf(inbuf,"%*s %d %*s %*s %lg %*s %lg\n",
				&dttplots,&plmax,&plmin);
		  plsep = 0;
		}
                if (narg<3) {
		    plmax =  .04;
		    plmin = -.07;
		}
		if (numpl) ttplots = numpl;
                else ttplots = dttplots;
		if (setplsep) plsep = !plsep;
                break;
      case 'e':
                narg = sscanf(inbuf,"#e %*s %lg %*s %lg %*s %d %*s %s\n",
				&plotxmin,&plotxmax,&rseed,(char*)&vbuf);
		ncversion = scanversion(vbuf);

		if (narg < 4) {				/* old version */

                   narg = sscanf(inbuf,"#e %lg %*s %*s %d %*s %lg\n",
				&plotxmax,&rseed,&ncversion);
		   plotxmin = 0;
		}
		endexp = plotxmax;
		simtime = plotxmin;
                break;
      case 'x':
	      np = npl + 1;		/* y value increments plot number */
              narg = sscanf(inbuf,"#x %*s %d %d %d %d %*s %s %*s %lg %*s %lg\n",
                  &plotnod[np].cnod1,&plotnod[np].cnod2,
                  &plotnod[np].cnod3,&plotnod[np].cnod4,
		  pmodc,&plotnod[np].pxmax,&plotnod[np].pxmin);

               if (narg != 5)
		if (narg<7) {
                 narg = sscanf(inbuf,"#x %*s %d %d %d %*s %s %*s %lg %*s %lg\n",
                  &plotnod[np].cnod1,&plotnod[np].cnod2, &plotnod[np].cnod3,
		  pmodc,&plotnod[np].pxmax,&plotnod[np].pxmin);
		}

		  if (!notinit(setxmax) && !notinit(plotnod[np].pxmax)) 
				 plotnod[np].pxmax = setxmax;
		  if (!notinit(setxmin) && !notinit(plotnod[np].pxmin)) 
				 plotnod[np].pxmin = setxmin;

		if (*(pmodc+1)) {
			char pmodchar2;

		  pmodchar2 = *(pmodc+1);
		  if (pmodchar2=='I') pmod2 = IREC;
		  else if (pmodchar2=='G') pmod2 = NRECG0;
		  else if (pmodchar2>='0' && pmodchar2 <= '9') {
		 	sscanf (pmodc+1,"%d",&pval);
		  }
		}
		else { pval = 0; pmod2 = 0; }

                  switch (*pmodc) {
                    case 'V':   pmod = VREC;   break;      
                    case 'W':   pmod = WREC;   break;      
                    case 'M':   pmod = MREC;   break;      
                    case 'I':   pmod = IREC;   break;      
                    case 'L':   pmod = LREC;   break;      
                    case 'G':   pmod = NRECG; break;      
                    case 'C':   pmod = CREC_CAV+pval; break;      
                    case 'F':   pmod = NRECA0; break;      
                    case 'R':   pmod = NRECA9; break;      
                    case 'N':   pmod = CREC_GLU+pval; break;      
                    case 'U':   pmod = FREC;   break;      
                    case 'S':   pmod = SREC;   break;      
                    case 'X':   pmod = GREC;   break;      
                    case 'Y':   pmod = CREC_CACONC; break;
                  }
                plotnod[np].pmod  = pmod;
                plotnod[np].pmod2 = pmod2;
                plotnod[np].pval  = pval;
                if (plotnod[np].cnod1 == NULND) plotnod[np].cnod1 = NULLVAL;
                if (plotnod[np].cnod2 == NULND) plotnod[np].cnod2 = NULLVAL;
                if (plotnod[np].cnod3 == NULND) plotnod[np].cnod3 = NULLVAL;
                if (plotnod[np].cnod4 == NULND) plotnod[np].cnod4 = NULLVAL;
                break;
      case 'y':
         if (npl <  dttplots-1) npl++;		/* increment plot number */
         plotn = 0;
	 plotsiz = 1;
         narg=sscanf(inbuf, "#y fr %s ", plframe);

	if (narg==1) {
         strcpy(plotnod[npl].plframe,plframe);  /* copy frame name */
						/* look for char mode */
         narg=sscanf
       (inbuf,
  "#y fr %*s pl %d siz %lg pen %d %*s %c %*s '%c' %*s %lg %*s %lg %*s %lg\n",
			&plotn,&plotsiz,&pen,
			&plotnod[npl].charmode,&plotnod[npl].charfl,
			&plotnod[npl].csize,
			&plotnod[npl].pymax, &plotnod[npl].pymin);

	if (narg < 8)			/* no char */
        narg = sscanf(inbuf,
 		 "#y fr %*s pl %d siz %lg pen %d styl %c ymax %lg ymin %lg\n",
			&plotn, &plotsiz, &pen, &plotnod[npl].charmode,
			&plotnod[npl].pymax, &plotnod[npl].pymin);
	if (narg < 6)			/* no char */
           narg=sscanf(inbuf,
		"#y fr %*s pl %d siz %lg pen %d mode %c ymax %lg ymin %lg\n",
			&plotn, &plotsiz, &pen, &plotnod[npl].charmode,
			&plotnod[npl].pymax, &plotnod[npl].pymin);

	if (narg < 5)			/* no max min */
                  narg = sscanf(inbuf,"#y fr %*s pl %d siz %lg pen %d styl %c\n",
			&plotn, &plotsiz, &pen, &plotnod[npl].charmode);

		if (narg < 2) {			/* old version */
                  narg = sscanf(inbuf, "#y fr %*s pen %d ymax %lg ymin %lg\n",
                      &pen, &plotnod[npl].pymax, &plotnod[npl].pymin);
		  plotnod[npl].charmode = PLINES;
		 }
	   }

		else {		/* no frame */

	 	narg=sscanf
	    (inbuf, 
	"#y pl %d siz %lg pen %d %*s %c %*s '%c' %*s %lg %*s %lg %*s %lg\n",
			&plotn,&plotsiz,&pen,
			&plotnod[npl].charmode,&plotnod[npl].charfl,
			&plotnod[npl].csize,
			&plotnod[npl].pymax, &plotnod[npl].pymin);

		if (narg < 8)			/* no char */
               narg = sscanf(inbuf, 
			"#y pl %d siz %lg pen %d styl %c ymax %lg ymin %lg\n",
			&plotn, &plotsiz, &pen, &plotnod[npl].charmode,
			&plotnod[npl].pymax, &plotnod[npl].pymin);

		if (narg < 6)			/* no char */
               narg = sscanf(inbuf, 
			"#y pl %d siz %lg pen %d mode %c ymax %lg ymin %lg\n",
			&plotn, &plotsiz, &pen, &plotnod[npl].charmode,
			&plotnod[npl].pymax, &plotnod[npl].pymin);

		if (narg < 5)			/* no max min */
                  narg = sscanf(inbuf, "#y pl %d siz %lg pen %d styl %c\n",
			&plotn, &plotsiz, &pen, &plotnod[npl].charmode);

		if (narg < 2)			/* old version, no pl, no char */
               narg = sscanf(inbuf, "#y pen %d styl %c ymax %lg ymin %lg\n",
			&pen, &plotnod[npl].charmode,
			&plotnod[npl].pymax, &plotnod[npl].pymin);

		if (narg < 4)			/* no max min */
                  narg = sscanf(inbuf, "#y pl %d pen %d styl %c\n",
			&plotn, &pen, &plotnod[npl].charmode);

		if (narg < 4) {			/* old, old version */
                  narg = sscanf(inbuf, "#y pen %d ymax %lg ymin %lg\n",
                      &pen, &plotnod[npl].pymax, &plotnod[npl].pymin);
		  plotnod[npl].charmode = PLINES;
		 }

		}   /* no frame */
		/*
		if (!notinit(setymax) && notinit(plotnod[npl].pymax)) 
				 plotnod[npl].pymax = setymax;
		if (!notinit(setymin) && notinit(plotnod[npl].pymin)) 
				 plotnod[npl].pymin = setymin;
		*/
		plotnod[npl].plotn = plotn;
		plotnod[npl].plotsiz = plotsiz;
		if (pen== -1) pen = NULLVAL;
		plotnod[npl].ppen = pen;
                break;

      case 'r':                         /* restart */
                plotrst (ttplots);
                break;
      case 'n':                         /* node label */
		switch (inbuf[2]) {
		   char ext[20];
#define ST 4
		 case 'o':		/* "node" line */
                  plotinit (ttplots);
		  automode = 0;
                  break;
		 case ' ':		/* node name line */
		  ext[0] = 0;
		  np = npl + 1;
		  if (inbuf[3] == '"') {
		     for (i=ST; i<PNAMSIZ; i++) {
			if (inbuf[i] == '"') break;
		     }
		     strncpy(plotnod[np].plname,(const char *)&inbuf[ST],i-ST);
		  } else {
		    sscanf (inbuf,"%*s  %16s %16s\n", (char*)&plotnod[np].plname,(char*)&ext);
		    if (ext[0]) {
		       strncat (plotnod[np].plname," ",PNAMSIZ);
		       strncat (plotnod[np].plname,ext,
				PNAMSIZ - sizeof(plotnod[np].plname) - 1);
		    }
		  }
		  break;
		}
		break;
#undef ST

    }  /* switch */ 
   }  /* if (inbuf[]==) */ 

  else {
    int j;
    static double plotval[PLOTNODSIZ];
    static int firstline = 1;

  if (automode && firstline) {
      firstline = 0;
      autograph(instream,inbuf,numpl);
  }
  else {
    str=ptr=inbuf;
    val=strtod(str,&ptr);
    xval = val;
    plotnum = npl;

    for (j=k=0; str != ptr && k < PLOTNODSIZ ; k++) {
      str = ptr;
      val=strtod(str,&ptr);
      if (str != ptr) {
         plotval[k] = val;	/* read in plots from one line */
	 j++;
      }
    }
    if (rowmode) {		 /* put out rows, not columns */

      if (numpl)
        for (k=0; k<numpl; k++) {       /* display in possibly different order */
           n = looklist(k);   		/* set by command line */
	   if (n<j) storeplot(plotval[n],xval,ttplots,k,xpos,plotarr); 
      }
      else
        for (k=0; k<j; k++) {	       /* otherwise orig order */
	   storeplot(plotval[k],xval,ttplots,k,xpos,plotarr); 
     }  
     xpos++;				/* next column */

    }
    else {
     if (numpl)
       for (k=0; k<numpl; k++) {       /* display in possibly different order */
           n = looklist(k);   		/* set by command line */
	   if (n<j) mplot(plotval[n],xval,ttplots,k,1); 
     }
     else
       for (k=0; k<j; k++) {	       /* otherwise orig order */
	   mplot(plotval[k],xval,ttplots,k,1); 
       }  
    } /* else ! rowmode */
   }  /* else (!automode) */
  }  /* else */
 }  /* for fgets() */

 if (rowmode) {
     double xval, yval;
     int x, y, xmax, ymax;

   ymax = ttplots; 
   xmax = xpos;
   for (y=0; y<ymax; y++ ) {
     for (x=0; x<xmax; x++ ) {
       xval = recallx (ymax, y, x, plotarr);
       yval = recally (ymax, y, x, plotarr);
       mplot(yval,xval,ymax,y,0); 
    } /* for (; x<xmax; ) */
   } 
 } /* if (rowmode) */
}

/*------------------------------------*/

void storeplot(double yval, double xval, 
		int ymax, int y, int x, double *arr)

{
    if (y==0) *(arr+x*(ymax+1)) = xval;
     *(arr+x*(ymax+1)+ y + 1) = yval;
}

/*------------------------------------*/

double recallx (int ymax, int y, int x, double *arr)

{
   return (*(arr+x*(ymax+1)));
}

/*------------------------------------*/

double recally (int ymax, int y, int x, double *arr)

{
   return (*(arr+x*(ymax+1) + y + 1));
}
/*------------------------------------*/

void autograph(FILE *instream, char *inbuf, int numpl)
                  
/* Plot a graph from a file with multiple lines, 
    each line containing one x value and one or 
    more y values.  Also works for only one number
    per line, in which case the number is assumed to
    be a y value.  In this case, the x value is assumed
    to be the line number.  The "numpl" parameter
    controls how many of the columns (or rows) are plotted.

 */

{
   int i,j,jstart,k,cols,tcols,lines,nox,totplots,pl;
   int graphmem, graphsiz;
   char *str,*ptr;
   double val,xval,yval,ymax,ymin,xmax,xmin;
   double *plotarr,*plotp;

#define GRAPHMEMA (PLOTNODSIZ*100000)

 graphmem = GRAPHMEMA;

 if (!(plotarr = (double *)emalloc (graphmem * sizeof(double)))) {
    graphmem /= 4;
    if (!(plotarr = (double *)emalloc (graphmem * sizeof(double)))) {
	ncfprintf (stderr,"Error, autograph: can't allocate memory\n");
        return;
    }
 }
 ptr=inbuf;
 str = 0;
 plotp = plotarr;
 for (cols=k=0; str != ptr && k < PLOTNODSIZ ; k++) {
   str = ptr;
   val=strtod(str,&ptr);
   if (str != ptr) {
      *plotp++ = val;			/* read in plots from first line */
      cols++;
   }
 }  
 if (k>=PLOTNODSIZ) {
	ncfprintf (stderr,"Error: getplot: %d too many columns  file\n",k);
	return;
 } 
 if (cols==0) {
	ncfprintf (stderr,"Error: getplot: invalid input file, no cols\n");
	return;
 }
 graphsiz = graphmem / cols;

 for (lines=1; fgets(inbuf,INBUFSIZ,instream) && lines<graphsiz; lines++) {
   if (inbuf[0] == '#') { lines--; continue; }
    ptr=inbuf;
    str = 0;
    for (k=0; str != ptr && k < cols ; k++) {
      str = ptr;
      val=strtod(str,&ptr);
      if (str != ptr) {
         *plotp++ = val;		/* read in plots from one line */
      }
      else
         *plotp++ = 0.0;		/* if not enough, then insert a 0 */
    }  
  }  /* for (lines=0; fgets; ) */

  if (lines>=graphsiz) {
	ncfprintf (stderr,"Error: getplot: too many lines %d\n",lines);
	return;
  }

 if (rowmode) {				/* if "rowmode", swap cols with lines */
	tcols = cols;
	cols = lines;
	lines = tcols;
 }
 if (cols>0 && lines>1) {
   if (cols>1) {
      totplots = cols - 1;
      nox = 0;
   }
   else {
     totplots = 1;
     nox = 1;
   } 
   if (setnox) {
      nox = 1;		/* here we want no y val, multiple x values */
      totplots = cols;
   }

   dttplots = totplots;
   if (numpl) totplots = numpl;		/* override default number of plots */

				/* find max, min */
   if (nox) jstart = 0;
   else     jstart = 1;
   for (j=0; j<totplots; j++) {		/* do one column at at time */
     k = j + jstart;
     pl = looklist(j);
     xmax=ymax = -1e30;
     xmin=ymin = 1e30;
     for (i=0; i<lines; i++) {
       if (rowmode) {				/* if rowmode,swap i,j,row,col*/
         if (nox) xval = i+1;			/* xval is line number */
         else     xval = *(plotarr+i);		/* find max, min of x */
         yval = *(plotarr+i+k*lines);		/* find max, min of y */
       }
       else {
         if (nox) xval = i+1;			/* xval is line number */
         else     xval = *(plotarr+i*cols);	/* find max, min of x */
         yval = *(plotarr+i*cols+k);		/* find max, min of y */
       }
       if (xmax<xval) xmax = xval;
       if (xmin>xval) xmin = xval;
       if (ymax<yval) ymax = yval;
       if (ymin>yval) ymin = yval;
     }
    plotnod[pl].pxmax = xmax;
    plotnod[pl].pxmin = xmin;
    plotnod[pl].pymax = ymax;
    plotnod[pl].pymin = ymin;
    plotnod[pl].xrange = xmax - xmin;
    plotnod[pl].yrange = ymax - ymin;
    plotnod[pl].ppen = colortab[pl];
    plotnod[pl].pmod = GREC;
    plotnod[pl].automax = 1;
   }

   if (!nox) {
     xmax=ymax = -1e30;
     xmin=ymin = 1e30;
     for (j=0; j<totplots; j++) {
       pl = looklist(j);         	/* check plots actually used */

       xval = plotnod[pl].pxmax;	/* find overall maxes */
       yval = plotnod[pl].pymax;
       if (xmax<xval) xmax = xval;
       if (ymax<yval) ymax = yval;

       xval = plotnod[pl].pxmin;	/* find overall mins */
       yval = plotnod[pl].pymin;
       if (xmin>xval) xmin = xval;
       if (ymin>yval) ymin = yval;
     }
     if (!notinit(setxmax)) xmax = setxmax; /* set defaults from overrides */
     if (!notinit(setxmin)) xmin = setxmin;
     if (!notinit(setymax)) ymax = setymax;
     if (!notinit(setymin)) ymin = setymin;

     plotnod[0].pxmax = xmax; 
     plotnod[0].pxmin = xmin; 

     for (j=0; j<totplots; j++) {	/* set all ymin, ymax to overall vals */
       pl = looklist(j);         	/* check plots actually used */
       plotnod[pl].pymin = ymin; 
       plotnod[pl].pymax = ymax; 
     }
   }    /* if (!nox) */

   nonodes = 1;

   plotinit (totplots);		/* make and label axes of graph */

   if (nox) tcols = totplots;
   else     tcols = totplots+1;
   if (dashfl) {
     for (k=jstart; k<tcols; k++) {
       if (nox) {
         pl = k;
         j = looklist(pl);             /* possibly modify plot order */
       }
       else  {
         pl = k-1;  
         j = looklist(pl) + 1;         /* possibly modify plot order */
       }
       dash ((long)(pl & 7)); /* */
       for (i=0; i<lines; i++) {
         if (rowmode) {			/* if rowmode,swap i,j,row,col */
          if (nox) 
           mplot (*(plotarr+i+j*lines),(double)(i+1),totplots,pl,1);
	  else
           mplot (*(plotarr+i+j*lines),*(plotarr+i),totplots,pl,1);
        }
        else {
         if (nox)
           mplot (*(plotarr+i*cols+j),(double)(i+1),totplots,pl,1);
	 else
           mplot (*(plotarr+i*cols+j),*(plotarr+i*cols),totplots,pl,1);
      } /* else */
     }  /* for (i=1;;) */
    }  /* for (k=1;;) */
   }
   else {	/* not dashfl */
      for (i=0; i<lines; i++) {
       for (k=jstart; k<tcols; k++) {
        if (nox) {
          pl = k;
          j = looklist(pl);             /* possibly modify plot order */
        }
        else  {
          pl = k-1;  
          j = looklist(pl) + 1;         /* possibly modify plot order */
        }
         if (rowmode) {			/* if rowmode,swap i,j,row,col */
          if (nox) 
           mplot (*(plotarr+i+j*lines),(double)(i+1),totplots,pl,1);
	  else
           mplot (*(plotarr+i+j*lines),*(plotarr+i),totplots,pl,1);
        }
        else {
         if (nox)
           mplot (*(plotarr+i*cols+j),(double)(i+1),totplots,pl,1);
	 else
           mplot (*(plotarr+i*cols+j),*(plotarr+i*cols),totplots,pl,1);
      } /* else */
     }  /* for (k=1;;) */
    }  /* for (i=1;;) */
  }  /* if (dashfl) {} else {} */

 } /* if (cols) */
/* ncfprintf (stderr,"lines %d cols %d\n",lines,cols); */
}

