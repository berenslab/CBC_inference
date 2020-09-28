/* program plotmod */

/* Program to plot output of "nc" on screen */

#include <stdio.h>
#include "ncsub.h"
#include "ncelem.h"
#include "ncval.h"
#include "ncomp.h"
#include "ncplot.h"
#include "controlx.h"
#include "ncio.h"
#include "nc.h"

#include "nclist.h"

plotfr plotnod[PLOTNODSIZ]={0};	/* node numbers, zero sample voltage */

int numpl=0;			/* number of command-line plots */
extern int numlaxis;		/* number of left axes  */
extern int numraxis;		/* number of right axes  */

int setvid=0;			/* flag to set video mode on output */
int unsetvid=0;
int drseed=0;
int ncerror=0;			/* error flag used by warning for screen disp */
int tabmode=0;			/* print tabs on output */
int plotnum=-1;
int nonodes=0;			/* =1 -> no node numbers in video display */ 
int rowmode=0;			/* =1 -> numbers for plot in a row, not col */
int nofilfl=0;			/* =1 -> no filename on plot */
int finfl=0;			/* =1 -> final plot, no filename or labels */
int nolabels=0;			/* =1 -> no labels in display */
int setauto=0;			/* =1 -> set max and min for graph from data */
int charind=0;			/* index for plots with -c option */
int setnox=0;			/* =1 -> no x value, multiple y vals in input */
int yprecision=6;		/* precision of y value (in mplot()) */
int disp_ray=0;			/* make .pov file (not in plotmod) */
double charsiz=CHARSIZ;		/* size of label chars */

double ssetxmax = LARGENUM;	/* override for xmax */
double ssetxmin = LARGENUM;	/* override for xmin */

int setplsep = 0;		/* override plsep (separate plots) */

extern int listplots[PLOTNODSIZ];	/* plot order from command line */
extern int leftaxes[PLOTNODSIZ];	/* left axes set from command line */
extern int rightaxes[PLOTNODSIZ];	/* right axes set from command line */

char *infile=0;			/* the input file */
char *srcfile=0;		/* the source file that generated data */
char *stfile=0;			/* stimulus file */
char *progname=0;		/* name of program from shell. */
char *plotlabel=0;		/* label for plot in "vid" */

FILE *pictin;

#ifdef __cplusplus
extern "C" {
#endif
  #include "stdplt.h"
  double atof(const char *);
  double strtod(const char *, char **);
  FILE *freopen(const char *, const char *, FILE *);
  void exit (int);
#ifdef __cplusplus
}
#endif

void plotrun(void);
void initpl(int);
void getplot(FILE *instream, int numpl, int setauto);
void printphelp();
char *printversion(int ncversion);
FILE *openfilz (char *filnam, char *infilm, const char *fm);

#define FILENAMLEN 200

/*------------------------------------*/

int main(int argc, char **argv)
{
   char *cptr,*strp;
   FILE *temp;
   int i,pstart,pstop,ptemp;
   static double pcsize=0.018;
	 
 pictin = stdin;
 stdplt = stdout;
 for (i=0; i<PLOTNODSIZ; i++) {
   listplots[i] = 0;
   leftaxes[i]  = 0;
   rightaxes[i] = 0;
 }
 listplots[0] = 0;			/* reset plot order table */
 vidmode = 1;				/* turn on video mode by default */
 setxmax = LARGENUM;
 setxmin = LARGENUM;
 setymax = LARGENUM;
 setymin = LARGENUM;
 progname = argv[0];			/* program name */
 do					/* if there are any switches */
  {
   argc--; argv++;
   cptr = *argv;
   if (argc)
    {
     if (*cptr == '-')
      {
       cptr++;
       switch (*cptr)
        {
	  case 'a': setauto = !setauto; break;

	  case 'c': argv++; argc--;
 		    plotnod[charind].charfl = **argv; /* char for plotting */
 		    plotnod[charind].charmode = LINES;    
 		    if (**argv==' ') 		     /* reset if turned off */
     		       plotnod[charind].charmode = PLINES;
 		    plotnod[charind].csize = pcsize;    
		    if (++charind >= PLOTNODSIZ) charind = PLOTNODSIZ-1;
		    break;

	  case 'C': argv++; argc--;
 		    plotnod[charind].charfl = **argv; /* char for plotting */
 		    plotnod[charind].charmode = NOLINES;    
 		    if (**argv==' ') 		      /* reset if turned off */
     		       plotnod[charind].charmode = PLINES;
 		    plotnod[charind].csize = pcsize;    
		    if (++charind >= PLOTNODSIZ) charind = PLOTNODSIZ-1;
		    break;

	  case 'h': printphelp();		/* print help table */
		    break;

	  case 'l': argv++; argc--;		/* set output file label */
 		    plotlabel = *argv;
		    break;

	  case 'w': argv++; argc--;		/* set plot char size */
 		    pcsize = atof(*argv);	/* use this switch before -c */
		    break;

	  case 'W': argv++; argc--;		/* set label char size */
 		    charsiz = atof(*argv);
		    break;

	  case 'd': dashfl = !dashfl; break;

	  case 'S': setplsep = 1; break;

	  case 'f': nofilfl = !nofilfl; break;

	  case 'F': finfl = !finfl; break;

	  case 't': vidmode = !vidmode;
	  	    unsetvid = 1;
 		    break;
 
	  case 'x': tabmode = !tabmode; break;
	  
	  case 'n': nonodes = !nonodes; break;	/* no node numbers in display */

	  case 'r': rowmode = !rowmode; break;	/* numbers for plot in a row */

	  case 's': argv++; argc--;
		    yprecision = atof(*argv); /* precision of y values in printout */
		    if (yprecision<=0) yprecision = 6;
		    break;
	  
	  case 'y': setnox = !setnox; break;	/* no x values in input */

 	  case 'L': 				/* specify indiv.  L axes */
		argv++; argc--;
		pstart = (int)strtod(*argv,&strp);
		if (*strp=='-') {		/*  or specify range: a-b */
		  pstop = (int)strtod((strp+1),&strp);
		  if (pstop  > PLOTNODSIZ) pstop  = PLOTNODSIZ-1;
		  if (pstop  < 0) pstop  = 0;
		  if (pstart > PLOTNODSIZ) pstart = PLOTNODSIZ-1;
		  if (pstart  < 0) pstart  = 0;
		  if (pstop < pstart) {		/*  if they're backwards */
		    for (i=pstart; i>=pstop; i--) {
			leftaxes[numlaxis++] = i;  /* make range of plots */
		    }
		  }
		  else {
		    for (i=pstart; i<=pstop; i++) {
			leftaxes[numlaxis++] = i;  /* make range of plots */
		    }
		  }
		}
		else leftaxes[numlaxis++]  = pstart;
		break;

 	  case 'R': 				/* specify indiv.  R axes */
		argv++; argc--;
		pstart = (int)strtod(*argv,&strp);
		if (*strp=='-') {		/*  or specify range: a-b */
		  pstop = (int)strtod((strp+1),&strp);
		  if (pstop  > PLOTNODSIZ) pstop  = PLOTNODSIZ-1;
		  if (pstop  < 0) pstop  = 0;
		  if (pstart > PLOTNODSIZ) pstart = PLOTNODSIZ-1;
		  if (pstart  < 0) pstart  = 0;
		  if (pstop < pstart) {			/* if they're backwards*/
		    for (i=pstart; i>=pstop; i--) {
			rightaxes[numraxis++] = i;	/* make range of plots */
		    }
		  } else {
		    for (i=pstart; i<=pstop; i++) {
			rightaxes[numraxis++] = i;	/* make range of plots */
		    }
		  }
		}
		else rightaxes[numraxis++]  = pstart;
		break;

	  case 'B': allowrev = !allowrev; break;/* allow backwards plotting */

	  case 'E': 
		argv++; argc--;
		ssetxmax = atof(*argv);		/* override xmax */
		break;
 
	  case 'e': 
		argv++; argc--;
		ssetxmin = atof(*argv);		/* override xmin */
		break;
 
	  case 'M': 
		argv++; argc--;
		setymax = atof(*argv);		/* override ymax */
		break;
 
	  case 'm': 
		argv++; argc--;
		setymin = atof(*argv);		/* override ymin */
		break;
 
 	  case 'p': 				/* specify indiv. plots */
		argv++; argc--;
		pstart = (int)strtod(*argv,&strp);
		if (*strp=='-') {		/*  or specify range: a-b */
		  pstop = (int)strtod((strp+1),&strp);
		  if (pstop  > PLOTNODSIZ) pstop  = PLOTNODSIZ-1;
		  if (pstop  < 0) pstop  = 0;
		  if (pstart > PLOTNODSIZ) pstart = PLOTNODSIZ-1;
		  if (pstart  < 0) pstart  = 0;
		  if (pstop < pstart) {		/*  if they're backwards */
		     for (i=pstart; i>=pstop; i--) {
			listplots[numpl++] = i;	/* make range of plots */
		     }
		  }
		  else {
		    for (i=pstart; i<=pstop; i++) {
			listplots[numpl++] = i;	/* make range of plots */
		    }
		 }
		}
		else listplots[numpl++]  = pstart;
		break;

 	/*     	  case 'e': 
		argv++; argc--;
		stop = atof(*argv);
		stopfg = T;
		break;
	*/

	  default:
		ncfprintf (stderr,"plotmod: unknown switch '%s'\n",*argv);
		exit(1);

        }  /* switch */
      }	   /* if (*cptr) */
     else
      {
          char infil[FILENAMLEN]={0};
       if((temp=openfilz(cptr,infil,"r"))==NULL)
         {
           ncfprintf(stderr,"plotmod: cannot open %s\n",infil);
	   fflush (stderr);
           continue;
         }
       else  {
	   pictin = temp;
	   infile = infil;
       };
       plotrun();
       if (argc <= 1) break;
     }	/* else */

    } /* if (argc) */

   else plotrun(); 
  }
 while (argc > 0);
}

/*------------------------------------*/

void printphelp()

{
  ncfprintf (stdout,"## %s version #%-4.3g\n",progname,4.04);
  ncfprintf (stdout,"plotmod -a    auto max min mode.\n");
  ncfprintf (stdout,"        -B    allow backwards plotting.\n");
  ncfprintf (stdout,"        -c x  plot with char 'x', connected with lines.\n"); 
  ncfprintf (stdout,"        -C x  plot with char 'x', no lines.\n"); 
  ncfprintf (stdout,"        -d    set dash mode.\n");
  ncfprintf (stdout,"        -E n  override Xmax.\n"); 
  ncfprintf (stdout,"        -e n  override Xmin.\n"); 
  ncfprintf (stdout,"        -f    no file label.\n");
  ncfprintf (stdout,"        -F    no labels in plot.\n");
  ncfprintf (stdout,"        -M n  override Ymax.\n"); 
  ncfprintf (stdout,"        -m n  override Ymin.\n"); 
  ncfprintf (stdout,"        -n    no node numbers on plot.\n");
  ncfprintf (stdout,"        -p n  specify plots in different order.\n");
  ncfprintf (stdout,"           n-m specify range of plots.\n");
  ncfprintf (stdout,"        -r    row mode (numbers in rows instead of columns).\n");
  ncfprintf (stdout,"        -L n  add Y axes on left side of graph. (Like -p above)\n");
  ncfprintf (stdout,"        -R n  add Y axes on right side of graph. (Like -p above)\n");
  ncfprintf (stdout,"        -S    make separate plot for each trace.\n"); 
  ncfprintf (stdout,"        -t    text mode (prints out same numbers as input)\n");
  ncfprintf (stdout,"        -w n  set plot char size.\n"); 
  ncfprintf (stdout,"        -W n  set label char size.\n"); 
  ncfprintf (stdout,"        -y    no X values in input file.\n");
  exit(0);
}

/*------------------------------------*/

void plotrun(void)
{
  initpl(charind);
  getplot(pictin,numpl,setauto);
  listplots[0] = 0;
  numpl = 0;
}

/*------------------------------------*/

double ncabs(double arg)
{
  return ((arg<0) ? -(arg) : arg);
}

/*------------------------------------*/

void rprintf (char*, double, double, double, double, int) {}

