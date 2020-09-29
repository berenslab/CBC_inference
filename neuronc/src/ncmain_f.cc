/* segment ncmain in program ncmod */

/* Main routines for interpreter */

extern "C" {

#include <ctype.h>
#include <unistd.h>
#include <stdio.h>

#ifdef CPML
#include <cpml.h>
#endif

#include <math.h>
#include <signal.h>
#include <setjmp.h>
#include <sys/times.h>
#include <sys/types.h>

void doexit(int code);

}

#include "nc.h"
#include "ncval.h"
#include "y.tab.h"
#include "ncmain.h"
#include "ncio.h"

jmp_buf begin;

#define SECPERMIN 60.0		/* number of seconds per minute */

/* #define TIMECALIB 60.0	/* times() system call gives 1/60  sec */
/* #define TIMECALIB 100.0	/* times() system call gives 1/100 sec */
/* #define TIMECALIB 1024.0	/* times() system call gives 1/1024 sec */

#ifdef DEC_ALPHA
#define TIMECALIB 1024.0	/* times() system call gives .9765625 msec */
#else
#define TIMECALIB 100.0		/* times() system call gives 1/100 sec */
#endif

struct tms timebuf;
struct tms *timepnt;
double runmin=0.0;
double timecalib=TIMECALIB;
double startclk,startu,totmin;

extern char	*progname;
extern char 	*runfile;	/* interactive file name */
extern int makestim;		/* = 1 -> says we're running "stim", not "nc" */
extern int prcomps;
extern int vidmode;		/* video mode set (defined in control.h) */
extern int silent;		/* =1 -> no display of info about simulation */
extern int info;		/* info level for printouts */
int setplsep=0;			/* separate plots on output (def in control.h)*/

double ssetxmax = LARGENUM;	/* override for "endexp" in plotinit() */
double ssetxmin = LARGENUM;	/* override for "xtime" in plotinit() */


int	pipefl=0;		/* =1 -> pipe to vid, wait on exit */
int	nofilfl=0;		/* =1 -> no filename in plot */
int	finfl=0;		/* =1 -> no labels in plot */
int     setvid = 0;     	/* set vid mode on output */
int     unsetvid = 0;   	/* unset vid mode on output */
int     setdebug = 0;   	/* set debug level from command line */
int     setdebgz = 0;   	/* set debug mode from command line */
int     setdisp = -1;   	/* set disp mode from command line */
int     setprmap = -1; 		/* set print mode from command line */
double  setlamcrit = -1;	/* set lambda crit from command line */
double	vidsiz = 1.0;		/* size of vid window for "ncv" or "ndv" */
int     drseed = 0;    		/* random seed set from command line */
int     tabmode = 0;   		/* set tabs on output from command line */
int	nonodes = 0; 		/* =1 -> no nodes in video display */
int	rowmode = 0; 		/* =1 -> numbers for plot in row */
int	setnox = 0; 		/* =1 -> no x value, multiple y vals on input*/
int	cppfl = 0;		/* =1 -> pre-process with cpp */
int	disp_ray = 0;		/* =1 -> render 3D with ray-tracing */
int     interp_only = 0;       	/* =1 -> run interpreter, no simulator */
int	comndfl = 0;		/* =1 -> look for command args after file name */

extern int runyet;		/* resets time, in ncsub.c */
double charsiz=CHARSIZ;		/* char size for plot labels */
extern double ncversion;	/* nc version number */
double linewidth=LINEWIDTH;  /* line width for ncdisp.c */

double varval[VVSIZE] = {0};	/* temp values for setting variables */
char *varpnt[VVSIZE] = {0};	/* temp pointers for setting variables */
int vartyp[VVSIZE] = {0};	/* temp type for setting variables */
int varset = 0;			/* index for setting variables */

char    *srcfile = 0;		/* source file name */
char    *infile = "stdin";	/* input file name, used by "plotmod" */
char    *einfile;       	/* input file name for error reports */

int istat;	/* interrupt status upon entering program */
extern int set_tty;		/* defined in ncm.cc, =1 -> keyboard input */

extern FILE *stimout;
extern FILE *stimin;
extern FILE *compout;
FILE *fin;

void ncrun(int hlp);
void ncleanup();
void initpl(int charind);
void init();
void run_interp();	
void prtime(FILE *stimout);
void sortstimfile();
void bexit();
void printhelp();
char *print_version(double version);
datum printsym();
int setint (double val);

/*------------------------------------*/

int ncmain(int argc, char **argv)
{
   char *cptr;

 if (setjmp(begin) > 0) return 0;
 fin = stdin;
 runfile = "ncfile";
 einfile = infile;
 ncversion = 5.758;
 progname = argv[0];
 if (argc==1)			/* if user needs help */
   ncrun(1);
 else {
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
        case 'c': comndfl = !comndfl;		/* comnd flags after filename */
                  break;

        case 'd': argc--; argv++;    		 /* display comps and conns */
		  setdisp = (int) atof(*argv);
                  break;

        case 'e': argc--; argv++;    		 /* xmin override */
		  ssetxmin = atof(*argv);
                  break;

        case 'E': argc--; argv++;    		 /* xmax override */
		  ssetxmax = atof(*argv);
                  break;

        case 'h': printhelp();			/* print help table */
                  break;

        case 'l': argc--; argv++;    		/* set lambda crit. */
		  setlamcrit = atof(*argv);
                  break;

        case 'f': nofilfl = !nofilfl;      	 /* =1 -> no filename in disp */
                  break;

        case 'F': finfl = !finfl;      	 	/* =1 -> no labels in disp */
                  break;

        case 'K': { 				/* print all symbols */ 
		     datum d;
		     init();	/* init variables, symbol table, in "init.c" */
		     d = printsym();
		     doexit(0);
		  }
                  break;

        case 'n': nonodes = !nonodes;       /* =1 -> no node numbers in disp */
                  break;

        case 'p': argc--; argv++;
		  setprmap = (int)atof(*argv);      /* print comps and conns */
                  break;

        case 'q': silent = 1;      		/* don't print info statements */
                  break;

        case 'r': argc--; argv++;
		  drseed = setint(atof(*argv));    /* set random seed */
		  if (drseed ==0) drseed = 1;
                  break;

        case 's': argc--;argv++;		/* set variable */ 
		  varpnt[varset] = *argv;       /* remember variable name */
		  argc--; argv++;
                  if (argc) {
		   if (isalpha(**argv)) {
		     varval[varset] = (double)(long int)*argv; /*remember str*/
		     vartyp[varset] = STRING;
		   }
		   else {
		     varval[varset] = atof(*argv); /* remember value */
		     vartyp[varset] = NUMBER;
		   }
		  }
		  varset++;
		  if (varset >= VVSIZE) varset = VVSIZE-1;
                  break;

        case '-': cptr++;                      /* set variable */
                  varpnt[varset] = cptr;       /* remember variable name */
                  argc--; argv++;
                  if (argc) {
		   if (isalpha(**argv)) {
                     varval[varset] = (double)(long int)*argv; /*remember str*/
                     vartyp[varset] = STRING;
                   }
                   else {
                     varval[varset] = atof(*argv); /* remember value */
                     vartyp[varset] = NUMBER;
                   }
		  }
                  varset++;
                  if (varset >= VVSIZE) varset = VVSIZE-1;
                  break;

        case 'S': setplsep = 1;		/* separate plots on output */
                  break;

        case 't': unsetvid = 1;         /* unset vid mode on output */
                  setvid = 0;
                  break;

        case 'v': setvid = 1;           /* set vid mode on output */
                  unsetvid = 0;
                  break;

        case 'V': { 
		     ncfprintf (stdout,"## nc version %s\n", 
					print_version(ncversion));
		     doexit(0);
		  }
                  break;

        case 'L': argc--; argv++;
		  linewidth = atof(*argv);      /* set line width for ncdisp */
                  break;

        case 'w': argc--; argv++;
		  vidsiz = atof(*argv);		/* vid size for "ncv" */
                  break;

        case 'W': argc--; argv++;
		  charsiz = atof(*argv);	/* char size for plot labels */
                  break;

        case 'x': tabmode = !tabmode;		/* print tabs on output */
                  break;

        case 'y': argc--; argv++;    		/* set debug level */
		  setdebug = (int)atof(*argv);
                  break;

        case 'z': argc--; argv++;    		/* set debug mode */
		  setdebgz = (int)atof(*argv);
                  break;

	case 'C': cppfl = !cppfl;		/* pre-process intput file */
		  break;

	case 'R': disp_ray = !disp_ray;		/* render 3D ray-tracing */
		  break;

	case 'I': interp_only = !interp_only;	/* no simulator */
		  break;

	default:
		ncfprintf (stderr,"nc: unknown switch '%s'\n",*argv);
		doexit(1);

        }  /* switch */
      }	   /* if (*cptr=) */
     else
      {
       if((fin=fopen(cptr,"r"))==NULL)
         {
           ncfprintf(stderr, "%s: can't open %s\n", progname, cptr);
           continue;
         }
       srcfile = cptr;
       infile = srcfile;
       if (!comndfl) ncrun(0);
       if (argc <= 1) break;
      }			/* else (*cptr != '-') */
    }			/* if (argc) */
  }			/* do */
 while (argc > 0);
 if (fin==NULL) doexit(3);
 if (!srcfile) ncrun(0);	/* arguments but no filename */
 else if (comndfl) ncrun(0);	/* arguments after filename */
 }
}

/*------------------------------------*/

void doexit(int code)

{
 longjmp (begin,0);		/* return from ncmain() */
}

/*------------------------------------*/

void printhelp()

{
  ncfprintf (stdout,"## %s version %s\n",
			progname,print_version(ncversion));
  ncfprintf (stdout,"nc -v      video mode to stdout (makes graphics).\n");
  ncfprintf (stdout,
	"   --var n  set variable from command line. (same as -s var n).\n");
  ncfprintf (stdout,
	"   -c      run from inside script with first line = #! nc -c\n");
  ncfprintf (stdout,"   -d 1    display 'neural elements' with 'display' statement.\n");
  ncfprintf (stdout,"   -d 2            'compartments'   \n");
  ncfprintf (stdout,"   -d 4            'connections'\n");
  ncfprintf (stdout,"   -d 8            'nodes'\n");
  ncfprintf (stdout,"   -d 16           'stimulus'\n");
  ncfprintf (stdout,"   -d 32           'movie (vcolor,cacolor)'\n");
  ncfprintf (stdout,"   -p 1    print out compartments as conductances.\n");
  ncfprintf (stdout,"   -p 2    print out compartments as spheres, chan densities.\n");
  ncfprintf (stdout,"   -q      quiet, don't print extra info (eqiv to \"--info 0\").\n");
  ncfprintf (stdout,"   -E n -e n override xmax, xmin on plots\n");
  ncfprintf (stdout,"   -l n    set 'lamcrit' variable (0=no condensation).\n");
  ncfprintf (stdout,"   -f      no file name on graph.\n");
  ncfprintf (stdout,"   -F      no labels on graph.\n");
  ncfprintf (stdout,"   -K      print out all predefined symbols.\n");
  ncfprintf (stdout,"   -n      no node numbers on graph.\n");
  ncfprintf (stdout,"   -r n    random number seed. Negative = different each time.\n");
  ncfprintf (stdout,"   -C      run input file through CPP preprocessor .\n"); 
  ncfprintf (stdout,"   -R      output symbolic display for \"povray\" ray-tracer.\n"); 
  ncfprintf (stdout,"   -L n    set line width for display statement.\n");
  ncfprintf (stdout,"   -w n    set vid window size for \"ncv\" or \"ndv\".\n");
  ncfprintf (stdout,"   -W n    tic label char with in terms of screen size.\n");
  ncfprintf (stdout,"   -y n    debug level.\n");
  ncfprintf (stdout,"   -z n    debug category.\n");
  doexit(0);
}

/*------------------------------------*/

void ncrun(int hlp)
{
/* if(fin==stdin) system ("stty erase \n"); 		/* */

  if (cppfl) {
          FILE *ftemp;
	if (srcfile) {
	    if ((ftemp=freopen (srcfile,"r",stdin))==NULL) {
              ncfprintf (stderr,"Error: nc: cpp can't open input file %s\n",srcfile);
	    }
	}
        if ((ftemp=popen("cpp -P","r"))==NULL) {
            ncfprintf (stderr,"Error: nc: cpp can't open input file\n");
        }
        else fin = ftemp;
  }
  timepnt = &timebuf;
  startclk = times(timepnt);
  startu = timepnt->tms_utime;
  if (hlp && !silent) 
     ncfprintf (stdout,"## %s version %s\n",progname,print_version(ncversion));
  runyet = 0;		/* resets time */
  if (srcfile) einfile = srcfile;
  ncleanup();		/* erase lists  (in "ncmak.c")  */
  initpl(0);		/* init plot stuff, in "ncplot.c" */
  init();			/* init variables, symbol table, in "init.c" */
  run_interp();			/* in "ncm.c"    */
  if (makestim) {
     ncfprintf(stimout,"## %s version %s\n",progname,print_version(ncversion));
     prtime(stimout);	       /* print out time */
     fclose (stimout);
     if (stimout!=stdout) sortstimfile();
  }
  if (prcomps) {
     ncfprintf(compout,"## %s version %s\n",progname,print_version(ncversion));
     prtime(compout);	       /* print out time */
     fclose (compout);
  }
  else bexit();
  if (set_tty) system ("stty sane");  
}

void bexit(void)

/* things to do before exit */

{
  if (!vidmode && !disp_ray) prtime(stdout);  /* print run times */
  if (pipefl) sleep(0); 		      /* if pipe to vid, wait */
}


void prtime(FILE *fout)

/* print user and total time from start */

{
  totmin = (times(timepnt) - startclk) / SECPERMIN / timecalib;
  runmin = (timepnt->tms_utime - startu) / SECPERMIN / timecalib;
  if (info >= 1) {
    ncfprintf (fout,"## time elapsed: %s %.4g  total %.4g minutes.\n",
					progname,runmin,totmin);
    fflush (fout);
  }
}

datum elap_time(void)

/* print compute time since start */

{
   datum d;

  totmin = (times(timepnt) - startclk) / SECPERMIN / timecalib;
  runmin = (timepnt->tms_utime - startu) / SECPERMIN / timecalib;
  d.val = runmin;
  d.vtype = NUMBER;
  return (d);			/* return elapsed run time */
}

datum system_speed(void)

/* Print CPU speed from simple benchmark. */
/*  to calibrate: */
/*  1) Compile this file with -O3 ... as in "makefile". */
/*  2) Run script that contains "print system_speed();" */
/*  3) Put elapsed time from nc printout into "runtimetest" below. */
/*  4) Put known CPU speed into MHz below. */
/*  5) Re-compile */
/*  6) Re-run script -- it should print out correct CPU speed */
/*  7) Run script on different machine -- it should print correct CPU speed */

{
   int i;
   double a,b,c,x,y,z;
   double totmin2,runmin2;
   double runtime, runtimetest,MHz;
   datum d;
  
  totmin = (times(timepnt) - startclk) / SECPERMIN / timecalib;
  runmin = (timepnt->tms_utime - startu) / SECPERMIN / timecalib;
  x = 2.134e20;
  y = 1.000000001;
  for (i=0; i<4000000; i++) {
    a = x + y;
    b = x - y;
    c = a + b + a / b;
    z = x * y * c * sqrt(x);
    x = z;
  }; 
  totmin2 = (times(timepnt) - startclk) / SECPERMIN / timecalib;
  runmin2 = (timepnt->tms_utime - startu) / SECPERMIN / timecalib;
  runtime = runmin2 - runmin;
  if (runtime==0) runtime = 1e-6;
  MHz = 1800;
  runtimetest=0.00256773;
  d.val = runtimetest / runtime * MHz;
  d.vtype = NUMBER;
  return (d);			/* return elapsed run time */
}

