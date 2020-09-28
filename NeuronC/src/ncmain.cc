/* segment ncmain in program nc */

/* Main routines for interpreter */

#ifdef __cplusplus
extern "C" {
#endif

#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>

#ifdef CPML
#include <cpml.h>
#endif
#include <math.h>
#include <setjmp.h>
#include <sys/times.h>
#include <sys/types.h>

void doexit(int code);

#ifdef __cplusplus
}
#endif


#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;

#include "nc.h"
#include "ncval.h"
#include "y.tab.h"
#include "ncmain.h"
#include "ncio.h"
#include "ncsetvar.h"

jmp_buf begin_ncmain;

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

extern char *progname;
extern char *runfile;		/* interactive file name */
extern int makestim;		/* = 1 -> says we're running "stim", not "nc" */
extern int prcomps;		/* = 1 -> print compartments when displaying them */
extern int vidmode;		/* video mode set (defined in control.h) */
extern int silent;		/* =1 -> no display of info about simulation */
extern int info;		/* info level for printouts */
int setplsep=0;			/* separate plots on output (def in control.h)*/

double ssetxmax = LARGENUM;	/* override for "endexp" in plotinit() */
double ssetxmin = LARGENUM;	/* override for "time" in plotinit() */


int	pipefl=0;		/* =1 -> pipe to vid, wait on exit */
int	nofilfl=0;		/* =1 -> no filename in plot */
int	finfl=0;		/* =1 -> no labels in plot */
int	nolabels=0;		/* =1 -> no labels in display */
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
int	yprecision = 6; 	/* precision of y values in printout (in mplot) */
int	cppfl = 0;		/* =1 -> pre-process with cpp */
int	disp_ray = 0;		/* =1 -> render 3D with ray-tracing */
int     interp_only = 0;       	/* =1 -> run interpreter, no simulator */
int     interp = 0;       	/* =1 -> interpreter is running, need to run varcopy */
int	comndfl = 0;		/* =1 -> look for command args after file name */
int	prcmdlin = 0;		/* =1 -> print command line to stdout */

extern int runyet;		/* resets time, in ncsub.c */
double charsiz=CHARSIZ;		/* char size for plot labels */
extern double ncversion;	/* nc version number */
double linewidth=LINEWIDTH;  /* line width for ncdisp.c */

stype varval[VVSIZE] = {0};	/* temp values for setting variables */
int varset = 0;			/* index for setting variables */

const char *srcfile = 0;	/* source file name */
const char *einfile = NULL;   	/* input file name for error reports */
const char *infile = "stdin";	/* input file name, used by "plotmod" */
const char *plotlabel = NULL;	/* label for output in "vid" */

int *istat;			/* interrupt status upon entering program */
extern int set_tty;		/* defined in ncm.cc, =1 -> keyboard input */

extern FILE *stimout;
extern FILE *stimin;
FILE *compout = stdout;
istream *fin;

extern std::ostream *nc_stdout;
extern std::ostream *nc_stderr;

void change_stdout(char *fname);
void change_stderr(char *fname);

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
double errcheck(double d, const char *s);
void simwait(double nsecs);
void set_argv_array(int argc, char **argv);

void ncmdlin(int argc, char **argv);

char **dargv; 			/* for print_cmdlin and set_argv_array() */
int    dargc;


#define NCVERSION 6.314

#ifdef XSTIM
  #define STIMCODE 1
#else
  #define STIMCODE 0
#endif

/*------------------------------------*/

int ncmain(int argc, char **argv)
{
 interp = 1;			/* running the interpreter here, must run varcopy */
 				/*  but can run simulator from external program */
 if (setjmp(begin_ncmain) > 0) return 0;
 fin = &cin;
 ncversion = NCVERSION;
 progname = basename(argv[0]);
 makestim = STIMCODE;
 prcomps = 0;
 runfile = (char *)"ncfile";
 einfile = (char *)infile;
 if (argc==1)			/* if user needs help */
   ncrun(1);
 else {
  ncmdlin(argc,argv); 		/* get command line switches and files */
  if (fin==NULL) doexit(3);
  if (!srcfile) ncrun(0);	/* arguments but no filename */
  else if (comndfl) ncrun(0);	/* arguments after filename */
 }
}

/*------------------------------------*/

void change_stdout(char *fname)
{
  fstream *fout = new fstream();
  fout->open(fname, ios::out);
  nc_stdout = fout;
}

void change_stderr(char *fname)
{
  fstream *ferr = new fstream();
  ferr->open(fname, ios::out);
  nc_stderr = ferr;
}

/*------------------------------------*/

void print_cmdlin(int argc, char **argv)

{
  printf ("# ");
  for (argc; argc>0; argc--,argv++) 
    printf ("%s ",*argv);
  printf ("\n");
}

/*------------------------------------*/
   
void ncmdlin(int argc, char **argv)

{
   char *cptr;

 ncversion = NCVERSION;
 progname = basename(argv[0]);
 dargc = argc;				/* save the command line for printout */ 
 dargv = argv;
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
        case 'c': comndfl = 1;			 /* "-c" after filename in first line of batch file */
                  fin = new ifstream(*(argv+1)); /* set interp to read from rest of batch file */
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

        case 'f': nofilfl = !nofilfl;      	 /* =1 -> no filename in plot */
                  break;

        case 'F': finfl = !finfl;      	 	/* =1 -> no filename or labels in plot */
                  break;

        case 'K': { 				/* print all symbols */ 
		     datum d;
		     init();	/* init variables, symbol table, in "init.c" */
		     d = printsym();
		     doexit(0);
		  }
                  break;

        case 'n': nonodes = !nonodes;         /* =1 -> no node numbers in plot */
                  break;

        case 'N': nolabels = !nolabels;       /* =1 -> no labels in morphology display */
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

        case 's': argc--; argv++;
		  yprecision = setint(atof(*argv));    /* precision of y values in printout */
		  if (yprecision<=0) yprecision = 6;
                  break;

        case '-': cptr++;                      /* set variable */
                  varval[varset].name = cptr;  /* remember variable name */
                  argc--; argv++;
                  if (argc) {
		   if (!isdigit(**argv) && ! (**argv == '-' && isdigit(*((*argv)+1)))
		   			&& ! (**argv == '.' && isdigit(*((*argv)+1)))
		   			&& ! (**argv == '-' && *((*argv)+1)=='.' &&
								isdigit(*((*argv)+2)))
			) {
                     varval[varset].cptr = *argv;	/*remember str*/
                     varval[varset].type = VSTRING;
                   }
                   else {
                     varval[varset].val  = atof(*argv); /* remember value */
                     varval[varset].type = VNUMBER;
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

	case 'D': prcmdlin = !prcmdlin;		/* print command line on stdout */
		  break;

	case 'R': disp_ray = !disp_ray;		/* render 3D ray-tracing */
		  break;

	case 'I': interp_only = 1;		/* no simulator */
		  break;

	case '1': argc--; argv++;               /* redirect stdout */
		  change_stdout(*argv);
		  break;

	case '2': argc--; argv++;               /* redirect stderr */
		  change_stderr(*argv);
		  break;

	default:
		ncfprintf (stderr,"nc: unknown switch '%s'\n",*argv);
		doexit(1);

        }  /* switch */
      }	   /* if (*cptr=) */
     else
      {
       if (!comndfl) {
       if((fin = new ifstream(cptr))->fail())
         {
           ncfprintf(stderr, "%s: can't open %s\n", progname, cptr);
	   fin = NULL;
           continue;
         }
       srcfile = cptr;
       infile = srcfile;
       if (!comndfl) ncrun(0);
       if (argc <= 1) break;
      }			/* if (!comndfl) */
     }			/* else (*cptr != '-') */
    }			/* if (argc) */
  }			/* do */
 while (argc > 0);
 if (prcmdlin) print_cmdlin(dargc,dargv);
}


/*------------------------------------*/

void doexit(int code)

{
  exit(code);
 // longjmp (begin_ncmain,0);		/* return from ncmain() */
}

/*------------------------------------*/

void printhelp()

{
  ncfprintf (stdout,"## %s: nc version %s\n",
			progname,print_version(ncversion));
  ncfprintf (stdout,"nc -v      video mode to stdout (makes graphics).\n");
  ncfprintf (stdout,
	"   --var n  set variable from command line.\n");
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
  ncfprintf (stdout,"   -M n -m n override ymax, ymin on plots\n");
  ncfprintf (stdout,"   -l n    set 'lamcrit' variable (0=no condensation).\n");
  ncfprintf (stdout,"   -f      no file name on graph.\n");
  ncfprintf (stdout,"   -F      no labels on graph.\n");
  ncfprintf (stdout,"   -I      run as interpreter only.\n");
  ncfprintf (stdout,"   -K      print out all predefined symbols.\n");
  ncfprintf (stdout,"   -n      no node numbers on graph.\n");
  ncfprintf (stdout,"   -N      no labels in morphology display.\n");
  ncfprintf (stdout,"   -r n    random number seed. Negative = different each time.\n");
  ncfprintf (stdout,"   -s n    set precision of output numbers. def=6.\n");
  ncfprintf (stdout,"   -C      run input file through CPP preprocessor .\n"); 
  ncfprintf (stdout,"   -R      output symbolic display for \"povray\" ray-tracer.\n"); 
  ncfprintf (stdout,"   -L n    set line width for display statement.\n");
  ncfprintf (stdout,"   -w n    set vid window size for \"ncv\" or \"ndv\".\n");
  ncfprintf (stdout,"   -W n    tic label char with in terms of screen size.\n");
  ncfprintf (stdout,"   -y n    debug level.\n");
  ncfprintf (stdout,"   -z n    debug category.\n");
  ncfprintf (stdout,"   -1 fn   redirect stdout to file \"fn\"\n");
  ncfprintf (stdout,"   -2 fn   redirect stderr to file \"fn\"\n");

  doexit(0);
}

/*------------------------------------*/

void ncrun(int hlp)
{
/* if(fin== &cin) system ("stty erase \n"); 		/* */

  /* if (cppfl) {
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
  } /* */

  timepnt = &timebuf;
  startclk = times(timepnt);
  startu = timepnt->tms_utime;
  runyet = 0;		/* resets time */
  if (srcfile) einfile = srcfile;
  ncleanup();		/* erase lists  (in "ncmak.c")  */
  initpl(0);		/* init plot stuff, in "ncplot.c" */
  init();			/* init variables, symbol table, in "init.c" */
  if (hlp && info!=0) 
     ncfprintf (stdout,"## %s version %s\n",progname,print_version(ncversion));
  if (comndfl) set_argv_array(dargc, dargv);
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

/*------------------------------------*/

void bexit(void)

/* things to do before exit */

{
  if (!vidmode && !disp_ray) prtime(stdout);  /* print run times */
  if (pipefl) sleep(0); 		      /* if pipe to vid, wait */
}

/*------------------------------------*/

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

/*------------------------------------*/

double elap_time(void)

/* print compute time since start */

{

  totmin = (times(timepnt) - startclk) / SECPERMIN / timecalib;
  runmin = (timepnt->tms_utime - startu) / SECPERMIN / timecalib;
  return(runmin);
}

/*------------------------------------*/

datum elap_timee(void)

{
   datum d;

  d.val = elap_time();
  d.vtype = NUMBER;
  return (d);			/* return elapsed run time */
}

/*------------------------------------*/

double system_speed(void)

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
   static double retval;
   static int runyet = 0;

  if (!runyet) { 
    runyet = 1;
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
    MHz = 3000;
    runtimetest=0.00256773;
    runtimetest=0.00438;
    runtimetest=0.00340;
    retval = runtimetest / runtime * MHz;
  }
  return retval;
}

/*------------------------------------*/

void simwait (double nsecs)

{
   int i, j;
   long unsigned int endw;

   // fprintf (stderr,"simwait speed %g\n",system_speed());

   if (nsecs<0) nsecs = 0;
   
   endw = (long unsigned int)(nsecs  * 100);
   for (i=0; i<endw; i++)
      usleep(10000);
//  endw = (long unsigned int)(nsecs * 4.8 * system_speed());
//  for (i=0; i<endw; i++)
//    for (j=0; j<100000; j++) ;
//
}
      

/*------------------------------------*/
datum system_speedd(void)

{
   datum d;

  d.val = system_speed();
  d.vtype = NUMBER;
  return d;			/* return elapsed run time */
}

/*------------------------------------*/

datum simwait(datum &d)

{
   simwait(d.val);
   return d;
}

/*------------------------------------*/

datum xgetpid(void)
{
        datum d;

        d.val=errcheck(getpid(), "getpid");
        d.vtype = NUMBER;
        return d;
}

