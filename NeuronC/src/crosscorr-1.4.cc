/* Program crosscorr */

/* file to compute crosscorrelation between stimulus and response */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <arpa/inet.h>

#ifdef __cplusplus
}
#endif
#include <stdlib.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/times.h> 

/* #ifndef _SYS_TIMES_H_
 long times(struct tms *);
#endif				/* */

void exit(int code);

#ifdef __cplusplus
}
#endif

#define islower(c)       ('a' <= (c) && (c) <= 'z')
#define isupper(c)       ('A' <= (c) && (c) <= 'Z')
#define isalpha(c)       (islower(c) || isupper(c))

#define abs(x)           ((x) < 0 ? -(x) : (x))

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

char *progname;

char *stimfile = 0;		/* stimulus file */
char *respfile = 0;		/* response file */


#define STRING 2		/* used in stype below */
#define NUMBER 3
#define INT    4
#define DOUBLE 5
#define VVSIZE 30
#define TIME 0


struct stype {
	char *name;
	int type;
	union {
	    double val;
	    int    *iptr;
	    char   *cptr;
	    char   **sptr;
	    double *dptr;
	};
};

stype varval[VVSIZE] = {0};	/* temp values for setting variables */
int varset = 0;			/* index for setting variables */

int version[3];

#define TIMCOL 0	/* column number for time in data file */
#define   VCOL 1	/* column number for voltage in data file */

#define ROWADDR 40
#define COLADDR 41
#define MATBUFSIZ 48
#define FRAMERATE 119.78  /*59.82*/     /* 119.78, Narenn's Rig */
#define NBOX 15
#define DURATION 5000
#define TOTFRAMES int(FRAMERATE*DURATION)
#define CCFRAMES  60

short int stimlist[TOTFRAMES][NBOX][NBOX];

double ccout[CCFRAMES][NBOX][NBOX] = {0};

int nbox;
int ccframes;
double mean_bg;
double rest_v;
double framerate;
double timestep, t1, t2;
int info=0;
int swap=0;

/* command line variables */

FILE *sfile;
FILE *rfile;

void run(int hlp);
void prtime(FILE *stimout);
void printhelp();
char *print_version(int *version);

/*------------------------------------*/

main(int argc, char **argv)
{
   char *cptr;

 rfile = stdin;
 version[0] = 1;
 version[1] = 0;
 version[2] = 4;
 progname = argv[0];
 if (argc==1)			/* if user needs help */
   run(1);
 else
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
        case 'h': printhelp();			/* print help table */
                  break;

        case 's': argc--; argv++;	        /* set variable */
		  varval[varset].name = *argv;  /* remember variable name */
		  argc--; argv++;
		  if (isalpha(**argv) || ispunct(**argv)) {
		     varval[varset].cptr = *argv; /*remember str*/
		     varval[varset].type = STRING;
		  }
		  else {
		     varval[varset].val = atof(*argv); /* remember value */
		     varval[varset].type = NUMBER;
		  }
		  varset++;
		  if (varset >= VVSIZE) varset = VVSIZE-1;
                  break;

        case '-': cptr++;	        	/* set variable */
		  varval[varset].name = cptr;   /* remember variable name */
		  argc--; argv++;
		  if (isalpha(**argv) || ispunct(**argv)) {
		     varval[varset].cptr = *argv; /*remember str*/
		     varval[varset].type = STRING;
		  }
		  else {
		     varval[varset].val = atof(*argv); /* remember value */
		     varval[varset].type = NUMBER;
		  }
		  varset++;
		  if (varset >= VVSIZE) varset = VVSIZE-1;
                  break;

	default:
		fprintf (stderr,"avg: unknown switch '%s'\n",*argv);
		exit(1);

        }  /* switch */
      }	   /* if (*cptr=) */
     else
      {
       respfile=cptr;
       argv++; argc--;
       cptr = *argv;
       stimfile = cptr;

   /*       if((fin=fopen(respfile=cptr,"r"))==NULL)
         {
           fprintf(stderr, "%s: can't open %s\n", progname, cptr);
           continue;
         }
   */
       run(0);
       if (argc <= 1) break;
      }			/* else (*cptr != '-') */
     if (argc<=1) run(1);		/* there were arguments but no files */
    }			/* if (argc) */
  }			/* do */
 while (argc > 0);
}

/*------------------------------------*/

void printhelp()

{
  printf ("## %s %s\n",progname,print_version(version));
  printf ("\n");
  printf ("   Usage: %s options filename\n",progname);
  printf ("\n");
  printf ("   --stimfile   s  stimulus file (\"binlist.mat\")\n");
  printf ("   --respfile   s  response file (\"??\")\n");
  printf ("   --nbox       n  size of checkerboard (15, set auto)\n");
  printf ("   --ccframes   n  number of correlation frames (60)\n");
  printf ("   --framerate  n  stimulus frame rate (119.78 Hz)\n");
  printf ("   --timestep   n  time step in response file (auto)\n");
  printf ("   --mean_bg    n  = mean background (2048)\n");
  printf ("   --rest_v     n  = avg resting potential (0)\n");
  printf ("   --swap       x  = 1, swap big-litte endian (0)\n");
  printf ("\n");
  exit(0);
}

/*------------------------------------*/

void run(int hlp)
{
   void ccrun();

  timepnt = &timebuf;
  startclk = times(timepnt);
  startu = timepnt->tms_utime;
  if (hlp  >= 3) printf ("## %s %s\n",progname,print_version(version));
  ccrun();
}

void prtime(FILE *fout)
{
  totmin = (times(timepnt) - startclk) / SECPERMIN / timecalib;
  runmin = (timepnt->tms_utime - startu) / SECPERMIN / timecalib;
  fprintf (fout,"## time elapsed: %s %.4g  total %.4g minutes.\n",
					progname,runmin,totmin);
  fflush (fout);
}


/*------------------------------------------------------------*/

static stype variables[100] = {0};

void setiptr (int i, char *name, int type, int *iptr)
{
   variables[i].name = name;
   variables[i].type = type;
   variables[i].iptr = iptr;
}

void setcptr (int i, char *name, int type, char *cptr)
{
   variables[i].name = name;
   variables[i].type = type;
   variables[i].cptr = cptr;
}

void setsptr (int i, char *name, int type, char **sptr)
{
   variables[i].name = name;
   variables[i].type = type;
   variables[i].sptr = sptr;
}

void setdptr (int i, char *name, int type, double *dptr)
{
   variables[i].name = name;
   variables[i].type = type;
   variables[i].dptr = dptr;
}


/*------------------------------------------------------------*/

void initvar ()

/* Set the pointers to the variables, which are different types */

{
   int i = 0;

  setsptr(i++, "stimfile",  STRING, &stimfile);
  setsptr(i++, "respfile",  STRING, &respfile);
  setiptr(i++, "nbox",      INT,    &nbox);
  setiptr(i++, "ccframes",  INT,    &ccframes);
  setdptr(i++, "framerate", DOUBLE, &framerate);
  setdptr(i++, "timestep",  DOUBLE, &timestep);
  setdptr(i++, "mean_bg",   DOUBLE, &mean_bg);
  setdptr(i++, "rest_v",    DOUBLE, &rest_v);
  setiptr(i++, "swap",      INT,    &swap);
  setiptr(i++, "info",      INT,    &info);
/*
  setiptr(i++, "nspikbin",  INT,    &nspikbin);
  setiptr(i++, "nbins",     INT,    &nbins);
  setiptr(i++, "ntrials",   INT,    &ntrials);
  setdptr(i++, "temp_freq", DOUBLE, &temp_freq);
  setdptr(i++, "dvdtthresh",DOUBLE, &dvdtthresh);
  setdptr(i++, "spikthresh",DOUBLE, &spikthresh);
  setiptr(i++, "spike_info",INT,    &spike_info);
  setdptr(i++, "ztime",     DOUBLE, &ztime);
*/
}

/*---------------------------------------------------------*/

void setvariable (stype *var)

/* procedure to set variables from command line */

{
    char *npt;
    int found,i;

  for (found=i=0; variables[i].name; i++) {
     if (strcmp(variables[i].name,var->name)==0) {
	found = 1;
	break;
     }
  }
  if (found) {
     if (var->type==STRING) *variables[i].sptr = var->cptr;
      else {
       if      (variables[i].type==DOUBLE) *variables[i].dptr = var->val;
       else if (variables[i].type==INT)    *variables[i].iptr = int(var->val);
      }
   }
}

/*---------------------------------------------------------*/

void setvar (void)

{
    int i;

  for (i=0; i<varset; i++) {
    setvariable(&varval[i]); 
  }
}

/*---------------------------------------------------------*/

double timeval;		/* time value from input file */
double v;		/* voltage value from input file */

#define STRSIZ 400
char linbuf[STRSIZ] = {0};

/*---------------------------------------------------------*/

void readline(FILE *fil)

{
    char *strp;
 while ((strp=fgets(linbuf,STRSIZ,fil)) && linbuf[0]=='#'); /* get 1 line */
 if (!strp) linbuf[0] = 0;
}

/*---------------------------------------------------------*/

int scanline (void) 

   /* read time and volt value */

{
   int n;

  n=sscanf(linbuf,"%lg %lg\n",&timeval,&v);
  return (n);
}

/*---------------------------------------------------------*/

char *print_version (int *ncversion)

{
     static char vbuf[20];

  sprintf (vbuf,"version %d.%d.%d",
                version[0],version[1],version[2]);
  return vbuf;
}

/*---------------------------------------------------------*/

void read_swap (char *pnum, char *p, int size)

/* Swaps number in big-endian format to little-endian format: 
   Format: high byte first. Works for double, float, int, or short int.

   pnum = pointer to output binary number.
   p    = pointer to input number.
   size = size of number in bytes.
 */

{
  int i;

    for (pnum+=size-1,i=0; i<size; i++) {
      *pnum-- = *p++; 		/* reverse order of bytes */
    }
}

/*---------------------------------------------------------*/

void ccrun(void)

{
    int t,vh,oldvh,totspikes,spikenum,nsamples,sframe;
    int i,j,k,ccf,x,y,nstimframes, nstimboxes, respline;
    int filbuf[MATBUFSIZ];
    double vsum, sx, sy, timesteps_frame;
    double mincc, maxcc, ccpeak, ccpeakmax, ccpeakmin, totvol;
    double bigframe, bigvol;
   
  initvar();

   /* set default values */

  if (!stimfile) stimfile = "binlist.mat";
  if (!respfile) respfile = "02419020";
  nbox	    = 0;
  ccframes  = CCFRAMES;
  framerate = FRAMERATE;
  timestep  = 0;		/* set to auto */
  mean_bg   = 2048.0;		/* mean background */
  rest_v    = 0;		/* mean resting potential for memb pot */

  setvar();			/* set variables from command line */

  if (ccframes > CCFRAMES) ccframes = CCFRAMES;

  sfile = fopen (stimfile,"r");
  fread (&filbuf,MATBUFSIZ,sizeof(int), sfile);   /* read in Matlab header, up to data vals */
  if (swap) read_swap((char *)&nstimboxes,(char *)&filbuf[ROWADDR],4);
  else nstimboxes = filbuf[ROWADDR];             /* find number of stimulus boxes */
  if (!nbox) nbox = int(sqrt(nstimboxes)); 
  if (nbox > NBOX) nbox = NBOX;

  if (swap) read_swap((char *)&nstimframes,(char *)&filbuf[COLADDR],4);
  else nstimframes = filbuf[COLADDR];	   /* find number of stimulus frames */
  // fprintf (stderr,"C %d %d\n",nstimboxes,nstimframes);

  for (i=0; i<nstimframes; i++) {
    for (x=0; x<nbox; x++) {
      for (y=0; y<nbox; y++) {
       k = 0;
       fread (&k, sizeof(short int),1, sfile);
       if (swap) read_swap ((char *)&stimlist[i][x][y], (char *)&k, 2);
       else stimlist[i][x][y] = k; 
       // fprintf (stderr,"%d ",stimlist[i][x][y]);
      }
    }
  }
  if (info >=2) fprintf (stderr,"# Stimlist rows %d cols %d data %hd %hd\n", 
				nstimboxes, nstimframes, 
				stimlist[0][0][0], stimlist[0][0][1]);

  for (i=0; i<ccframes; i++) {		/* zero the array */
    for (x=0; x<nbox; x++) {
      for (y=0; y<nbox; y++) {
       ccout[i][x][y] = 0;
      }
    }
  }

  if ((rfile=fopen(respfile,"r")) == NULL) {
      fprintf (stderr,"%s: can't open file '%s'.\n",progname,respfile);
      return;
  }

  readline(rfile);
  if (scanline() < 2) {
	fprintf (stderr,"%s: empty file '%s'.\n",progname,respfile);
	return;
  }
  t1 = timeval;
  readline(rfile);
  if (scanline() < 2) {
	fprintf (stderr,"%s: empty file '%s'.\n",progname,respfile);
	return;
  }
  t2 = timeval;
  if (timestep == 0) timestep = t2 - t1;
  rewind(rfile);
  readline(rfile);
  if (framerate <=0) framerate = 1;
  if (timestep  <=0) {
     fprintf (stderr,"%s: zero timestep.\n",progname);
     return;
  }
  timesteps_frame = 1.0/framerate/timestep;	/* timesteps per frame */
  if (timesteps_frame <= 0) timesteps_frame = 1;
  if (info>=2) fprintf (stderr,"timesteps_frame %g\n",timesteps_frame);

   for (respline=0; scanline() >= 2; readline(rfile),respline++) {

    vsum = (v - rest_v);   

    sframe = int(respline/timesteps_frame);
    if (sframe > TOTFRAMES) sframe = TOTFRAMES;

    vsum /= mean_bg * nstimframes;

    //vsum = stimlist[sframe-5][2][2];              /* test pattern */

    for (j=0; j<ccframes; j++) {
        if ((ccf=sframe-j) < 0) ccf += 10000; /* wrapping up in beginning of file */ 
        for (x=0; x<nbox; x++) {
          for (y=0; y<nbox; y++) {
               ccout[j][x][y] += (stimlist[ccf][x][y]-mean_bg) * vsum; /* output looks back */
             }
	  }
    }

  } /* for (line=0; ; ) */

  maxcc = -1e20;
  mincc =  1e20;
  for (j=0; j<ccframes; j++) {
     for (x=0; x<nbox; x++) {
        for (y=0; y<nbox; y++) {
          if (maxcc < ccout[j][x][y]) maxcc = ccout[j][x][y];
          if (mincc > ccout[j][x][y]) mincc = ccout[j][x][y];
        }
     }
  }

  ccpeakmax = abs(maxcc);
  ccpeakmin = abs(mincc);

  ccpeak = (ccpeakmax>ccpeakmin ? ccpeakmax : ccpeakmin);

  if (info>=2) fprintf (stderr,"ccpeak %g\n",ccpeak);

  for (j=0; j<ccframes; j++) {
     for (x=0; x<nbox; x++) {
        for (y=0; y<nbox; y++) {
	  ccout[j][x][y] /= ccpeak;
	}
     }
  }	  

/*
  fprintf (stderr,"# totvol %g before norm (/scontrast/volt/sec).\n",totvol/(((double)ccframes)/framerate));
  fprintf (stderr,"# ccpeak %g after norm.\n",ccpeak);
*/
  /* print out cross-correlation frames */
 
  for (j=0; j<ccframes; j++) {
     int prframe;
     prframe = ccframes-1-j; 
     for (x=0; x<nbox; x++) {
        for (y=0; y<nbox; y++) {
//		printf("%9.3g ", ccout[prframe][x][y]); /* output looks forward */
		printf("%9.3g ", ccout[j][x][y]); 	/* output looks backward */
	}
     //printf("\n");
     }
     printf("\n");
  }	  


  /* find frame with biggest volume response */

if (0) {
  bigvol = -1e20;
  bigframe = 0;
  for (j=0; j<ccframes; j++) {
     totvol = 0;
     for (x=0; x<nbox; x++) {
        for (y=0; y<nbox; y++) {
          totvol += abs(ccout[j][x][y]);
        }
     }
     if (bigvol < totvol) {	/* remember biggest one */
	 bigvol = totvol;
         bigframe = j;
     }
  }
    /* find peak height of frame with biggest volume */

  for (j=0; j<ccframes; j++) {
  maxcc = -1e20;
  mincc =  1e20;
  for (x=0; x<nbox; x++) {
     for (y=0; y<nbox; y++) {
       if (maxcc < ccout[j][x][y]) maxcc = ccout[j][x][y];
       if (mincc > ccout[j][x][y]) mincc = ccout[j][x][y];
     }
  }

  if (maxcc > mincc) ccpeak = maxcc;
  else               ccpeak = mincc;

/* print out biggest frame, normalized to 255 */

/*
     fprintf(stderr,"P2\n");
     fprintf(stderr,"%d %d\n",nbox,nbox);
     fprintf(stderr,"%d\n",255);
*/  
   for (x=0; x<nbox; x++) {
        for (y=0; y<nbox; y++) {
/*		fprintf(stderr,"%d ", int(ccout[j][x][y]/ccpeak*255.0)); */
	}
     fprintf(stderr,"\n");
     }
  } 
 }
}

/*---------------------------------------------------------*/
