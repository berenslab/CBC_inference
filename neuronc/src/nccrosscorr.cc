
/* module crosscorr */

/* Computes crosscorrelation between spatio-temporal stimulus and response */

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

void ccorr(double *stimarr, int xn, int yn, int nframes, double timestep,
		double *response_arr, int nresp, int ccframes, double **ccout)

{
    int t,vh,oldvh,totspikes,spikenum,nsamples,sframe;
    int i,j,k,n,ccf,x,y,nstimframes, nstimboxes, respline;
    int framesiz, ccout_size;
    double vsum, sx, sy, timesteps_frame;
    double mincc, maxcc, ccpeak, ccpeakmax, ccpeakmin, totvol;
    double bigframe, bigvol;
    double *parr;


  if (timestep  <=0 || nresp<=0) {
     ncfprintf (stderr,"ccorr: zero timestep.\n");
     return;
  }

   /* set default values */

  if (info >=2) fprintf (stderr,"# Stimlist rows %d cols %d data %hd %hd\n", 
				nstimboxes, nstimframes, 
				stimlist[0][0][0], stimlist[0][0][1]);

  framesiz = xn * yn;
  ccout_size = ccframes * framesiz;

  if ((ccout = (double*)emalloc(ccout_size * sizeof(double)))==NULL) {
       ncfprintf (stderr,"nc: crosscorr, can't init output array\n");
       return;
  }

  for (parr=ccout,i=0; i<ccframes; i++) {	/* zero the array */
    for (y=0; y<yn; y++) {
      for (x=0; x<xn; x++) {
        *parr++ = 0;
      }
    }
  }

  timesteps_frame = 1.0/framerate/timestep;	/* timesteps per frame */
  if (timesteps_frame <= 0) timesteps_frame = 1;
  if (info>=2) fprintf (stderr,"timesteps_frame %g\n",timesteps_frame);

   for (n=0; n<nresp; n++) {

    vnorm = (v - norm_v);
    vnorm /= mean_bg * nstimframes;
    sframe = n / timesteps_frame;

    //vnorm = stimarr[sframe-5][2][2];              /* test pattern */

    for (j=0; j<ccframes; j++) {
        for (x=0; x<nbox; x++) {
          for (y=0; y<nbox; y++) {
		ccf = sframe - j;
               ccout[j][y][x] += (stimarr[ccf][y][x]-mean_bg) * vnorm; /* output looks back */
             }
	  }
    }

  } /* for (line=0; ; ) */

  maxcc = -1e20;
  mincc =  1e20;
  for (j=0; j<ccframes; j++) {
     for (y=0; y<nbox; y++) {
        for (x=0; x<nbox; x++) {
          if (maxcc < ccout[j][y][x]) maxcc = ccout[j][y][x];
          if (mincc > ccout[j][y][x]) mincc = ccout[j][y][x];
        }
     }
  }

  ccpeakmax = abs(maxcc);
  ccpeakmin = abs(mincc);

  ccpeak = (ccpeakmax>ccpeakmin ? ccpeakmax : ccpeakmin);

  if (info>=2) fprintf (stderr,"ccpeak %g\n",ccpeak);

  for (j=0; j<ccframes; j++) {
     for (y=0; y<nbox; y++) {
        for (x=0; x<nbox; x++) {
	  ccout[j][y][x] /= ccpeak;
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
     for (y=0; y<nbox; y++) {
        for (x=0; x<nbox; x++) {
//		printf("%9.3g ", ccout[prframe][y][x]); /* output looks forward */
		printf("%9.3g ", ccout[j][y][x]); 	/* output looks backward */
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
     for (y=0; y<nbox; y++) {
        for (x=0; x<nbox; x++) {
          totvol += abs(ccout[j][y][x]);
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
  for (y=0; y<nbox; y++) {
     for (x=0; x<nbox; x++) {
       if (maxcc < ccout[j][y][x]) maxcc = ccout[j][y][x];
       if (mincc > ccout[j][y][x]) mincc = ccout[j][y][x];
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
     for (y=0; y<nbox; y++) {
        for (x=0; x<nbox; x++) {
/*		fprintf(stderr,"%d ", int(ccout[j][y][x]/ccpeak*255.0)); */
	}
     fprintf(stderr,"\n");
     }
  } 
 }
}

/*---------------------------------------------------------*/
