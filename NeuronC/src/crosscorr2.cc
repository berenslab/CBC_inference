/* Program crosscorr2 */

/* Computes crosscorrelation between 2 columns of a data file */
/* Autocorrelation added for either column */
/* R.G.Smith Nov, 2009 */

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


#define STRING 2		/* used in stype below */
#define NUMBER 3
#define INT    4
#define DOUBLE 5
#define VVSIZE 30
#define TIME 0


struct stype {
	const char *name;
	int type;
	union {
	    double val;
	    int    *iptr;
	    char   *cptr;
	    const char  **sptr;
	    double *dptr;
	};
};

stype varval[VVSIZE] = {0};	/* temp values for setting variables */
int varset = 0;			/* index for setting variables */

int version[3];

#define TIMCOL 0	/* column number for time in data file */
#define   VCOL 1	/* column number for voltage in data file */

#define MAXTRACES 50
double traces[MAXTRACES] = {0};

/* command line variables */

double framerate;
double correllen = 0.5;		/* length of correlation */
double timestep=0;
double vscale = 1;		/* voltage scale, mV -> V */
double dvdtthresh = 60.0;	/* dv/dt threshold (V/sec) */
double spikthresh = 10.0;	/* voltage threshold (V) */
double tstart = 0.005;		/* time to start at */
double tend = 10000;		/* time to end at */
int arg1;			/* number of first y-axis trace to be correlated */
int arg2;			/* number of second y-axis trace to be correlated */
int arg1m=0;			/* arg 1 mode: 0=voltage, 1=spike */
int arg2m=0;			/* arg 2 mode: 0=voltage, 1=spike */
int autocorr=0;			/* autocorrelogram mode, use as arg number */
int info;
int binsiz=1;			/* size of bins in counts of timestep */
double binlen=0;		/* size of bins in time */


#define FILENAMLEN 100

static const char *infile = {NULL};		/* file containing 2 colums for crosscorr */
char infilem [FILENAMLEN] = {0};

FILE *ifile;

void run(int hlp);
void prtime(FILE *stimout);
void printhelp();
char *print_version(int *version);

/*------------------------------------*/

main(int argc, char **argv)
{
   char *cptr;

 ifile = stdin;
 version[0] = 2;
 version[1] = 0;
 version[2] = 5;
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
		  if (isalpha(**argv)) {
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
		  if (isalpha(**argv)) {
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
       infile=cptr;
       argv++; argc--;
       cptr = *argv;

   /*       if((fin=fopen(infile=cptr,"r"))==NULL)
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
  printf ("   --infile     s  input file (1 time, 2 data cols)\n");
  printf ("   --timestep   n  time step (auto)\n");
  printf ("   --binlen     n  size of time bins (timestep)\n");
  printf ("   --correllen  n  length of correlation array (0.5 s)\n");
  printf ("   --arg1       n  col for first trace to be correlated  (1)\n");
  printf ("   --arg1m      n  arg1 mode, 0->voltage, 1->spikes (0)\n");
  printf ("   --arg2       n  col for second trace to be correlated (2)\n");
  printf ("   --arg2m      n  arg2 mode, 0->voltage, 1->spikes (0)\n");
  printf ("   --autocorr   n  autocorrel with trace (0)\n");
  printf ("   --vscale     n  voltage scale, mV -> V (1)\n");
  printf ("   --dvdtthresh n  dv/dt threshold for spikes (60 V/s)\n");
  printf ("   --spikthresh n  voltage threshold for spikes (10 V)\n");
  printf ("   --tstart     n  time to start with (0.005)\n");
  printf ("   --tend       n  time to end with (1000)\n");
  printf ("   --info       n  info printout (0)\n");
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
  if (hlp==0) ccrun();
  else printhelp();
  prtime(stdout);
}

void prtime(FILE *fout)
{
  totmin = (times(timepnt) - startclk) / SECPERMIN / timecalib;
  runmin = (timepnt->tms_utime - startu) / SECPERMIN / timecalib;
  fprintf (fout,"## time elapsed: %s %.4g  total %.4g minutes.\n",
					progname,runmin,totmin);
  fflush (fout);
}


/*---------------------------------------------------------*/

char *print_version (int *ncversion)

{
     static char vbuf[20];

  sprintf (vbuf,"version %d.%d.%d",
                version[0],version[1],version[2]);
  return vbuf;
}

/*------------------------------------------------------------*/

static stype variables[100] = {0};

void setiptr (int i, const char *name, int type, int *iptr)
{
   variables[i].name = name;
   variables[i].type = type;
   variables[i].iptr = iptr;
}

void setcptr (int i, const char *name, int type, char *cptr)
{
   variables[i].name = name;
   variables[i].type = type;
   variables[i].cptr = cptr;
}

void setsptr (int i, const char *name, int type, const char **sptr)
{
   variables[i].name = name;
   variables[i].type = type;
   variables[i].sptr = sptr;
}

void setdptr (int i, const char *name, int type, double *dptr)
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

  setsptr(i++, "infile",  STRING,   &infile);
  setdptr(i++, "timestep",  DOUBLE, &timestep);
  setdptr(i++, "correllen", DOUBLE, &correllen);
  setiptr(i++, "arg1",      INT,    &arg1);
  setiptr(i++, "arg1m",     INT,    &arg1m);
  setiptr(i++, "arg2",      INT,    &arg2);
  setiptr(i++, "arg2m",     INT,    &arg2m);
  setiptr(i++, "autocorr",  INT,    &autocorr);
  setdptr(i++, "vscale",    DOUBLE, &vscale);
  setdptr(i++, "dvdtthresh",DOUBLE, &dvdtthresh);
  setdptr(i++, "spikthresh",DOUBLE, &spikthresh);
  setdptr(i++, "tstart",    DOUBLE, &tstart);
  setdptr(i++, "tend",      DOUBLE, &tend);
  setiptr(i++, "info",      INT,    &info);
  setdptr(i++, "binlen",    DOUBLE, &binlen);
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
double v1;		/* voltage value from input file */
double v2;		/* voltage value from input file */

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

// int scanline__ (FILE *fd, double *args) 
// 
//    /* read time and volt value */
// 
// {
//    int n;
// 
//   readline(fd);
//   n=sscanf(linbuf,"%lg %lg %lg\n",&timeval,&v1, &v2);
//   return (n);
// }

/*---------------------------------------------------------*/

int scanline (FILE *fd, double *args) 

   /* Read time and all the y values from next line */
   /* return the number of y values */

{
    int i,n;
    double val;
    static double oldtimeval=0;
    static char *str,*ptr;

  str = linbuf;
  ptr = (char *)NULL;
  readline(fd);
  args[0] = timeval = strtod(str,&ptr);
  for (n=1; n < MAXTRACES; n++) {
    str = ptr;
    val = strtod(str,&ptr);
    if (str != ptr) {
      args[n] = val;                /* read in plots from one line */
    }
    else
      args[n] = 0.0;                /* if not enough, then insert a 0 */

    if (str==ptr) break;

  } /* for (n;;) */

  if (n < 1) {
    timeval = oldtimeval;
  }
 
 oldtimeval = timeval;
 return (n-1);
}

/*---------------------------------------------------------*/

FILE *openfil (const char *filnam, char *infilm)

/* Open file for reading.  Check if file is .gz or .bz2, 
   or has a .gz or .bz2 version in same directory.  If so,
   open pipe to read the file with zcat or bzcat. Store filename
   actually used in '*infilm'.
*/

{
   static char infilz[FILENAMLEN];
   FILE *fildesc;

     /* check if filename is .gz */

   if (info>=3) printf ("# opening file '%s'\n",filnam); 

   strcpy(infilm,filnam);

   if (strstr(filnam,".gz") != NULL) {
       strcpy(infilm,filnam);
       strcpy(infilz,"zcat ");
       strncat(infilz,filnam,FILENAMLEN-5);
       if ((fildesc=popen(infilz,"r"))==NULL) {
          fprintf (stderr,"%s: can't open file '%s'\n",progname,filnam);
          fildesc = NULL;
       }
   }
      /* check if filename is .bz2 */

   else if (strstr(filnam,".bz2") != NULL) {
       strcpy(infilm,filnam);
       strcpy(infilz,"bzcat ");
       strncat(infilz,filnam,FILENAMLEN-5);
       if ((fildesc=popen(infilz,"r"))==NULL) {
          fprintf (stderr,"%s: can't open file '%s'\n",progname,filnam);
          fildesc = NULL;
       }
   }
      /* Try to open non .gz file, then check if .gz version exists */
 
   else if ((fildesc=fopen(filnam,"r")) == NULL) {
      strcpy(infilz,filnam);
      strncat(infilz,".gz",FILENAMLEN-5);
      if ((fildesc=fopen(infilz,"r"))!=NULL){ /* look for .gz file */
         if (info>=3) printf ("# trying  file '%s'\n",infilz); 
         fclose(fildesc);
         strcpy(infilm,infilz);
         strcpy(infilz,"zcat ");          /* check for .gz file */
         strncat(infilz,filnam,FILENAMLEN-5);
         strncat(infilz,".gz",FILENAMLEN-5);
         if ((fildesc=popen(infilz,"r"))==NULL) {
            fprintf (stderr,"%s: can't open file '%s'\n", progname,filnam);
            fildesc = NULL;
         }
      }
      else {   /* Try to open non .bz2 file, then check if .bz2 version exists */
        strcpy(infilz,filnam);
        strncat(infilz,".bz2",FILENAMLEN-5);
        if ((fildesc=fopen(infilz,"r"))!=NULL){ /* look for .gz file */
           if (info>=3) printf ("# trying  file '%s'\n",infilz); 
           fclose(fildesc);
           strcpy(infilm,infilz);
           strcpy(infilz,"bzcat ");          /* check for .gz file */
           strncat(infilz,filnam,FILENAMLEN-5);
           strncat(infilz,".bz2",FILENAMLEN-5);
           if ((fildesc=popen(infilz,"r"))==NULL) {
              fprintf (stderr,"%s: can't open file '%s'\n", progname,filnam);
              fildesc = NULL;
           }
        }
        else {
          fprintf (stderr,"%s: can't find file '%s'\n",progname,filnam);
          fildesc=NULL;
        }
      }
   }    /* if (fopen(filnam,"r") == NULL) */

   //if (info>=3) printf ("# saving filename '%s' for rewind.\n",infilm); 

   return fildesc;
}


/*---------------------------------------------------------*/

FILE *rewindfil (FILE *fildesc, const char *filnam)

/* Rewind file that could have been a pipe. */
/* If filnam contains ".gz" or ".bz2", close with "pclose". */

{
   static char infilz[FILENAMLEN];

   if (strstr(filnam,".gz") != NULL) {
       pclose(fildesc);
       strcpy(infilz,"zcat ");
       strncat(infilz,filnam,FILENAMLEN-5);
       if ((fildesc=popen(infilz,"r"))==NULL) {
        fprintf (stderr,"%s: rewind popen can't open file '%s'\n",progname,filnam);
        fildesc = NULL;
       }
   }
      /* check if filename is .bz2 */

   else if (strstr(filnam,".bz2") != NULL) {
       pclose(fildesc);
       strcpy(infilz,"bzcat ");
       strncat(infilz,filnam,FILENAMLEN-5);
       if ((fildesc=popen(infilz,"r"))==NULL) {
        fprintf (stderr,"%s: rewind popen can't open file '%s'\n",progname,filnam);
        fildesc = NULL;
       }
   }
      /* filename is not .gz */
 
   else rewind(fildesc);
   return fildesc;
}

/*---------------------------------------------------------*/

void closefil (FILE *fildesc, char *filnam)

/* Close file that could have been a pipe. */
/* If filnam contains ".gz" or ".bz2", close with "pclose". */

{
   static char infilz[FILENAMLEN];

   if (strstr(filnam,".gz") != NULL) {
       pclose(fildesc);
   }
      /* check if filename is .bz2 */

   else if (strstr(filnam,".bz2") != NULL) {
       pclose(fildesc);
   }
      /* filename is not .gz */
 
   else fclose(fildesc);
}

/*---------------------------------------------------------*/

void reset_file ()

{
   int stim;

  ifile = rewindfil (ifile,infile);	/* rewind file back to start */
}

/*---------------------------------------------------------*/

double oldv[3]={0};
int oldvh[3]={0};
int nspikes[3]={0};

double getspik(double v, int arg)

{
    int vh;
    double vout;

  vh = (((v - oldv[arg])/timestep/vscale > dvdtthresh) || (v > spikthresh));
  if (vh && !oldvh[arg]) {           /* spike here */
    vout = 1;
    nspikes[arg]++;
  }
  else vout = 0;
  oldv[arg]  = v;
  oldvh[arg] = vh;
  return vout;
}

/*---------------------------------------------------------*/

void resetspik(void)
{
   int i;
 for (i=0; i<3; i++) {
    oldv[i]    = 0;
    oldvh[i]   = 0;
    nspikes[i] = 0;
 }
}

/*---------------------------------------------------------*/

void ccrun(void)

{
    int i,j,k,l,x,y;
    int maxpoints,ncols,npoints,obufsiz;
    int bincnt,correll,indx;
    double sum1, sum2, sumsq1, sumsq2;
    double mean1, mean2, val1, val2;
    double var1, var2;
    double meanvar, scale;
    double *buf1, *buf2, *obuf;
    double t,ts,oldtimeval;

  initvar();

   /* set default values */

  timestep  = 0;		/* set to auto */
  info = 0;
  arg1 = 1;
  arg2 = 2;
  
  setvar();			/* set variables from command line */

  ifile = openfil(infile,infilem);
  if (timestep==0) {
    scanline(ifile,traces);
    oldtimeval = timeval;
    scanline(ifile,traces);
    timestep = timeval - oldtimeval;
    rewindfil(ifile,infile);
    if (binlen==0) binlen = timestep;
  }
  if (timestep==0) timestep = 0.001;			/* prevent divide by zero */
  ts = tstart*timestep;
  for (t=0; t<ts; t+=timestep) {			/* skip time at beginning */
    scanline(ifile,traces);
  }
  if (info>=2) fprintf (stderr,"# correllen %g timestep %g\n",correllen, timestep);
  if (info>=2) fprintf (stderr,"# arg1 %d arg2 %d\n",arg1,arg2);

  binsiz = binlen / timestep + 0.5;			/* size of bins in terms of timestep */
  correll = correllen / binlen + 0.5;				/* size of correll array */

  if (info>=2) fprintf (stderr,"# binlen %g binsiz %d\n",binlen, binsiz);

  sum1=sum2=sumsq1=sumsq2=0;
  val1=val2=0;
  resetspik();
  for (i=indx=bincnt=0;scanline(ifile,traces)>0 && t<tend;i++,t+=timestep) { /* count number of points in file*/

    if (autocorr==2) v1 = traces[arg2];
    else             v1 = traces[arg1];
    if (arg1m)       v1 = getspik(v1,1);	     /* translate into spike: 0-1 */

    if (autocorr==1) v2 = traces[arg1];
    else             v2 = traces[arg2];
    if (arg2m)       v2 = getspik(v2,2);	     /* translate into spike: 0-1 */

    val1 += v1;
    val2 += v2;
    if (++bincnt >= binsiz) {
	  bincnt = 0;
	  indx++;
	  sum1   += val1;
	  sumsq1 += val1*val1;
	  sum2   += val2;
	  sumsq2 += val2*val2;
	  val1 = val2 = 0;
    }
    //if (info>=2) fprintf (stderr,"arg1 %g arg2 %g\n",traces[arg1],traces[arg2]);
    //if (info>=2) fprintf (stderr,"arg1 %g arg2 %g\n",v1,v2);
  }
  npoints = indx;
  mean1 = sum1 / npoints;
  mean2 = sum2 / npoints;
  var1 = (sumsq1 - (sum1 * mean1)) / npoints;
  var2 = (sumsq2 - (sum2 * mean2)) / npoints;
  meanvar = sqrt (var1*var2);
  if (info>=2) fprintf (stderr,"# mean1 %g mean2 %g\n",mean1, mean2);
  if (info>=2) fprintf (stderr,"# mean var %g\n",meanvar);

  obufsiz  = correll*2;
  buf1 = (double *) malloc (npoints*sizeof(double));
  buf2 = (double *) malloc (npoints*sizeof(double));
  obuf = (double *) malloc (obufsiz*sizeof(double));

  for (i=0; i<npoints; i++) {				/* zero the buffers to allow building average */
    buf1[i] = 0;
    buf2[i] = 0;
  }
  for (i=0; i<obufsiz; i++) {
    obuf[i] = 0;
  }

  rewindfil(ifile,infile);
  ts = tstart*timestep;
  for (t=0; t<ts; t+=timestep) {			/* skip time at beginning */
    scanline(ifile,traces);
  }
  resetspik();
  for (i=indx=bincnt=0;scanline(ifile,traces)>0 && t<tend; i++,t+=timestep) {    /* read in points minus mean */
    if (autocorr==2) v1 = traces[arg2];
    else             v1 = traces[arg1];
    if (autocorr==1) v2 = traces[arg1];
    else             v2 = traces[arg2];
    if (arg1m) {
	v1 = getspik(v1,1);	     /* translate into spike: 0-1 */
        buf1[indx] += v1;
    } else buf1[indx] += v1-mean1;
    if (arg2m) {
	v2 = getspik(v2,2);	     /* translate into spike: 0-1 */
        buf2[indx] += v2;
    }
    else buf2[indx] += v2-mean2;
    if (++bincnt >= binsiz) {
	  bincnt = 0;
          indx++;
    }
  } 
  // if (info>=2) fprintf (stderr,"dvdtthresh %g spikthresh %g\n",dvdtthresh,spikthresh);
  if (info>=2) if (arg1m || arg2m) 
	  fprintf (stderr,"# nspikes1 %d nspikes2 %d\n",nspikes[1],nspikes[2]);

  for (i=0; i<obufsiz; i++)	/* zero the array */ 
    obuf[i]= 0;
  
  for (j=0; j<npoints; j++) {
    for (k=0,x=j-correll; k<obufsiz; k++,x++) {
      if (x>=0 && x<npoints) { obuf[k] += buf1[j] * buf2[x]; }
    }
  }

  scale = npoints*meanvar;
  //scale = npoints;
  if (scale==0) scale = 1; 
  for (i=0,x=-correll; i<obufsiz; i++,x++) {	/* normalize and print the array */ 
    printf ("%8.4g %5.4g\n",x*binlen, obuf[i]/scale);
  }
  //if (info>=2) fprintf (stderr,"scale %g %g %d\n",scale, meanvar,npoints);

}

/*---------------------------------------------------------*/
