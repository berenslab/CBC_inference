
/* Stimulus functions for C++ NC script language */
/* Calls funcs in "stim" or "nc" */

#include "ncsub.h"
#include "ncelem.h"
#include "control.h"
#include "nc.h"
#include "y.tab.h"
#include "ncomp.h"
#include "ncplot.h"
#include "convarr.h"
#include "ndef.h"
#include "nconst.h"
#include "drand.h"


#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <signal.h>
#include <setjmp.h>
#include <sys/times.h>
#include <sys/types.h>
#include <math.h>
//#include "stdplt.h"

#ifdef __cplusplus
}
#endif

#include "ncio.h"

extern node *nodepnt;           /* pointer to node list */
extern photrec *recpnt;        /* pointer to photrec list */
extern elem *elempnt;           /* pointer to current element */

extern plotfr plotnod[PLOTNODSIZ];  /* holds nodes to be recorded */
extern int numplots;                /* number of plots to be displayed */
extern int runyet;
extern char *rootframe;             /* name of root plot frame */
extern int debug;
extern int debugz;
extern int setdebug;
extern int setdebgz;
extern int setdisp;
extern int istat;
extern int makestim;            /* = 1 -> says we're running "stim", not "nc" */
extern char *progname;
extern int stimfflg;
extern char *stfile;                 /* stimulus file for plotinit */
 
extern FILE *stimout;
extern FILE *stimin;                 /* stimulus file */

#define SECPERMIN 60.0          /* number of seconds per minute */

/* #define TIMECALIB 60.0       /* times() system call gives 1/60  sec */
/* #define TIMECALIB 100.0      /* times() system call gives 1/100 sec */
/* #define TIMECALIB 1024.0     /* times() system call gives 1/1024 sec */

#ifdef DEC_ALPHA
#define TIMECALIB 1024.0        /* times() system call gives .9765625 msec */
#else
#define TIMECALIB 100.0         /* times() system call gives 1/100 sec */
#endif

extern struct tms timebuf;
extern struct tms *timepnt;
extern double runmin;
extern double timecalib;
extern double startclk,startu,totmin;

recstim *makvstim(double start, int nodnm1, int nodnm2, int nodnm3, int nodnm4,
                        double value, int bgend, const char *action);
recstim *makrstim(double start, int recnm1, int recnm2, int recnm3, int recnm4,
                double inten, double wavel, double mask, int stimchan, int bgend, const char *action);
recstim *makrstim(double start, double dur, photrec *rpnt,
                double inten, double wavel, double mask, int stimchan, const char *action);
recstim *makrstim(double start, double dur, photorec *ept,
                double inten, double wavel, double mask, int stimchan, const char *action);

extern elem *elempnt;           /* pointer to current element */
void findconnect(int donefl);

void stcomment (void);
int makblur(double rad, double blur_ampl, double scale, double scatter_rad,
                        double scatter_pow, double scatter_ampl);
node *maknod(nodeint nodea, nodeint nodeb, nodeint nodec,
        nodeint noded);
node *findnode(nodeint node1a, nodeint node1b, nodeint node1c,
        nodeint node1d, const char *s);
elem *makelem(int etype, elem *epnt);

int inrect(double x, double y, double maxx, double minx,
                        double maxy, double miny, double orient);
int incirc(double x, double y, double circx,
                        double circy, double rad);
/* functions for stim */
void makcheckerboard(int arr, int arrsize, double xsize, double ysize, int xn, int yn,
		double xloc, double yloc, double xcent, double ycent,
		double scale, double orient, double tfreq,
		double inten, double contrast, double start, double dur,
		double **stim_rndarr, int *stim_nfr, int ckrseed);

void makgabor(int arr, int arrsize, double contrast, double inten_add, double inten_mult,
		double xloc, double yloc, double xcent, double ycent,
		double xenv, double yenv, int makenv, int sq,
		double tfreq, double drift, double speriod, double sphase,
		double scale, double orient, double wavel,
		double mask, int stimchan, double start, double dur);

void maksineann(int arr, int arrsize, double contrast, double inten_add, double inten_mult,
		double xloc, double yloc, double xcent, double ycent,
		double renv, int makenv, int sq,
		double tfreq, double drift, double speriod, double sphase,
		double scale, double wavel, double mask, int stimchan,
		double start, double dur);

void makwindmill(int arr, int arrsize, double contrast, double inten_add, double inten_mult,
		double xloc, double yloc, double xcent, double ycent,
		double renv, int makenv, int sq,
		double tfreq, double drift, double speriod, double sphase,
		double scale, double start, double dur, double wavel, double mask, int stimchan);

void makimage(int arr, int arrsize, double width, double length,
		double xoff, double yoff, double orient, double scale,
		double inten, int type, const char *imgfil);

void makspot(int arr, int arrsize, double dia,
                double xloc, double yloc, double xcent, double ycent,
                double scale, double inten, double start, double dur,
                double wavel, double mask, int stimchan, int invert);

void makbar(int arr, int arrsize, double width, double length, 
                double xloc, double yloc, double xcent, double ycent,
                double scale, double orient, double inten,              
                double start, double dur, double wavel, double mask, int stimchan);

void stimclampout(double time, node *npt, double value, int bgend, const char *action);

void find_photarrsiz(double *xmax,double *xmin,double *ymax,double *ymin);
void execerror (const char *str1, const char *str2);
double axy(double x, double y);
char *emalloc(unsigned int n);
char *makrand (int rndsiz, const char *mesg, int num);
void recback(int array, double inten, double wavel);
void abslist(double ratio, double time, double mask, int stimchan);
void initrand (int ngen, int nrseed);

/*---------------------------------------------------------------------*/

void vclamp (node *npnt, double inten, double start, double dur)
{
 if (npnt==NULL) { ncfprintf (stderr,"vclamp: null node\n"); return; }
 if (makestim) {
   stimclampout(start,     npnt,inten,0,"evclamp");
   stimclampout(start+dur, npnt,inten,1,"evclamp");
 }
 else if (dur > 0 && !makestim) {
   makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,inten,0,"evclamp");
   makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,inten,1,"evclamp");
 }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void cclamp (node *npnt, double inten, double start, double dur)
{
 if (npnt==NULL) { ncfprintf (stderr,"cclamp: null node\n"); return; }
 if (makestim) {
   stimclampout(start,     npnt,inten,0,"gcclamp");
   stimclampout(start+dur, npnt,-inten,1,"gcclamp");
 }
 else if (dur > 0 && !makestim) {
  makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,inten,0,"gcclamp");
  makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-inten,1,"gcclamp");
 }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void puff (node *npnt, int puffmsg, double puffconc, double start, double dur)

/* puff neurotrans/2nd msng */

{
 if (dur > 0 && !makestim) {
  switch (puffmsg) {
    case GLU:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"i");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"i");
    break;
    case AMPA:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"j");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"j");
   break;
    case KAINATE:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"k");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"k");
   break;
   case NMDA:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"l");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"l");
   break;
   case CNQX:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"m");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"m");
   break;
   case GABA:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"n");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"n");
   break;
   case BIC:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"o");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"o");
   break;
   case PTX:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"p");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"p");
   break;
   case GLY:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"q");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"q");
   break;
   case STRY:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"r");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"r");
   break;
   case CAMP:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"s");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"s");
   break;
   case CGMP:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"t");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"t");
   break;
   case CA:
      makvstim(start,     npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, puffconc,0,"u");
      makvstim(start+dur, npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-puffconc,1,"u");
   break;
   }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_node (node *npnt, double inten, double start, double dur)
{
   const char *action;
   double wavel=1;
   double mask=0;

  action = "d";			// additive
  makrstim(start,    npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, inten,wavel,mask,0,0,action);
  makrstim(start+dur,npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-inten,wavel,-mask,0,1,action);
}
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_node (node *npnt, double inten, double start, double dur, double wavel)
{
   const char *action;
   double mask=0;

  action = "d";			// additive
  makrstim(start,    npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4, inten,wavel,mask,0,0,action);
  makrstim(start+dur,npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4,-inten,wavel,-mask,0,1,action);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_cone (node *npnt, double inten, double start, double dur)
{
   stim_node (npnt, inten, start, dur);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_cone (node *npnt, double inten, double start, double dur, double wavel)
{
   stim_node (npnt, inten, start, dur, wavel);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_rod (node *npnt, double inten, double start, double dur, double wavel)
{
   stim_node (npnt, inten, start, dur, wavel);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_rod (node *npnt, double inten, double start, double dur)
{
   stim_node (npnt, inten, start, dur);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

extern int stimhlp;		/* in stimsub.cc */

void stim_file (const char *filename)

{
#define  STFILNAM 200
    int compress_type;
    static char stimfile[STFILNAM]={0};
    static char stimfilz[STFILNAM]={0};
    FILE *ftemp;

  if (makestim) {			/* "stim" running, make stim file */
	  				/* User can compress to make .gz, bz2, or .xz version */
     if (*stimfile) fclose (stimout);
     *stimfile = 0;
     strcpy (stimfile,filename);	/* copy filename */
     if (stimhlp) ncfprintf (stderr,"stimfil '%s'\n",stimfile);
     if (*stimfile) {
          stimfflg = 1;
          if ((ftemp=fopen(stimfile,"w")) == NULL) {
            ncfprintf (stderr,"%s: can't open file '%s'\n", progname,stimfile);
            return;
          }
          else stimout = ftemp;
     }
     else {
           stimfflg = 0;
           stimout = stdout;
     }
     stcomment();			/* print comment in stim file */
  }
  else	/* nc running */

        /* Look up the stim file.  
	   If it doesn't exist, look for .gz, .bz2, or .xz compressed versions.  */
  {
    if (stimfflg && *stimfile) fclose(stimin); /* close previous */
    strcpy (stimfile,filename);  	/* copy filename */
    stfile = stimfile;			/* for plotinit */
    if (! *filename) { stimfflg = 0; return; }
    if ((ftemp=fopen(stimfile,"r"))==NULL) { /* open file */
      strcpy(stimfilz,stimfile);
      strncat(stimfilz,".gz",STFILNAM-5);
      compress_type = 1;
      if ((ftemp=fopen(stimfilz,"r"))==NULL){ /* look for .gz file */
         strcpy(stimfilz,stimfile);
         strncat(stimfilz,".bz2",STFILNAM-5);
         compress_type = 2;
         if ((ftemp=fopen(stimfilz,"r"))==NULL){ /* look for .bz2 file */
            strcpy(stimfilz,stimfile);
            strncat(stimfilz,".xz",STFILNAM-5);
            compress_type = 3;
            if ((ftemp=fopen(stimfilz,"r"))==NULL){ /* look for .xz file */
               ncfprintf (stderr,"%s: can't find stim file.\n",progname);
               stimfflg = 0;		/* if not, don't use file */
               return;
	    }
	 }
      }
      else {				/* .gz file exists */

        fclose(ftemp);
	switch (compress_type) {
		case 1: strcpy(stimfilz,"gzip -dc ");		/* check for .gz file */
			break;
		case 2: strcpy(stimfilz,"bzip2 -dc ");		/* check for .bz2 file */
			break;
		case 3: strcpy(stimfilz,"xz -dc ");		/* check for .xz file */
			break;
	}
        strncat(stimfilz,stimfile,STFILNAM-8);
	switch (compress_type) {
		case 1: strncat(stimfilz,".gz",STFILNAM-8);	/* .gz file */
			break;
		case 2: strncat(stimfilz,".bz2",STFILNAM-8);	/* .bz2 file */
			break;
		case 3: strncat(stimfilz,".xz",STFILNAM-8);	/* .xz file */
			break;
	}
        if ((ftemp=popen(stimfilz,"r"))==NULL) {
          ncfprintf (stderr,"%s: can't open file '%s'\n",progname,stimfile);
          stimfflg = 0;
          return;
        }
      }
    }
    stimin = ftemp;
    stimfflg = 1;
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#define BACKARR 0
#define STIMARR 1

extern int arraymade;				/* in modcode.cc */

void stim_backgr (double backgr, double wavel, double mask, int stimchan, double start)

{
    const char *action;
    photrec *rpnt;
    photorec *ept;
    int backarr, stimarr;

  if (makestim) {
      if (arraymade) {
          recback(BACKARR,backgr,wavel); /* set photrec backgr inten */
          recback(STIMARR,0.0,wavel);   /* zero photrec stim inten */
          abslist(0.0, start, mask, stimchan);  /* make action list, ratio 0.0 sets to backgr */
      } else {
	  execerror ("Please run stim_blur() before stim_backgr(),",
			"can't continue, stopping...");
          _exit(2);
      }
      return;
  }

  if (mask) action = "c";	/* c -> mask, a -> additive */
  else      action = "a";
	/* look for existing photoreceptors: */
  if (runyet) {
    for (rpnt=recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
        if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR ||
	      rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER)
        makrstim (start,rpnt->recnm1,rpnt->recnm2,rpnt->recnm3,
			rpnt->recnm4,backgr,wavel,mask,stimchan,0,action);
    }
  }
	/* look for new rod, code elements, */
        /*  even if run yet: */

  for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
    if (ept->ctype == ROD || ept->ctype == CONE || ept->ctype == CHR ||
        ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER)
    makrstim (start,ept->node1a,ept->node1b,ept->node1c,
		ept->node1d,backgr,wavel,mask,stimchan,0,action);
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_backgr (double backgr, int stimchan, double start)

{
  stim_backgr(backgr,1,0, stimchan, start);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_backgr (double backgr, double start)

{
  stim_backgr(backgr,1,0,0,start);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_backgr (double backgr, int stimchan)

{
  stim_backgr(backgr,1,0,stimchan,simtime);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_backgr (double backgr)

{
  stim_backgr(backgr,1,0,0,simtime);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

static int stim_centerfl = 0;			/* center specified */
static double stim_xcent=0, stim_ycent=0;	/* center of photoreceptor array */
extern int blurfl;				/* blur array made */
extern int cumstim;
extern int blursize;
extern int stimarrsize;				/* in modcode.cc */

#define max(x, y)	(((x) < (y)) ? (y) : (x))
#define ARREDGEMUL 1.0

void stim_blur (double blurrad, double blur_ampl, double scatter_ampl, double scatter_rad, 
			double scatter_pow, double sscale)

{
    extern int blursize;
    double xrange, yrange;
    double xmax,xmin,ymax,ymin;

   if (! makestim) return;

   blursize=makblur(blurrad,blur_ampl,sscale,  /* make a blur function */
        scatter_rad, scatter_pow, scatter_ampl);
   blurfl = 1;
   ncfprintf (stimout,"## %d blur     array size %d\n", ++cumstim, blursize*2); /* */
   if (arraymade) {                     /* if backgr array made */
       delc(2);                       /* delete existing arrays */
       arraymade = 0;
   }
   if (!arraymade) {
      find_photarrsiz(&xmax,&xmin,&ymax,&ymin);  /* find maxmin of photrecs */
//      if (!stim_centerfl) {
         stim_xcent = (xmax+xmin) * 0.5;
         stim_ycent = (ymax+ymin) * 0.5;
//      }
      xrange = xmax-xmin;
      yrange = ymax-ymin;
      if (sscale <=0) sscale = 1.0;
      stimarrsize = (long int)(max(xrange,yrange) * ARREDGEMUL / sscale);
      if (stimarrsize==0) stimarrsize = 1;
      if (stimarrsize&1) stimarrsize++;
      if (blurrad>0) stimarrsize += blursize*2;

      ncfprintf (stimout,"## %d stimulus array size %d\n",
                           ++cumstim,stimarrsize); /* */

      if (initc(2,stimarrsize)<=0) {           /* make stimulus to fit */
          ncfprintf (stderr,"Stimulus array is too big = %d\n",stimarrsize);
          ncfprintf (stderr,"Blur array dia %d\n",blursize*2);
          ncfprintf (stderr,"Stimulus xmax %g xmin %g ymax %g ymin %g\n",
                   xmax,xmin,ymax,ymin);
          execerror ("can't continue"," stopping...");
      }
      arraymade = 1;
   }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_blur (double blurrad, double scatter_ampl, double scatter_rad, 
			double scatter_pow, double sscale)
{
     double blur_ampl;

    stim_blur (blurrad, blur_ampl=1.0, scatter_ampl, scatter_rad, scatter_pow, sscale);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void make_barx (double width, double height, double xloc, double yloc, double orient,
	  double inten, double start, double dur, double wavel, double mask, int stimchan)
{
    double halfx, halfy;
    const char *action;
    photrec *rpnt;
    photorec *ept;

  if (stimfflg) return;

  if (mask>0) action = "c";	/* c -> abs mask, d -> additive */
  else        action = "d";
  halfx = width / 2;
  halfy = height / 2;
  if (runyet) {
    for (rpnt=recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
      if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR ||
	      rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER)
       if (inrect(rpnt->xloc,rpnt->yloc, xloc+halfx, xloc-halfx,
		yloc+halfy, yloc-halfy, orient)) {
	makrstim(start,dur,rpnt,inten,wavel,mask,stimchan,action);
      }
    }
  }
  else {
    for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
      if (ept->ctype == ROD || ept->ctype == CONE || ept->ctype == CHR ||
	      ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER)
      if (inrect(ept->xpos,ept->ypos,
		xloc+halfx, xloc-halfx,
		yloc+halfy, yloc-halfy, orient)) {
	makrstim(start,dur,ept,inten,wavel,mask,stimchan,action);
      }
    }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_bar (double width, double length, double xloc, double yloc, 
		double xcent, double ycent, double scale, double orient, 
		double inten, double start, double dur, double wavel, double mask, int stimchan)
{
   int arr;

   if (makestim) 
   makbar(arr=STIMARR, stimarrsize, width, length, xloc, yloc, xcent, ycent, scale,
				orient, inten, start, dur, wavel, mask, stimchan);
   else make_barx (width, length, xloc, yloc, orient, inten, start, dur, wavel, mask, stimchan);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_bar (double width, double height, double xloc, double yloc, double orient,
	  double inten, double start, double dur, double wavel, double mask, int stimchan)
{
   double xcent, ycent, scale;

   stim_bar (width, height, xloc, yloc, xcent=0, ycent=0, scale=1.0,
			orient, inten, start, dur, wavel, mask, stimchan);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_bar (double width, double height, double xloc, double yloc, double orient,
	  double inten, double start, double dur, double mask, int stimchan)
{
   double xcent, ycent, scale;

   stim_bar (width, height, xloc, yloc, xcent=0, ycent=0, scale=1.0,
			orient, inten, start, dur, 1.0, mask, stimchan);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_bar (double width, double height, double xloc, double yloc, double orient,
	  double inten, double start, double dur, double mask)
{
   double xcent, ycent, scale;

   stim_bar (width, height, xloc, yloc, xcent=0, ycent=0, scale=1.0,
			orient, inten, start, dur, 1.0, mask, 0);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_bar (double width, double height, double xloc, double yloc, double orient,
	  double inten, double start, double dur)
{
   double xcent, ycent, scale;

   stim_bar (width, height, xloc, yloc, xcent=0, ycent=0, scale=1.0,
			orient, inten, start, dur, 1.0, 0, 0);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void make_spotx (double dia, double xloc, double yloc, double inten, 
			double start, double dur, double wavel, double mask, int stimchan, int invert)

/* invert makes spot into inverted spot: useful for masking to central spot */

{
    double radius;
    const char *action;
    photrec *rpnt;
    photorec *ept;

  if (stimfflg) return;

  if (mask>0) action = "c";	/* c -> abs mask, d -> additive */
  else        action = "d";
  radius = dia / 2;
  if (runyet) {			/* look for existing photophotrecs */
    for (rpnt=(photrec*)recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
      if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR ||
	      rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER) {
        if (invert ^ incirc(rpnt->xloc,rpnt->yloc, xloc, yloc,radius)) {
	  makrstim(start,dur,rpnt,inten,wavel,mask,stimchan,action);
        }
      }
    }
  }
  else {
    for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
       if (ept->ctype == ROD || ept->ctype == CONE || ept->ctype == CHR ||
		      ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER) {
	if (invert ^ incirc(ept->xpos,ept->ypos, xloc, yloc,radius)) {
		makrstim(start,dur,ept,inten,wavel,mask,stimchan,action);
        }
      }
    }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_spot (double dia, double xloc, double yloc, double xcent, double ycent,
		double scale, double inten, double start, double dur, double wavel, double mask, int stimchan)
{
   int arr, invert;

   if (makestim) makspot(arr=STIMARR, stimarrsize, dia, xloc, yloc, xcent, ycent, 
		   	scale, inten, start, dur, wavel, mask, stimchan, invert=0);
   else make_spotx (dia, xloc, yloc, inten, start, dur, wavel, mask, stimchan, invert=0);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_spot (double dia, double xloc, double yloc, double inten, double start, 
		double dur, double wavel, double mask, int stimchan)
{
   double xcent, ycent, scale;

   stim_spot(dia, xloc, yloc, xcent=0, ycent=0, scale=1.0, inten, start, dur, wavel, mask, stimchan);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_spot (double dia, double xloc, double yloc, double inten, double start, 
		double dur, double wavel, double mask)
{
   double xcent, ycent, scale;
   int stimchan;

   stim_spot(dia, xloc, yloc, xcent=0, ycent=0, scale=1.0, inten, start, dur, wavel, mask, stimchan=0);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_spot (double dia, double xloc, double yloc, double inten, double start, double dur, double mask)
{
   stim_spot (dia, xloc, yloc, inten, start, dur, 1.0, mask);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_spot (double dia, double xloc, double yloc, double inten, double start, double dur)
{
   stim_spot (dia, xloc, yloc, inten, start, dur, 1.0, 0);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_ispot (double dia, double xloc, double yloc, double xcent, double ycent,
		double scale, double inten, double start, double dur, double wavel, 
		double mask, int stimchan, int invert)
{
   int arr;

   if (makestim) makspot(arr=STIMARR, stimarrsize, dia, xloc, yloc, xcent, ycent, 
		   	scale, inten, start, dur, wavel, mask, stimchan, invert);
   else make_spotx (dia, xloc, yloc, inten, start, dur, wavel, mask, stimchan, invert);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_ispot (double dia, double xloc, double yloc, double inten, double start, 
		double dur, double wavel, double mask, int stimchan, int invert)
{
   double xcent, ycent, scale;

   stim_ispot(dia, xloc, yloc, xcent=0, ycent=0, scale=1.0, inten, start, dur, wavel, mask, stimchan, invert);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_annulus (double idia, double odia, double xloc, double yloc, 
		   double inten, double start, double dur, double wavel, double mask)
{
   int arr, invert, stimchan;
   double xcent, ycent, scale;

 if (makestim) {
     makspot(arr=STIMARR, stimarrsize, odia, xloc, yloc, xcent=0, ycent=0, 
	   	scale=1.0,  inten, start, dur, wavel, mask, stimchan=0, invert=0);
     makspot(arr=STIMARR, stimarrsize, idia, xloc, yloc, xcent=0, ycent=0, 
	   	scale=1.0, -inten, start, dur, wavel, mask, stimchan=0, invert=0);
 }
 else {
   make_spotx (odia, xloc, yloc,  inten, start, dur, wavel, mask, stimchan=0, invert=0);
   make_spotx (idia, xloc, yloc, -inten, start, dur, wavel, mask, stimchan=0, invert=0);
 }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_annulus (double idia, double odia, double xloc, double yloc, 
		   double inten, double start, double dur)
{
     double wavel, mask;

  stim_annulus (idia, odia, xloc, yloc, inten, start, dur, wavel=1.0, mask=0);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void stim_image (char *imgfil, double width, double height, double xloc, double yloc, 
		double orient, double inten, double start, double dur)

{
    double val;
    double halfx, halfy;
    const char *action;
    double wavel=1, mask=0;
    photrec *rpnt;
    photorec *ept;
    FILE *fimg;
   
  if ((fimg=fopen(imgfil,"r"))==NULL) { /* open file */
    ncfprintf (stderr,"%s: can't open image file '%s'\n", progname,imgfil);
  }

  halfx = width / 2;
  halfy = height / 2;
  if (mask>0) action = "c";	/* c -> abs mask, d -> additive */
  else        action = "d";
  if (runyet) {
    for (rpnt=(photrec*)recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
      if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR ||
	      rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER)
       if (inrect(rpnt->xloc,rpnt->yloc, xloc+halfx, xloc-halfx,
		yloc+halfy, yloc-halfy, orient)) {
        fread(&val,4,1,fimg);
	makrstim(start,dur,rpnt, inten+val,wavel,mask,0,action);
       }
    }
  }
  else {
    for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
      if (ept->ctype == ROD || ept->ctype == CONE || ept->ctype == CHR ||
	      ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER)
      if (inrect(ept->xpos,ept->ypos, xloc+halfx, xloc-halfx,
		yloc+halfy, yloc-halfy, orient)) {
        fread(&val,4,1,fimg);
	makrstim(start,dur,ept,inten+val,wavel,mask,0,action);
      }
    }
  }
}

/*---------------------------------------------------------------------*/

double sinsq (double theta, int sq)

/* change a sine wave into a square wave */

{
    double val;

  if (sq==0) val = sin(theta);
  else {
    while (theta >= 2*MPI) theta -= (2*MPI);
    while (theta < 0)      theta += (2*MPI);
    if ((theta>0) && (theta<=MPI)) val = 1; 
    else                           val = -1;
  } 
  return val;
}



/*---------------------------------------------------------------------*/

double gaussf1(double x, double xcent, double rad)

	/* Return a gaussian function of x,
	 *    centered in an array of given size. */

{
	   double distsq,tempx;

	      tempx = (double)(x - xcent);
	         distsq = tempx*tempx;
		    return (exp(-(distsq) / (rad*rad)));

}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double gabormod(double contrast, double inten,
		double x, double y,
		double xcent, double ycent,
		double xenv,  double yenv, int sq, 
		double sino, double coso,
		double sphasenorm, double speriod, 
		double sinttheta, double stoff)

/* Low level function to generate drifting, counterphasing, or */
/* static sine-wave and gabor gratings. */

{
    double sinten, sinstheta, stheta;
    double xr,yr;

  xr = (x-xcent)*coso - (y-ycent)*sino + xcent; 
  yr = (x-xcent)*sino + (y-ycent)*coso + ycent; 
  stheta = (xr-xcent)/speriod + stoff + sphasenorm;
  sinstheta = sinsq(2*MPI * stheta,sq) * gaussf1(xr,xcent,xenv) *
                                    gaussf1(yr,ycent,yenv);

  sinten = inten * (contrast * sinstheta * sinttheta);
  return sinten;
}

/*------------------------------------------------------*/

double sineannmod(double contrast, double inten,
		double x, double y,
		double xcent, double ycent,
		double renv, int sq,
		double sphasenorm, double speriod, 
		double sinttheta, double stoff)

/* Low level function to generate drifting, counterphasing, or */
/* static concentric sine-wave and gabor gratings. */

{
    double sinten, sinstheta, stheta;
    double r,xtemp,ytemp;

  xtemp = x - xcent;
  ytemp = y - ycent;
  r = sqrt(xtemp*xtemp + ytemp*ytemp);
  stheta = r/speriod + stoff + sphasenorm;
  sinstheta = sinsq(2*MPI * stheta,sq) * gaussf1(r,0,renv);
  sinten = inten * (contrast * sinstheta * sinttheta);
  return sinten;
}

/*------------------------------------------------------*/

double sineannmod(double contrast, double inten,
		double x, double y,
		double xcent, double ycent,
		double renv,  
		double sphasenorm, double speriod, 
		double sinttheta, double stoff)
{
   return  sineannmod(contrast, inten, x, y, xcent, ycent, renv, 0,
		sphasenorm, speriod, sinttheta, stoff);
}


/*------------------------------------------------------*/

double windmillmod(double contrast, double inten,
		double x, double y,
		double xcent, double ycent,
		double renv, int sq,
		double sphasenorm, double speriod, 
		double sinttheta, double stoff)

/* Low level function to generate drifting, counterphasing, or */
/* static concentric windmill gratings. "speriod" sets the number */
/*  of vanes */

{
    double sinten, sinstheta, stheta;
    double r,t,xtemp,ytemp;

  xtemp = x - xcent;
  ytemp = y - ycent;
  r = sqrt(xtemp*xtemp + ytemp*ytemp);
  t = axy(xtemp,ytemp); 
  stheta = t*speriod/(2*MPI) + stoff + sphasenorm;
  sinstheta = sinsq(2*MPI * stheta,sq) * gaussf1(r,0,renv);
  sinten = inten * (contrast * sinstheta * sinttheta);
  return sinten;
}

/*------------------------------------------------------*/

double windmillmod(double contrast, double inten,
		double x, double y,
		double xcent, double ycent,
		double renv,  
		double sphasenorm, double speriod, 
		double sinttheta, double stoff)
{
return windmillmod(contrast, inten, x, y, xcent, ycent, renv, 0,
		sphasenorm, speriod, sinttheta, stoff);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#define STIMINC dsintinc        /* Min time incr for sine waves (.002 sec) */
#define STIMRES dsintres        /* Time res for sine waves (0.002) */

void makgaborx(double speriod, double sphase, double orient,
		double xloc, double yloc,
                double tfreq, double drift,
                double inten_add, double inten_mult, double contrast, 
		double wavel, double xenv, double yenv, 
		int makenv, int sq, double mask, int stimchan,
		const char *action,
		double start, double dur)

/* Make a Gabor filter (sine wave inside Gaussian envelope)
   in the stimulus list.
*/

{
  int tmod;
  double x, y;
  double xr, yr, xcent, ycent;
  double xedgemax, xedgemin, yedgemax, yedgemin;
  double coso, sino, orientrad, stheta, ttheta, sinttheta, stoff;
  double toffset,tperiod,tincr,stime,sinten,sphasenorm;
  double stimend;
  photorec *ept;
  photrec *rpnt;


  if (stimfflg) return;

  xcent = xloc;				/* yloc is xmax; xloc is xmin */
  ycent = yloc; 
  if (makenv) {				/* make envelope function */
    xedgemax = xcent + xenv * 5;
    xedgemin = xcent - xenv * 5;
    yedgemax = ycent + yenv * 5;
    yedgemin = ycent - yenv * 5;
  }
  else {  				/* sine, make sharp border */
    xedgemax = xcent + xenv;		/* x,yenv is original sharp border */
    xedgemin = xcent - xenv;
    yedgemax = ycent + yenv;
    yedgemin = ycent - yenv;
    xenv = 1e30; 
    yenv = 1e30; 
  }

  if (speriod <= 0.0) speriod = 1;	/* spatial period from user */

  orientrad = MPI / 180.0 * orient;
  sphasenorm = sphase / 360;
  coso = cos(orientrad);
  sino = sin(orientrad);
  if (drift || tfreq) {			/* drifting or counterphase grating */
     tmod = 1;				/* temporal modulation */
     if (tfreq==0.0) tfreq = 1.0;
     tperiod = 1.0 / tfreq;		/* temporal period */
     tincr = STIMRES * tperiod;		/* time incr = time res * t period */
     tincr = abs(tincr);		/* allow negative temporal frequency */
     if (tincr < STIMINC) tincr = STIMINC; /* time incr = photrec timinc */
     stimend = start+dur;
  }
  else {				/* static grating */
     tmod = 0;				/* no temporal modulation */
     stoff = 0.0;			/* no spatio-temporal offset */
     sinttheta = 1.0;			/* no counterphase modulation */
     stimend = start+dur;		/* stop after one step */
     tincr = dur;			/* dummy tincr */
  }
  for (stime=start; stime<stimend; stime += tincr) {

   if (stime+tincr >= stimend) {
     tincr = stimend-stime; 		/* fix last step so it doesn't overflow */
   }
   if (tmod) {
     toffset = (start-stime)/tperiod;
     if (drift) {                      /* drifting grating */
       stoff = drift * toffset;        /* add temporal shift into spatial phase */
       sinttheta = 1.0;
     }
     else      {                       /* stationary counterphasing grating */
       stoff = 0.0;
       sinttheta = sinsq(2*MPI * toffset,sq); /* temporal shift is indep sine func */
     }
   }
   if (runyet) {
    for (rpnt=recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
      if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR ||
	      rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER) {
	 x = rpnt->xloc;
	 y = rpnt->yloc;

			/* if the photorec is under grating */
          if (inrect(x,y,xedgemax,xedgemin,yedgemax,yedgemin,orient)) {

            sinten = inten_add + gabormod(contrast, inten_mult, x, y, xcent, ycent, xenv, yenv, 
			sq, sino, coso, sphasenorm, speriod, sinttheta, stoff);
	  } 
	  else {		/* set mean intens if outside */
            sinten = inten_add;	
	  }
	  makrstim(stime,tincr,	rpnt,sinten,wavel,mask,stimchan,action);

       } /* if (rpnt->ctype==ROD... */
     } /* for (rpnt;;) */
  } else {
      for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
        if (ept->ctype == ROD || ept->ctype == CONE || ept->ctype == CHR ||
	      ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER) {
	  x = ept->xpos;
	  y = ept->ypos;

			/* if the photorec is under grating */
          if (inrect(x,y,xedgemax,xedgemin,yedgemax,yedgemin,orient)) {

              sinten = inten_add + gabormod(contrast, inten_mult, x, y, xcent,ycent, xenv, yenv, 
			sq, sino, coso, sphasenorm, speriod, sinttheta, stoff);
	   } 
	     else {		/* set mean intens if outside */
              sinten = inten_add;	
	   }
	   makrstim(stime,tincr,ept,sinten,wavel,mask,stimchan,action);

        } /* if (ept->ctype==ROD... */
      } /* for (ept;;)   */
    }  /* else (!runyet) */
  }  /* for (stime;;) */
} 

/*------------------------------------------------------*/

void makgaborx(double speriod, double sphase, double orient,
		double xloc, double yloc,
                double tfreq, double drift,
                double inten_mult, double contrast, 
		double wavel, double xenv, double yenv, 
		int makenv, int sq, double mask, int stimchan, 
		const char *action,
		double start, double dur)

{
     double inten_add;

   makgaborx(speriod, sphase, orient, xloc, yloc, tfreq, drift, inten_add=0, inten_mult, contrast, 
		   wavel, xenv, yenv, makenv, sq, mask, stimchan, action, start, dur);
}

/*------------------------------------------------------*/

void maksineannx(double speriod, double sphase,
		double xloc, double yloc,
                double tfreq, double drift,
                double inten_add, double inten_mult, double contrast, 
		double wavel, double renv, 
		int makenv, int sq, double mask, int stimchan,
		const char *action,
		double start, double dur)

/* Make a sine wave annulus (radial sine wave inside Gaussian envelope)
   in the stimulus list.
*/

{
  int tmod;
  double x, y;
  double r, xcent, ycent;
  double xtemp, ytemp, edgemax;
  double coso, sino, stheta, ttheta, sinttheta, stoff;
  double toffset,tperiod,tincr,stime,sinten,sphasenorm;
  double stimend;
  photorec *ept;
  photrec *rpnt;


  if (stimfflg) return;

  sphasenorm = sphase / 360;
  xcent = xloc;				/* yloc is xmax; xloc is xmin */
  ycent = yloc; 
  if (makenv) {				/* make envelope function */
    edgemax = renv * 5;
  }
  else {  				/* sine, make sharp border */
    edgemax = renv;
    renv = 1e30; 
  } if (speriod <= 0.0) speriod = 1;	/* spatial period from user */
  if (drift || tfreq) {			/* drifting or counterphase grating */
     tmod = 1;				/* temporal modulation */
     if (tfreq==0.0) tfreq = 1.0;
     tperiod = 1.0 / tfreq;		/* temporal period */
     tincr = STIMRES * tperiod;		/* time incr = time res * t period */
     tincr = abs(tincr);		/* allow negative temporal frequency */
     if (tincr < STIMINC) tincr = STIMINC; /* time incr = photrec timinc */
     stimend = start+dur;
  }
  else {				/* static grating */
     tmod = 0;				/* no temporal offset */
     stoff = 0.0;			/* no spatio-temporal offset */
     sinttheta = 1.0;			/* no temporal modulation */
     stimend = start+dur;		/* stop after one step */
     tincr = dur;			/* dummy tincr */
  }
  for (stime=start; stime<stimend; stime += tincr) {

   if (stime+tincr >= stimend) {
     tincr = stimend-stime; 		/* fix last step so it doesn't overflow */
   }
   if (tmod) {
     toffset = (start-stime)/tperiod;
     if (drift) {                      /* drifting grating */
       stoff = drift * toffset;        /* add temporal shift to spatial phase */
       sinttheta = 1.0;
     }
     else      {                       /* stationary counterphasing grating */
       stoff = 0.0;
       sinttheta = sin(2*MPI * toffset); /* temporal shift is indep sine func */
     }
   }
   if (runyet) {
    for (rpnt=recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
     if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR ||
	      rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER) {
	 x = rpnt->xloc;
	 y = rpnt->yloc;

			/* if the photorec is under grating */
         if (incirc(x,y,xcent,ycent,edgemax)) {

            sinten = inten_add + sineannmod(contrast, inten_mult, x, y, xcent, ycent, renv,  
			sq, sphasenorm, speriod, sinttheta, stoff);

	 } 
	 else {		/* set mean intens if outside */
            sinten = inten_add;	
	 }
	 makrstim(stime,tincr,rpnt,sinten,wavel,mask,stimchan,action);

       } /* if (rpnt->ctype==ROD... */
     } /* for (rpnt;;) */
  } else {
     for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
       if (ept->ctype == ROD || ept->ctype == CONE || ept->ctype == CHR ||
	      ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER) {
	  x = ept->xpos;
	  y = ept->ypos;

			/* if the photorec is under grating */
          if (incirc(x,y,xcent,ycent,edgemax)) {

             sinten = inten_add + sineannmod(contrast, inten_mult, x, y, xcent, ycent, renv,  
			sq, sphasenorm, speriod, sinttheta, stoff);
	  } 
	  else {		/* set mean intens if outside */
             sinten = inten_add;	
	  }
	  makrstim(stime,tincr,ept,sinten,wavel,mask,stimchan,action);

        } /* if (ept->ctype==ROD... */
      } /* for (ept;;)   */
    }  /* else (!runyet) */
  }  /* for (stime;;) */
} 

/*------------------------------------------------------*/

void maksineannx(double speriod, double sphase,
		double xloc, double yloc,
                double tfreq, double drift,
                double inten_mult, double contrast, 
		double wavel, double renv, 
		int makenv, int sq, double mask, int stimchan,
		const char *action,
		double start, double dur)
{
    double inten_add;

    maksineannx(speriod, sphase, xloc, yloc, tfreq, drift, inten_add=0, inten_mult, contrast, wavel, 
		    renv, makenv, sq, mask, stimchan, action, start, dur);

}

/*------------------------------------------------------*/

void makwindmillx(double speriod, double sphase,
		double xloc, double yloc,
                double tfreq, double drift,
                double inten_add, double inten_mult, double contrast, 
		double wavel, double renv, 
		int makenv, int sq, double mask, int stimchan,
		const char *action,
		double start, double dur)

/* Make a sine wave windmill in the stimulus list.  */

{
  int tmod;
  double x, y;
  double r, xcent, ycent;
  double xtemp, ytemp, edgemax;
  double coso, sino, stheta, ttheta, sinttheta, stoff;
  double toffset,tperiod,tincr,stime,sinten,sphasenorm;
  double stimend;
  photorec *ept;
  photrec *rpnt;


  if (stimfflg) return;

  sphasenorm = sphase / 360;
  xcent = xloc;				/* yloc is xmax; xloc is xmin */
  ycent = yloc; 
  if (makenv) {				/* make envelope function */
    edgemax = renv * 5;
  }
  else {  				/* sine, make sharp border */
    edgemax = renv;
    renv = 1e30; 
  } if (speriod <= 0.0) speriod = 1;	/* spatial period from user */
  if (drift || tfreq) {			/* drifting or counterphase grating */
     tmod = 1;				/* temporal modulation */
     if (tfreq==0.0) tfreq = 1.0;
     tperiod = 1.0 / tfreq;		/* temporal period */
     tincr = STIMRES * tperiod;		/* time incr = time res * t period */
     tincr = abs(tincr);		/* allow negative temporal frequency */
     if (tincr < STIMINC) tincr = STIMINC; /* time incr = photrec timinc */
     stimend = start+dur;
  }
  else {				/* static grating */
     tmod = 0;				/* no temporal offset */
     stoff = 0.0;			/* no spatio-temporal offset */
     sinttheta = 1.0;			/* no temporal modulation */
     stimend = start+dur;		/* stop after one step */
     tincr = dur;			/* dummy tincr */
  }
  for (stime=start; stime<stimend; stime += tincr) {

   if (stime+tincr >= stimend) {
     tincr = stimend-stime; 		/* fix last step so it doesn't overflow */
   }
   if (tmod) {
     toffset = (start-stime)/tperiod;
     if (drift) {                      /* drifting grating */
       stoff = -drift * toffset;        /* add temporal shift to spatial phase */
       sinttheta = 1.0;
     }
     else      {                       /* stationary counterphasing grating */
       stoff = 0.0;
       sinttheta = sin(2*MPI * toffset); /* temporal shift is indep sine func */
     }
   }
   if (runyet) {
    for (rpnt=recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
     if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR ||
	      rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER) {
	 x = rpnt->xloc;
	 y = rpnt->yloc;

			/* if the photorec is under grating */
         if (incirc(x,y,xcent,ycent,edgemax)) {

            sinten = inten_add + windmillmod(contrast, inten_mult, x, y, xcent, ycent, renv, sq,
			sphasenorm, speriod, sinttheta, stoff);

	 } 
	 else {		/* set mean intens if outside */
            sinten = inten_add;	
	 }
	 makrstim(stime,tincr,rpnt,sinten,wavel,mask,stimchan,action);

       } /* if (rpnt->ctype==ROD... */
     } /* for (rpnt;;) */
  } else {
     for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
       if (ept->ctype == ROD || ept->ctype == CONE || ept->ctype == CHR ||
	      ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER) {
	  x = ept->xpos;
	  y = ept->ypos;

			/* if the photorec is under grating */
          if (incirc(x,y,xcent,ycent,edgemax)) {

             sinten = inten_add + windmillmod(contrast, inten_mult, x, y, xcent, ycent, renv, sq,
			sphasenorm, speriod, sinttheta, stoff);
	  } 
	  else {		/* set mean intens if outside */
             sinten = inten_add;	
	  }
	  makrstim(stime,tincr,ept,sinten,wavel,mask,stimchan,action);

        } /* if (ept->ctype==ROD... */
      } /* for (ept;;)   */
    }  /* else (!runyet) */
  }  /* for (stime;;) */
} 

/*------------------------------------------------------*/

void makwindmillx(double speriod, double sphase,
		double xloc, double yloc,
                double tfreq, double drift,
                double inten_mult, double contrast, 
		double wavel, double renv, 
		int makenv, int sq, double mask, int stimchan,
		const char *action,
		double start, double dur)
{
    double inten_add;

  makwindmillx(speriod, sphase, xloc, yloc, tfreq, drift, inten_add=0, inten_mult, contrast, wavel, 
		   renv, makenv, sq, mask, stimchan, action, start, dur);
}

/*------------------------------------------------------*/

void makcheckerboardx(double width, double height, int xn, int yn, double orient,
		double xloc, double yloc, double tfreq, double inten, double contrast,
		const char *action, double start, double dur, double **rndarr, int *nfr,int ckrseed)

/* width, height = overall size */
/* xn, yn = size in pixels */

{
  int i,j,k,indx,nframes,framesiz,kindx;
  double x,y,t;
  double pixwidth, pixheight;
  double tperiod, tincr, stimend;
  double orientrad, coso, sino;
  double rand_inten, stime, sinten, wavel, mask;
  double xedgemax, xedgemin, yedgemax, yedgemin;
  double xoff, yoff;
  double *rand_arr;
  int stimchan;
  char *chpnt;
  photorec *ept;
  photrec *rpnt;
  unsigned int chseed;

  if (stimfflg) return;

  if (tfreq==0.0) tfreq = 1.0;
  tincr = 1.0 / tfreq;               /* temporal period */
  tincr = abs(tincr);                /* allow negative temporal frequency */
  if (tincr < stiminc) tincr = stiminc; /* time incr = photrec timinc */
  nframes = dur / tincr + 0.5;
  stimend = start + tincr*nframes; /* calc stim end from nframes */
 
  pixwidth  = width  / xn;
  pixheight = height / yn;
  framesiz = xn * yn;

  /* make a new array to contain the spatio-temporal stimulus */

  if ((rand_arr = (double*)emalloc(xn * yn * nframes * sizeof(double)))==NULL) {
     ncfprintf (stderr,"nc: checkerboard, can't init array\n");
     return;
  }
  if (rndarr != NULL && nfr != NULL) {
     *rndarr = rand_arr;
     *nfr = nframes;
  }

  /* make a new random number generator and initialize it */

  chseed = ckrseed^7924517;
  initrand (1,chseed);

  /* load the array with the random stimulus (-1 -> 1) */

  for (k=indx=0; k<nframes; k++) {
    for (j=0; j<yn; j++) {  
      for (i=0; i<xn; i++) {
        rand_arr[indx++] = contrast * (rrand(1) * 2.0 - 1.0);
      }
    }
  }

  orientrad = MPI / 180.0 * orient;
  coso = cos(orientrad);
  sino = sin(orientrad);
  wavel = 0;
  mask = 0;
  stimchan = 0;

  for (j=indx=0; j<yn; j++) {  // each region has unique inten
    for (i=0; i<xn; i++,indx++) {
      x = i*pixwidth+pixwidth/2-width/2;
      y = j*pixheight+pixheight/2-height/2;
      xoff = coso*x - sino*y;
      yoff = sino*x + coso*y;
      xedgemin = xoff + xloc;
      xedgemax = xedgemin + pixwidth;
      yedgemin = yoff + yloc;
      yedgemax = yedgemin + pixheight;

      /* make a rectangle pixel */

      if (runyet) {

      /* for each photoreceptor inside the checkerboard pixel, set stimulus for all frames */

       for (rpnt=recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
         if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR ||
	         rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER) {
	   x = rpnt->xloc;
	   y = rpnt->yloc;

	   for (k=kindx=0; k<nframes; k++,kindx+=framesiz) {

	     stime = start + k * tincr;

			/* if the photorec is under grating */
             if (inrect(x,y,xedgemax,xedgemin,yedgemax,yedgemin,orient)) {
                 sinten = inten + rand_arr[kindx+indx];
	     } 
	     else {		/* set mean intens if outside */
                 sinten = inten;	
	     }

	     makrstim(stime,tincr,rpnt,sinten,wavel,mask,stimchan,action);
	   } /* for (k;;) */

         } /* if (rpnt->ctype==ROD... */
       } /* for (rpnt;;) */
  } else {
      for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
        if (ept->ctype == ROD || ept->ctype == CONE || ept->ctype == CHR ||
	      ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER) {
	  x = ept->xpos;
	  y = ept->ypos;

	   for (k=kindx=0; k<nframes; k++,kindx+=framesiz) {

	     stime = start + k * tincr;

			/* if the photorec is under grating */
             if (inrect(x,y,xedgemax,xedgemin,yedgemax,yedgemin,orient)) {
                sinten = inten + rand_arr[kindx+indx];
	     } 
	     else {		/* set mean intens if outside */
                sinten = inten;	
	     }
	     makrstim(stime,tincr,ept,sinten,wavel,mask,stimchan,action);
	  } /* for (k;;) */

        } /* if (ept->ctype==ROD... */
      } /* for (ept;;)   */
    }  /* else (!runyet) */

    } /* for (i<xn;)  */
  } /* for (j<yn;)  */

}


/*---------------------------------------------------------------------*/
 
void stim_sine(double speriod, double sphase, double orient, 
		double xloc, double yloc, double xcent, double ycent,
		double tfreq, int drift, double scale, 
		double inten_add, double inten_mult, double contrast, 
		int sq, double start, double dur, double wavel, double mask, int stimchan)
{
    const char *action;
    int arr,makenv;
    double xenv, yenv;

  if (mask>0) action = "c";	/* b -> abs mask, d -> additive */
  else        action = "d";

  xenv = 1e30;
  yenv = 1e30;
  if (makestim)
	makgabor(arr=STIMARR, stimarrsize, contrast, inten_add, inten_mult,
		xloc, yloc, xcent, ycent=0,
		xenv, yenv, makenv=0, sq,
		tfreq, drift, speriod, sphase,
		scale, orient, wavel,
		mask, stimchan, start, dur);
  else 
	makgaborx(speriod, sphase, orient, xloc, yloc, tfreq, drift,
                inten_add, inten_mult, contrast, wavel=1, xenv, yenv, makenv=0, sq,
		mask=0, stimchan, action="d", start, dur);
}

/*---------------------------------------------------------------------*/

void stim_sine(double speriod, double sphase, double orient, 
		double xloc, double yloc, double xcent, double ycent,
		double tfreq, int drift, double scale, 
		double inten_mult, double contrast, 
		int sq, double start, double dur, double wavel, double mask, int stimchan)
{
    double inten_add;

  stim_sine(speriod, sphase, orient, xloc, yloc, xcent, ycent, tfreq, drift, scale, 
		   inten_add=0, inten_mult, contrast, sq, start, dur, wavel, mask, stimchan); 
}
/*---------------------------------------------------------------------*/

void stim_sine(double speriod, double sphase, double orient, 
		double xloc, double yloc, double tfreq, int drift, 
		double inten_add, double inten_mult, double contrast, int sq, double start, double dur)

{
   double xcent, ycent, scale;
   double wavel, mask;
   int stimchan;

   stim_sine(speriod, sphase, orient, xloc, yloc, xcent=0, ycent=0,
		tfreq, drift, scale=1, inten_add, inten_mult, contrast, sq, start, dur, wavel=1, mask=0, stimchan=0);
}

/*---------------------------------------------------------------------*/

void stim_sine(double speriod, double sphase, double orient, 
		double xloc, double yloc, double tfreq, int drift, 
		double inten_mult, double contrast, int sq, double start, double dur)

{
   double xcent, ycent, scale;
   double wavel, mask, inten_add;
   int stimchan;

   stim_sine(speriod, sphase, orient, xloc, yloc, xcent=0, ycent=0,
		tfreq, drift, scale=1, inten_add=0, inten_mult, contrast, sq, start, dur, wavel=1, mask=0,stimchan=0);
}

/*---------------------------------------------------------------------*/
   
void stim_gabor(double speriod, double sphase, double orient, 
		double xloc, double yloc, double xcent, double ycent,
		double tfreq, double drift, double scale,
		double inten_add, double inten_mult, double contrast, double xenv, double yenv, 
		int sq, double start, double dur, double wavel, double mask, int stimchan)
{
    const char *action;
    int arr,makenv;

  if (mask>0) action = "c";	/* b -> abs mask, d -> additive */
  else        action = "d";

  if (makestim) 
	makgabor(arr=STIMARR, stimarrsize, contrast, inten_add, inten_mult,
		xloc, yloc, xcent=0, ycent=0,
		xenv, yenv, makenv=1, sq,
		tfreq, drift, speriod, sphase,
		scale, orient, wavel,
		mask, stimchan, start, dur);
  else
  	makgaborx(speriod, sphase, orient, xloc, yloc, tfreq, drift,
                inten_add, inten_mult, contrast, wavel=1, xenv, yenv, makenv=1, sq,
		mask=0, stimchan, action="d", start, dur);

}
/*---------------------------------------------------------------------*/

void stim_gabor(double speriod, double sphase, double orient, 
		double xloc, double yloc, double xcent, double ycent,
		double tfreq, double drift, double scale,
		double inten_mult, double contrast, double xenv, double yenv, 
		int sq, double start, double dur, double wavel, double mask, int stimchan)
{
     double inten_add;

   stim_gabor(speriod, sphase, orient, xloc,  yloc,  xcent,  ycent,  tfreq,  drift,  scale,  inten_add=0,  inten_mult,  contrast,  xenv,  yenv,  sq,  start,  dur,  wavel,  mask, stimchan);

}

/*---------------------------------------------------------------------*/

void stim_gabor(double speriod, double sphase, double orient, 
		double xloc, double yloc, 
		double tfreq, double drift, 
		double inten_add, double inten_mult, double contrast, 
		double xenv, double yenv, int sq, double start, double dur)
{
   double xcent, ycent, scale;
   double wavel, mask;
   int stimchan;

  stim_gabor(speriod, sphase, orient, xloc, yloc, xcent=0, ycent=0, tfreq, drift, scale=1,
		inten_add, inten_mult, contrast, xenv, yenv, sq, start, dur, wavel=1, mask=0, stimchan=0);
}

/*---------------------------------------------------------------------*/

void stim_gabor(double speriod, double sphase, double orient, 
		double xloc, double yloc, 
		double tfreq, double drift, 
		double inten_mult, double contrast, 
		double xenv, double yenv, int sq, double start, double dur)
{
   double xcent, ycent, scale;
   double wavel, mask, inten_add;
   int stimchan;

  stim_gabor(speriod, sphase, orient, xloc, yloc, xcent=0, ycent=0, tfreq, drift, scale=1,
		inten_add=0, inten_mult, contrast, xenv, yenv, sq, start, dur, wavel=1, mask=0, stimchan=0);
}

/*---------------------------------------------------------------------*/
// 20
void stim_sineann(double speriod, double sphase, double xloc, double yloc, 
		  double xcent, double ycent, double tfreq, double drift, double scale,
		  double inten_add, double inten_mult, double contrast, double renv, int makenv, int sq, 
		  double start, double dur, double wavel, double mask, int stimchan)
{
   int arr;
   const char *action;

  if (mask>0) action = "c";	/* b -> abs mask, d -> additive */
  else        action = "d";

  if (makestim)
	maksineann(arr=STIMARR, stimarrsize, contrast, inten_add, inten_mult, xloc, yloc, xcent, ycent,
		renv, makenv, sq, tfreq, drift, speriod, sphase,
		scale, wavel, mask, stimchan, start, dur);
  else
	maksineannx(speriod, sphase, xloc, yloc, tfreq, drift, inten_add, inten_mult,
		contrast, wavel=1, renv, makenv, sq, mask=0, stimchan, action, start, dur);
}

/*---------------------------------------------------------------------*/
// 19
void stim_sineann(double speriod, double sphase, double xloc, double yloc, 
		  double xcent, double ycent, double tfreq, double drift, double scale,
		  double inten_mult, double contrast, double renv, int makenv, int sq, 
		  double start, double dur, double wavel, double mask, int stimchan)
{
    double inten_add;

  stim_sineann(speriod, sphase, xloc, yloc, xcent, ycent, tfreq, drift, scale, 
		inten_add=0, inten_mult, contrast, renv, makenv, sq, start, dur, wavel, mask, stimchan);
}

/*---------------------------------------------------------------------*/
// 18
void stim_sineann(double speriod, double sphase, double xloc, double yloc, 
		  double xcent, double ycent, double tfreq, double drift, double scale,
		  double inten_mult, double contrast, double renv, int sq, 
		  double start, double dur, double wavel, double mask, int stimchan)
{
    int makenv;
    double inten_add;

  stim_sineann(speriod, sphase, xloc, yloc, xcent, ycent, tfreq, drift, scale, 
		inten_add=0, inten_mult, contrast, renv, makenv=1, sq, start, dur, wavel, mask, stimchan);
}

/*---------------------------------------------------------------------*/
// 12
void stim_sineann(double speriod, double sphase, double xloc, double yloc, 
		  double tfreq, double drift, double inten_add, double inten_mult, 
		  double contrast, double renv, double start, double dur)
{
   double xcent, ycent, scale;
   double wavel, mask;
   int stimchan;

   stim_sineann(speriod, sphase, xloc, yloc, xcent=0, ycent=0, tfreq, drift, scale=1,
		  inten_add, inten_mult, contrast, renv, 0,  start, dur, wavel=1, mask=0, stimchan=0);
}

/*---------------------------------------------------------------------*/
// 13
void stim_sineann(double speriod, double sphase, double xloc, double yloc, 
		  double tfreq, double drift, double inten_add, double inten_mult, 
		  double contrast, double renv, int sq, double start, double dur)
{
   double xcent, ycent, scale;
   double wavel, mask;
   int stimchan;

   stim_sineann(speriod, sphase, xloc, yloc, stim_xcent, stim_ycent, tfreq, drift, scale=1,
		  inten_add, inten_mult, contrast, renv, sq,  start, dur, wavel=1, mask=0, stimchan=0);
}

/*---------------------------------------------------------------------*/
// 14
/* makenv sets hard outer limit: */

void stim_sineann(double speriod, double sphase, double xloc, double yloc, 
		  double tfreq, double drift, double inten_add, double inten_mult, 
		  double contrast, double renv, int makenv, int sq, double start, double dur)
{
   double xcent, ycent, scale;
   double wavel, mask;
   int stimchan;

   stim_sineann(speriod, sphase, xloc, yloc, stim_xcent, stim_ycent, tfreq, drift, scale=1,
		  inten_add, inten_mult, contrast, renv, makenv, sq,  start, dur, wavel=1, mask=0, stimchan=0);
}

/*---------------------------------------------------------------------*/
// 11
/* in these, set inten_add = 0 : */

void stim_sineann(double speriod, double sphase, double xloc, double yloc, 
		  double tfreq, double drift, double inten_mult, 
		  double contrast, double renv, double start, double dur)
{
   double xcent, ycent, scale;
   double wavel, mask, inten_add;

   stim_sineann(speriod, sphase, xloc, yloc, stim_xcent, stim_ycent, tfreq, drift, scale=1,
		  inten_add=0, inten_mult, contrast, renv, 0,  start, dur, wavel=1, mask=0);
}

/*---------------------------------------------------------------------*/
// 12, int sq
void stim_sineann(double speriod, double sphase, double xloc, double yloc, 
		  double tfreq, double drift, double inten_mult, 
		  double contrast, double renv, int sq, double start, double dur)
{
   double xcent, ycent, scale;
   double wavel, mask, inten_add;

   stim_sineann(speriod, sphase, xloc, yloc, stim_xcent, stim_ycent, tfreq, drift, scale=1,
		  inten_add=0, inten_mult, contrast, renv, sq,  start, dur, wavel=1, mask=0);
}

/*---------------------------------------------------------------------*/
// 13, int makenv, sq
void stim_sineann(double speriod, double sphase, double xloc, double yloc, 
		  double tfreq, double drift, double inten_mult, 
		  double contrast, double renv, int makenv, int sq, double start, double dur)
{
   double xcent, ycent, scale;
   double wavel, mask, inten_add;
   int stimchan;

   stim_sineann(speriod, sphase, xloc, yloc, stim_xcent, stim_ycent, tfreq, drift, scale=1,
		  inten_add=0, inten_mult, contrast, renv, makenv, sq,  start, dur, wavel=1, mask=0, stimchan=0);
}

/*---------------------------------------------------------------------*/

void stim_windmill(double speriod, double sphase, double xloc, double yloc, 
			double xcent, double ycent, double tfreq, double drift, 
			double scale, double inten_add, double inten_mult, double contrast, 
			double renv, int makenv, int sq, double start, double dur, 
			double wavel, double mask, int stimchan)
{
   int arr;
   const char *action;

  if (mask>0) action = "b";	/* b -> mask, d -> additive */
  else        action = "d";

  if (makestim)
	makwindmill(arr=STIMARR, stimarrsize, contrast, inten_add, inten_mult, xloc, yloc, xcent, ycent,
		renv, makenv, sq, tfreq, drift, speriod, sphase, scale, 
		start, dur, wavel, mask, stimchan);
  else
	makwindmillx(speriod, sphase, xloc, yloc, tfreq, drift, inten_add, inten_mult, 
    		contrast, wavel=1, renv, makenv, sq, mask, stimchan, action, start, dur);
}

/*---------------------------------------------------------------------*/

void stim_windmill(double speriod, double sphase, double xloc, double yloc, 
			double xcent, double ycent, double tfreq, double drift, 
			double scale, double inten_mult, double contrast, 
			double renv, int sq, double start, double dur, 
			double wavel, double mask, int stimchan)
{
     double inten_add;
     int makenv;

  stim_windmill(speriod, sphase, xloc, yloc, xcent, ycent, tfreq, drift, scale, 
		inten_add=0, inten_mult, contrast, renv, makenv=1, sq, start, dur, wavel, mask, stimchan);
}

/*---------------------------------------------------------------------*/

void stim_windmill(double speriod, double sphase, double xloc, double yloc, 
			double tfreq, double drift, double inten_add, double inten_mult,
			double contrast, double renv, double start, double dur)
{
   double xcent, ycent, scale;
   double wavel, mask;
   int stimchan,makenv;

   stim_windmill(speriod, sphase, xloc, yloc, xcent=0, ycent=0, tfreq, drift, scale=1, 
			inten_add, inten_mult, contrast, renv, makenv=1,0, start, dur, wavel=1, mask=0, stimchan=0);
}

/*---------------------------------------------------------------------*/

void stim_windmill(double speriod, double sphase, double xloc, double yloc, 
			double tfreq, double drift, double inten_add, double inten_mult, 
			double contrast, double renv, int sq, double start, double dur)
{
   double xcent, ycent, scale;
   double wavel, mask;
   int stimchan,makenv;

   stim_windmill(speriod, sphase, xloc, yloc, xcent=0, ycent=0, tfreq, drift, scale=1, 
			inten_add, inten_mult, contrast, renv, makenv=1, sq, start, dur, wavel=1, mask=0, stimchan=0);
}

/*---------------------------------------------------------------------*/

void stim_windmill(double speriod, double sphase, double xloc, double yloc, 
			double tfreq, double drift, double inten_mult, double contrast, 
			double renv, double start, double dur)
{
   double xcent, ycent, scale;
   double wavel, mask, inten_add;
   int stimchan,makenv;

   stim_windmill(speriod, sphase, xloc, yloc, xcent=0, ycent=0, tfreq, drift, 
			scale=1, inten_add=0, inten_mult, contrast, renv, makenv=0, 0, start, dur, wavel=1, mask=0, stimchan=0);
}

/*---------------------------------------------------------------------*/

void stim_windmill(double speriod, double sphase, double xloc, double yloc, 
			double tfreq, double drift, double inten_mult, double contrast, 
			double renv, int sq, double start, double dur)
{
   double xcent, ycent, scale;
   double wavel, mask, inten_add;
   int stimchan,makenv;

   stim_windmill(speriod, sphase, xloc, yloc, xcent=0, ycent=0, tfreq, drift, 
			scale=1, inten_add=0, inten_mult, contrast, renv, makenv=1, sq, start, dur, wavel=1, mask=0, stimchan=0);
}

/*---------------------------------------------------------------------*/

void stim_checkerboard(double width, double height, int xn, int yn, 
		double orient, double xloc, double yloc, 
		double xcent, double ycent, double scale,
		double tfreq, double inten, double contrast, 
		double start, double dur, double **stim_rndarr, int *stim_nfr, int ckrseed)

{
   int arr;
   const char *action;

  if (makestim)
	makcheckerboard(arr=STIMARR, stimarrsize, width, height, xn, yn,
		xloc, yloc, xcent, ycent, scale, orient, tfreq,
		inten, contrast, start, dur, stim_rndarr, stim_nfr, ckrseed);
  else
	makcheckerboardx(width, height, xn, yn, orient, xloc, yloc,
		tfreq, inten, contrast, action="d", start, dur, stim_rndarr,stim_nfr, ckrseed);
}

/*---------------------------------------------------------------------*/

void stim_checkerboard(double width, double height, int xn, int yn, 
		double orient, double xloc, double yloc, 
		double tfreq, double inten, double contrast, 
		double start, double dur,double **stim_rndarr, int *stim_nfr, int ckrseed)
{
   double xcent, ycent, scale;

  stim_checkerboard(width, height, xn, yn, 
		orient, xloc, yloc, 
		xcent=0, ycent=0, scale=1,
		tfreq, inten, contrast, 
		start, dur, stim_rndarr, stim_nfr, ckrseed);
}

/*---------------------------------------------------------------------*/

void stim_checkerboard(double width, double height, int xn, int yn, 
		double orient, double xloc, double yloc, 
		double tfreq, double inten, double contrast,double start, double dur, int ckrseed)
{
   double xcent, ycent, scale;

  stim_checkerboard(width, height, xn, yn, 
		orient, xloc, yloc, 
		xcent=0, ycent=0, scale=1,
		tfreq, inten, contrast, 
		start, dur, NULL, NULL, ckrseed);
}

/*---------------------------------------------------------------------*/

void stim_checkerboard(double width, double height, int xn, int yn, 
		double orient, double xloc, double yloc, 
		double tfreq, double inten, double contrast,double start, double dur)
{
   double xcent, ycent, scale;

  stim_checkerboard(width, height, xn, yn, 
		orient, xloc, yloc, 
		xcent=0, ycent=0, scale=1,
		tfreq, inten, contrast, 
		start, dur, NULL, NULL, rseed);
}

/*---------------------------------------------------------------------*/

void stim_grating (int type, double speriod, double sphase, double orient,
	      double xloc, double yloc, double tfreq, double drift,
	      double inten_add, double inten_mult, double contrast, double wavel,
		double xenv, double yenv, int sq, double mask, int stimchan,
		double start, double dur)
{
     double xcent, ycent;
     double scale;
     int makenv;

  switch (type) {
	case SINE:			/* stim for possibly all photrecs */
	  stim_sine(speriod, sphase, orient, xloc, yloc, xcent=0, ycent=0,
		tfreq, drift, scale=1, inten_add, inten_mult, contrast, sq, start, dur,
	       	wavel, mask, stimchan);
	  break;

	case GABOR:			/* stim for possibly all photrecs */
	    stim_gabor(speriod, sphase, orient, xloc, yloc, xcent=0, ycent=0,
		tfreq, drift, scale=1, inten_add, inten_mult, contrast, xenv, yenv, sq,
		start, dur, wavel, mask, stimchan);
	  break;

	case SINEANN:			/* stim for possibly all photrecs */

	  stim_sineann(speriod, sphase, xloc, yloc, xcent=0, ycent=0, 
		tfreq, drift, scale=1, inten_add, inten_mult, contrast, xenv, sq,
		start, dur, wavel, mask, stimchan);
	  break;

	case WINDMILL:			/* stim for possibly all photrecs */

	  stim_windmill(speriod, sphase, xloc, yloc, xcent=0, ycent=0,
		tfreq, drift, scale=1, inten_add, inten_mult, contrast, xenv, makenv=1, sq, 
		start, dur, wavel, mask, stimchan);

	  break;
  }
}

/*---------------------------------------------------------------------*/

void stim_grating (int type, double speriod, double sphase, double orient,
	      double xloc, double yloc, double tfreq, double drift,
	      double inten, double contrast, double wavel,
		double xenv, double yenv, int sq, double mask, int stimchan,
		double start, double dur)
{
     double xcent, ycent;
     double scale;

  switch (type) {
	case SINE:			/* stim for possibly all photrecs */
	  stim_sine(speriod, sphase, orient, xloc, yloc, xcent=0, ycent=0,
		tfreq, drift, scale=1, inten, contrast, sq, start, dur,
	       	wavel, mask, stimchan);
	  break;

	case GABOR:			/* stim for possibly all photrecs */
	    stim_gabor(speriod, sphase, orient, xloc, yloc, xcent=0, ycent=0,
		tfreq, drift, scale=1, inten, contrast, xenv, yenv, sq,
		start, dur, wavel, mask, stimchan);
	  break;

	case SINEANN:			/* stim for possibly all photrecs */

	  stim_sineann(speriod, sphase, xloc, yloc, xcent=0, ycent=0, 
		tfreq, drift, scale=1, inten, contrast, xenv, sq,
		start, dur, wavel, mask, stimchan);
	  break;

	case WINDMILL:			/* stim for possibly all photrecs */

	  stim_windmill(speriod, sphase, xloc, yloc, xcent=0, ycent=0,
		tfreq, drift, scale=1, inten, contrast, xenv, sq, 
		start, dur, wavel, mask, stimchan);

	  break;
  }
}

/*---------------------------------------------------------------------*/

void sector_mask(double width, double orient, double xloc, double yloc, double val, double stimtime, double dur, int stimchan)

/* Make mask that masks everything but a pie-shaped wedge at a specified angle = orient */
/* Make wedge with 2 bars that overlap at different angles */

#define BARLENGTH 1000
{
  double w,h,xl,yl,th1,th2;
  double mask;

  th1 = orient - width/2;
  if (th1 <  0)   th1 += 360;
  if (th1 >= 360) th1 -= 360;
  if (width < 0)   width = -width;
  th2 = th1 + width + 180;
  if (th2 > 360) th2 -= 360;
  w = h = BARLENGTH;
  stim_bar (w, h, xl= -w*0.5*cos(th1/180.0*MPI)+xloc, yl=w*0.5*sin(th1/180.0*MPI)+yloc, th1, val, stimtime, dur, mask=1,stimchan);
  stim_bar (w, h, xl= -w*0.5*cos(th2/180.0*MPI)+xloc, yl=w*0.5*sin(th2/180.0*MPI)+yloc, th2, val, stimtime, dur, mask=1,stimchan);
}
#undef BARLENGTH

/*---------------------------------------------------------------------*/

void sector_mask(double width, double orient, double val, double stimtime, double dur)

{
   sector_mask(width, orient, 0, 0, val, stimtime, dur, 0);
}

