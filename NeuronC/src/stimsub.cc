/* segment stimsub in program stim */

/* subroutines to make stimulus and blur it by a gaussian.
 Prints array on standard output */

extern "C" {
#include <stdio.h>
#include <unistd.h>
}

#include "nc.h"
#include "y.tab.h"
#include "ndef.h"
#include "ncsub.h"
#include "ncelem.h"
#include "ncomp.h"
#include "stim.h"
#include "drand.h"
#include "control.h"
#include "wave.h"
#include "ncio.h"
#include "nconst.h"

#if (CONVSIZE <= 256)			/* if 16 bit compiler */
#define BLURSIZE 89
#else
#define BLURSIZE 256
#endif

extern int convsize;			/* size of stimulus array */

#define DEBUG
#ifdef DEBUG
#include "ncdebug.h"
#endif

#define SCALE 1
#define BACKTIME 0.00001

double *blurarr = 0;
int blursize=0;
static int blurmade=0;

double stimdia  = 10.0;
double stimydist = 10.0;
int stimhlp = 0;
extern int runyet;				/* dummy def for ncmain.c */
extern Symbol *timeptr;

static int blurzero = 0;		/* =1 -> no blur */

double swavel[2]= {0.0,0.0};

recnod *reclist = 0;                    /* head of recnod list */
recnod *reclend = 0;                    /* tail of recnod list */
extern int reccum;                      /* number of recnod's  */

extern FILE *stimout;

#ifdef __cplusplus
extern "C" {
#endif
  extern double atof(const char *);
  extern double sqrt(double);
  extern double exp(double);
  extern double sin(double);
  extern double cos(double);
  extern double pow(double,double);
#ifdef __cplusplus
}
#endif

void execerror (const char *str1, const char *str2);
double getblur(int x, int y);

double readc(int arr, int x, int y);
char *emalloc(unsigned int n);
char *smalloc(unsigned int n);
void setc(int arr, double val);
void delout (recnod *npt, double delinten, double mask, int stimchan, double time, int seq);
void absout (recnod *npt, double inten, double mask, int stimchan, double time, int seq);
int inrect(double x, double y, double maxx, double minx, 
		double maxy, double miny, double orient);
void addc(int arr, int x, int y, double val);
int incirc(double x, double y, double circx, double circy, double rad);
void efree (void *ptr);
double axy (double x, double y);
char *makrand (int rndsiz, const char *mesg, int num);

double gabormod(double contrast, double inten, double x, double y,
                double xcent, double ycent, double xenv,  double yenv,
                double sino, double coso, double sphasenorm, double speriod,
                double sinttheta, double toff);

double sineannmod(double contrast, double inten, double x, double y,
                double xcent, double ycent, double renv, double sphasenorm, 
		double speriod, double sinttheta, double toff);

double windmillmod(double contrast, double inten, double x, double y,
                double xcent, double ycent, double renv, double sphasenorm, 
		double speriod, double sinttheta, double toff);


/*------------------------------------*/

void printa (int arr)

/* debugging function */

{
  int x,y;
  
  ncfprintf (stdout,"# slice through center of sim array, along X-axis:\n");
  y = convsize / 2;  
  for (x=0; x<convsize; x++)
     ncfprintf (stdout,"# %-4d %-8g\n", x, readc (arr, x, y));
}


/*------------------------------------*/

void printb (void)

/* debugging function */

{
  int x,y;
  double getblur(int x, int y);
 
  if (!blurmade) {
    ncfprintf (stderr,"printb: blur array not made yet\n");
    return;
  }
  y = 0;  
  for (x=0; x<blursize; x++)
     ncfprintf (stdout,"%d %g\n", x, getblur(x,y));
}

/*---------------------------------------------------*/

int setstimchan (int stimchan)

    /* set global variable numstimchan from user set stimchan */
	        
{       
  if (stimchan < 0) stimchan = 0;
  if (stimchan >= MAXSTIMCHAN) stimchan = MAXSTIMCHAN-1;
  if (stimchan > numstimchan) numstimchan = stimchan;
  return stimchan;
}

/*------------------------------------*/

void recpat(int array, int arrsize, double xoff, double yoff, double scale)

/* Calculate stimulus and background intensity events
   for the list of photoreceptor stimulus nodes, given 
   the stimulus and background intensity arrays,
   and an optional (x,y) offset.

   Each stimulus node consists of:

typedef struct  RECNOD {
        nodeint recnm1;		receptor number
        nodeint recnm2;
        nodeint recnm3;
        nodeint recnm4;
        double xpos;		position, move the stim and backgnd arrs
        double ypos;
        double actual;		actual intensity for stimulus file
        double curval;		most recent value computed
        double wavel;		wavelength of stimulus
        double stim;		photoreceptor's stimulus value
        double backgnd;		photoreceptor's background value
        struct RECNOD *next;
        } recnod;

The "actual" intensity is only changed when steps of at least a
minimum size are generated.
 
*/
 
{
  int j,x,y,bmid,xmid,ymid,amid,interp,xo,yo,iscale;
  recnod *npt;
  double blurval, arrtot, getblur(int x, int y);
  double xr, yr, xr1, yr1, xoffs, yoffs;

#ifdef DEBUG
  if (debug & NCSTIM && debugz & 1)
		 ncfprintf (stderr,"recpat: start, arr %d\n",array);
#endif

  if (scale==0.0) scale = 1.0;
  iscale = int(1.0 / scale + 0.5);
  bmid = 0;					/* getblur has 0,0 center */
  amid = arrsize / 2;				/* for 0,0     center */
/*  amid = 0;					/* for 128,128 center */
  xoffs = xoff * iscale;
  yoffs = yoff * iscale;
  xr = xoffs - (int)(xoffs);
  yr = yoffs - (int)(yoffs);
  for (npt=reclist,j=0; npt; npt=npt->next) {
    xmid = (int)(npt->xpos * iscale + amid - bmid - xoffs + 0.5);
    ymid = (int)(npt->ypos * iscale + amid - bmid - yoffs + 0.5);
    arrtot = 0.0;
    if (blurzero) 				/* no blur */
       arrtot = readc(array, xmid, ymid);   
    else {
     if (!blurmade) {
       execerror ("recpat: blur array not made yet,", 
		  "can't continue, stopping...");
       _exit(2);
       return;
     }						/* run convolution */
      for (x= -(blursize-1); x<blursize; x++)
       for (y= -(blursize-1); y<blursize; y++) {
	 interp = (xr != 0) + ((yr != 0) << 1);
	 if (xr) {
	    if (x < blursize-1) xo = 1;
	    else xo = 0;
	    xr1 = 1 - xr;
	 }
	 if (yr) {
	    if (y < blursize-1) yo = 1;
	    else yo = 0;
	    yr1 = 1 - yr;
	 }
	 switch (interp) {	/* get blur weighting value */
				/* interpolate between adjacent values */

	   case 0: blurval = getblur(x,y); 
		   break;

	   case 1: blurval = xr1 * getblur(x,y) +
			      xr * getblur(x+xo,y); 
		   break;

	   case 2: blurval = yr1 * getblur(x,y) +
			      yr * getblur(x,y+yo); 
		   break;

	   case 3: blurval = yr1 * (xr1 * getblur(x,y) +
			             xr * getblur(x+xo,y)) +
			     yr  * (xr1 * getblur(x,y+yo) +
				     xr * getblur(x+xo,y+yo)); 
		   break;
	 }
         arrtot += blurval * readc(array, xmid+x, ymid+y);   
      }
    }
#ifdef DEBUG
  if (debug & NCSTIM && debugz & 1)
	ncfprintf (stderr,"recep %d %d %d: inten %g\n",
		npt->recnm1,npt->recnm2,npt->recnm3,arrtot);
#endif
    if (array==BACKARR) npt->backgnd = arrtot;
    else		npt->stim    = arrtot;
    npt->wavel = swavel[array];
    if (stimhlp) {
        ncfprintf (stderr,"s");
        if (++j >= 50) {
            ncfprintf (stderr,"\n");
            j=0;
        }
    }
  }
  if (stimhlp) ncfprintf (stderr,"\n");
#ifdef DEBUG
  if (debug & NCSTIM && debugz & 1)
		 ncfprintf (stderr,"recpat: end.\n");
#endif
}

/*------------------------------------*/

extern int arraymade;                           /* in modcode.cc */

void recback(int array, double inten, double wavel)
             
/* set receptor nodes and array to a diffuse
   background intensity. */
 
{
  recnod *npt;
  int j;

  for (npt=reclist,j=0; npt; npt=npt->next) {
    if (array==BACKARR) npt->backgnd = inten;
    else		npt->stim    = inten;
    npt->wavel = wavel;
    if (stimhlp) {
        ncfprintf (stderr,"b");
        if (++j >= 50) {
            ncfprintf (stderr,"\n");
            j=0;
        }
    }
  }
  if (stimhlp) ncfprintf (stderr,"\n");
  if (arraymade) setc (array,inten);
}

/*------------------------------------*/

static double saved_time = -1e20;
static int stimseq = 0;

void stimlist(double delratio, double time, double mask, int stimchan)
                   
/* Make a delta action list for stimulus from the list
   of receptor stimulus nodes.
   An action is an individual event at a receptor
   which consists of:

        1) time of event
        2) receptor number
        3) receptor intensity

    This list is used by "nc" to define the stimuli
    for the photoreceptors during a modeling run.
    It is output to the file "stimout", which may be "stdout".

    If simulation "time" is reset backwards, set new 
    sequence number to allow file to be sorted. 
*/

{
  recnod *npt;
  double newval;

  if (delratio < -1.0 || delratio > 1.0) {
        ncfprintf (stderr,"stimulus ratio out of range: %g\n",delratio);
        delratio = 0.0;
  }

  if ((delratio==1) && (saved_time > (simtime+BACKTIME))) {
      stimseq++;
  } 
  saved_time = simtime;

  for (npt=reclist; npt; npt=npt->next) {
    newval = delratio * npt->stim;
    delout (npt,newval,mask,stimchan,time,stimseq);
  }
}

/*------------------------------------*/

void abslist(double ratio, double time, double mask, int stimchan)

/* Make an absolute action list for stimulus from the list
   of receptor stimulus nodes.
   An action is an individual event at a receptor
   which consists of:

        1) time of event
        2) receptor number
        3) receptor intensity

    This list is used by "nc" to define the stimuli
    for the photoreceptors during a modeling run.
    It is output to the file "stimout", which may be "stdout".
*/

{
  recnod *npt;
  double newval;

  if (ratio < 0.0 || ratio > 1.0) {
        ncfprintf (stderr,"stimulus ratio out of range: %g\n",ratio);
        ratio = 1.0;
  }
  if ((ratio>=0) && (saved_time > (time+BACKTIME))) {
      stimseq++;
  } 
  saved_time = time;
  for (npt=reclist; npt; npt=npt->next) {
    newval = ratio * npt->stim + (1-ratio) * npt->backgnd;
    absout (npt,newval,mask,stimchan,time,stimseq);
  }
}

/*------------------------------------*/

void stcomment (void)
{
 ncfprintf(stimout,"#time node  intensity   wavel mask action sequence\n");
}

/*------------------------------------*/

void stimfout(double time, recnod *npt, double inten, 
		double mask, int stimchan, const char *action, int seq)
                      
/* Make one line of the stimulus file. */
/* This makes one action at one time for one node. */

/* Make NULLVAL (-32768) into a more elegant number */
/*  to show unused node dimensions in file. */

{
  int node1,node2,node3,node4;

  node1 = npt->recnm1;
  node2 = npt->recnm2;
  node3 = npt->recnm3;
  node4 = npt->recnm4;
  if (node4 == NULLVAL) node4 = NULND;
  if (node3 == NULLVAL) node3 = NULND;
  if (node2 == NULLVAL) node2 = NULND;
  if (node1 == NULLVAL) node1 = NULND;

/*   ncfprintf(stimout,"%-6.4g %-5d %-5d %-5d %-5d %-12.5g %-5g  %s\n", */

     ncfprintf(stimout,"%-.4g %-d %-d %-d %-d %-12.6g %-g %-g %-d %s %d\n",
			time, node1,node2,node3,node4,inten,npt->wavel,
			mask, stimchan, action, seq);
     fflush (stimout);
#undef BACKTIME 
}

/*------------------------------------*/

void stimclampout(double time, node *npt, double value, 
		int bgend, const char *action)
                      
/* Make one line of the stimulus file. */
/* This makes one action at one time for one node. */

/* Make NULLVAL (-32768) into a more elegant number */
/*  to show unused node dimensions in file. */

{
  int node1,node2,node3,node4;
  int seq = 0;

  node1 = npt->nodenm1;
  node2 = npt->nodenm2;
  node3 = npt->nodenm3;
  node4 = npt->nodenm4;
  if (node4 == NULLVAL) node4 = NULND;
  if (node3 == NULLVAL) node3 = NULND;
  if (node2 == NULLVAL) node2 = NULND;
  if (node1 == NULLVAL) node1 = NULND;

     ncfprintf(stimout,"%-.4g %-d %-d %-d %-d %-12.6g %-d %s %d\n",
			time, node1,node2,node3,node4,value,bgend,
			action, seq);
     fflush (stimout);

#undef BACKTIME 
}

/*------------------------------------*/

void absout (recnod *npt, double inten, double mask, int stimchan, double time, int seq)

/* send out a stimulus action for a receptor node
    if the change is greater than 1e-3.
 */

{
   double oldval,delta;

 npt->curval = inten;
 oldval = npt->actual;
 if (oldval == 0.0) delta = npt->curval;
 else delta = (npt->curval - oldval) / oldval;  /* look for change */
 if (abs(delta) > 1e-3) {                       /* save and output new val */
     npt->actual = npt->curval;
     if (mask>0) stimfout(time,npt,inten,mask,stimchan,"c",seq);
     else        stimfout(time,npt,inten,mask,stimchan,"a",seq);
 }
}

/*------------------------------------*/

void delout (recnod *npt, double delinten, double mask, int stimchan, double time, int seq)
               
/* Send out a delta stimulus action for a receptor node */

{
 if (delinten == 0.0) return;
 npt->curval += delinten;
 npt->actual = npt->curval;
 if (mask>0) stimfout(time,npt,delinten,mask,stimchan,"b",seq);
 else        stimfout(time,npt,delinten,mask,stimchan,"d",seq);
}

/*------------------------------------*/

void vclist (double time, int node1, int node2, 
		int node3, int node4, double inten, 
		double wavel, char *action)
{
  recnod x;

 x.recnm1 = node1;
 x.recnm2 = node2;
 x.recnm3 = node3;
 x.recnm4 = node4;
 x.wavel = wavel;
 stimfout(time,&x,inten,0.0,0,action,0);
}

/*------------------------------------*/

double gauss2(int x, int y, int midpnt, double rad)
                                        /* size is odd so midpnt is middle */
/* Return a gaussian function of (x,y),
   centered in an array of given size. */

{
   double distsq,tempx,tempy;

   tempx = (double)(x - midpnt) * SCALE;
   tempy = (double)(y - midpnt) * SCALE;
   distsq = tempx*tempx + tempy*tempy;
   return (exp(-(distsq) / (rad*rad)));
}

/*------------------------------------*/

double pointsp(int x, int y, int midpnt, 
		double scatter_rad, double ppow, double scatter_ampl)
                                        /* size is odd so midpnt is middle */

/* Return the optical point-spread function of 
   the cat eye, according to Robson and Enroth-Cugell (1978).
  
   I (r) = 1 / (1 + (r/k)^(5/2))

 The original linespread function:

   I (r) = 1 / (1 + (r/k)^(3/2))

   The rationale here is to generate a scatter function that can
be used with a variety of Gaussian blur functions.  This works OK 
when the Gaussian blur is made a little narrower and is added
to a scatter function that is the original point spread function
but is weighted so that it does not dominate the central peak.
The scatter function radius is related to the blur radius but is
larger to allow it to match correctly.  

See curve fit in "linespr.c" -- can fit a variety of
line-spread functions:
Calibrated in um.

 From Robson & Enroth-Cugell (1978)
   For a 4 mm pupil, cat.

  scatter_ampl = 0.15;
  scatter_rad = 45;
  ppow = 2.5;
  blur_rad = 11; 

 From Campbell and Gubisch (1966) 
   For a 5.7 mm pupil, human.  See "cg58.txt". 

  scatter_ampl = 10.0;
  scatter_rad = 2.7;
  ppow = 1.95;
  blur_rad = 1.5; 

 For "best refracted" average human eye, from Guirao et al, (2001).
   For a 5.7 mm pupil.  See "refr.txt".

   scatter_ampl = 1000.0;
   scatter_rad = 1.3;
   ppow = 1.85;
   blur_rad = 1.5; 

 For "uncorrected" average human eye, from Guirao et al, (2001). 
      for a 5.7 mm pupil.  See "uncorr.txt". 

   scatter_ampl = 1000.0;
   scatter_rad = 1.5;
   ppow = 1.0;
   blur_rad = 1.5; 

*/

#define CAT_SFACT 0.15
#define CAT_SCATTER_RAD 2.5
#define CAT_PPOW 2.5
#define GAUSSCALE 0.82	/* because 1/(1+pow(r,2.5) is bigger than Gaussian */

{
   double tempx,tempy,radius;
   double r;

   tempx = (double)(x - midpnt) * SCALE;
   tempy = (double)(y - midpnt) * SCALE;
   radius = sqrt(tempx*tempx + tempy*tempy); 
   r = radius/scatter_rad;
   return (scatter_ampl / (1 + pow(r,ppow)));

#undef CAT_SFACT
#undef CAT_SCATTER_RAD
#undef CAT_PPOW
#undef GAUSSCALE
}

/*------------------------------------*/

double pointspg(int x, int y, int midpnt, double scatter_rad, 
			double ppow, double scatter_ampl)

       /* size is odd so midpnt is middle */

/* Return the optical point-spread function of the human eye in
   fovea, based on on fit of line-spread function of Campbell and
   Gubisch (1966) by Geisler (1984) Values are in terms of 
   minutes of visual angle.

   Original function is sum of 2 Gaussians:

   y = a1 / (2*b1) * exp( -0.5 * (radius/s1)^2 ) +
       a2 / (2*b2) * exp( -0.5 * (radius/s2)^2 );

Where:
   minutes_degree = 60;
   microns_degree = 300;
   a1 = 0.684;
   b1 = 0.443;
   s1 = 0.443 / minutes_degree * microns_degree ->  2.215 um radius

   a2 = 0.587;
   b2 = 2.035;
   s2 = 2.035 / minutes_degree * microns_degree ->  10.175 um radius

This is equivalent to:

   y = .783 * ( 1.0  * exp (-0.5*(r/blur_rad)^2) +
               .1868 * exp (-0.5*(r/scatter_rad)^2));
   where:
         blur_rad   = blur radius
        scatter_rad = 10.175 um

   See "linesp.c" for original curve fit.

Since the smaller Gaussian is normally "optical blur from
diffraction" we only define the larger Gaussian when "scatter=2"
is turned on.

In this case, PFACT makes this scatter function the correct
amplitude with respect to the "blur" Gaussian which we don't
specify here -- it's correct when the central Gaussian matches
the Campbell and Gubisch (Geisler fit) blur of 2.215 um radius.  
For this blur we get 10.175 um radius for the scatter.

Since a Gaussian function is (x,y) separable, the radius of the 
point-spread function is identical to the corresponding
line-spread function and only the relative amplitude (PFACT,
below) needs to be changed.

We don't have to get the overall amplitude correct since the
total blur function is normalized after it is generated.

Radius of scatter function is not directly related to blur
function so "rad" only contains the scale factor. 

 From Campbell and Gubisch (1966) 
   For a 2 mm pupil.  See "linesp.c" for Geisler's (1984) fit. 

  scatter_ampl = 0.0427;
  scatter_rad = 10.175; 
  ppow = 0;
  blur_rad = 2.339;	 ( includes 4.43 optical and 1.5 cone aperture ) 

*/

{

/*   This original amplitude factor for the scatter function was
     correct for the line spread function only. A smaller amplitude
     for a correct point spread function is required. */

/* #define PFACT 0.1868		   /* = (a2/b2) /  (a1/b1) */

#define SCATTER_RAD2   4.5936795	   /* 10.175/2.215 = s2/s1 */
#define PFACT (1/(SCATTER_RAD2*SCATTER_RAD2*1.1077382)) /* = 0.04278007 */ 
						      /* ratio of volumes */
   double tempx,tempy,r,radius;

   tempx = (double)(x - midpnt) * SCALE;
   tempy = (double)(y - midpnt) * SCALE;
   radius = sqrt(tempx*tempx + tempy*tempy); 
   r = radius / (scatter_rad);
   return (scatter_ampl * exp(-0.5*r*r));

#undef PFACT
#undef SCATTER_RAD2
}

/*------------------------------------*/

int makblur(double blur_rad, double blur_ampl, double scale, 
		double scatter_rad, double scatter_pow, double scatter_ampl)

/* make the blurring function array */
/*  Determine the size of the array automatically.
    Make it 5 times the standard deviation of the 
    Gaussian blur function.  By this point, the Gaussian 
    function is down to 1e-11 (hope that's safe). 

    The "blur_ampl" parameter scales both blur and 
    scatter functions. This is useful for making a
    center-surround blur function using transducers.
*/

{
   int i,j,irow;
   double g,h,gausstot,normval,nrad,nscatter_rad;
   static int bsize=0;
   extern double *blurarr;
   static int blursiz;
 
#define SAFEVAL 1e-3
 
  if (blurmade) { 
	efree (blurarr);	/* erase old blur array */
	blurmade = 0;
  }
  if (!blurmade) { 
     blursiz = int(blur_rad/scale * 5 + 0.5); /* make array = 5 std. dev */
     if (scatter_ampl > 0) {
	if (scatter_pow > 0)  blursiz = int(pow(1/SAFEVAL,1.0/scatter_pow) *
					scatter_rad/scale);
	else                  blursiz = int(scatter_rad/scale*5);
     }
     if (blursiz < 5) blursiz = 5; 
     if (CONVSIZE<=256) blursiz = BLURSIZE;	/* small machine, limit size */
     bsize = sizeof(double) * blursiz * blursiz; 
     if (!(blurarr = (double *)smalloc(bsize))) {
	ncfprintf (stderr,"makblur: can't allocate blur array, size %d\n",
				bsize);
        blursiz = 0;
        return blursiz;
    }
    blurmade = 1;
    blurzero = 0; 
  }

  if (blur_rad==0.0) {
        blurzero = 1; 
	blursiz = 0;
        return blursiz;
  }
 
  if (scale==0.0) scale = 1.0;
  nrad = blur_rad/scale;
  nscatter_rad = scatter_rad/scale;
  for (i=0; i<blursiz; i++) {                   /* make 1/8 of square */
    irow = i*blursiz;
    for (j=0; j<=i; j++) {
      g = gauss2(i,j,0,nrad);
      if (scatter_ampl > 0) {			/* if scatter turned on */
        if (scatter_pow == 0) 
		 h = pointspg(i,j,0,nscatter_rad,scatter_pow, scatter_ampl);
        else     h = pointsp(i,j,0,nscatter_rad,scatter_pow, scatter_ampl);
        g += h;
      }
      if (g < 1e-15) g = 0.0;
      *(blurarr+irow+j) = g;
    }
  }

  for (i=0; i<blursiz; i++) {                   /* copy into quadrant */
    irow = i*blursiz;
    for (j=0; j<i; j++) {
      *(blurarr+j*blursiz+i) = *(blurarr+irow+j);
    }
  }

  gausstot = 0.0;
  for (i=0; i<blursiz; i++) {                     /* total the quadrant */
    irow = i*blursiz;
    for (j=1; j<blursiz; j++)                   /*  leave out the j=0 line */
      gausstot += *(blurarr+irow+j);
  }

  gausstot *= 4.0;                              /* four quadrants = plane */
  gausstot += *(blurarr+0+0);                   /* add the center */

/* ncfprintf (stderr,"gausstot %g blursiz %d\n",gausstot,blursiz); /* */

  if (gausstot <= 0.0) gausstot = 1;            /* normalize to vol = 1 */
  normval = blur_ampl / gausstot;		/* multiply by blur amplitude factor */
  for (i=0; i<blursiz; i++) {
    irow = i*blursiz;
    for (j=0; j<blursiz; j++)
      *(blurarr+irow+j) *= normval;
  }

/*
  i = 0;
  irow = i * blursiz;
  for (j=0; j<blursiz; j++) {
     ncfprintf (stderr,"%d %g\n",j,*(blurarr+irow+j));
  }
     ncfprintf (stderr,"%d\n",blursiz);
*/

  return (blursiz);
}

/*------------------------------------*/

double getblur(int x, int y)
           
/* The blur array is stored in a one quadrant
   look-up table. 
   Map x and y into this quadrant, then return
   with the blur value. */

{
 if (x < 0) x = -x;
 if (x >= blursize) x = blursize-1;
 if (y < 0) y = -y;
 if (y >= blursize) y = blursize-1;
 return ((double) *(blurarr+x*blursize+y));
}

/*------------------------------------*/

void makrect(int arr, int arrsize, double width, double length, 
		double xoff, double yoff, double scale, double orient, 
		double inten, double wavel)

/* Make a rectangle in the stimulus array. 
   If it is a narrow bar, there may be some sampling
   error in the array.  For usual amounts of blur, this 
   is corrected by finding the ratio of the wanted to 
   actual bar width and weighting the intensity by this 
   factor.  Must make sure to find actual smallest width 
   (not necessarily the stated "width", could be "length").
*/

{
   int x,y;
   double halfx,halfy,xcenter,ycenter;
   double totrect,weight,smalldim;
   int xmax,xmin,ymax,ymin,iscale;
   double minx,maxx,miny,maxy;
   double coso,sino,ycoso,ysino;
   double orientrad,xr,yr;
 
  if (scale==0.0) scale = 1.0;
  iscale = int(1.0/scale+0.5);
  swavel[arr] = wavel;
  xcenter =  xoff * iscale + arrsize/2;
  ycenter =  yoff * iscale + arrsize/2;
  halfx = width / 2 * iscale;
  halfy = length / 2 * iscale;

  xmax = (int)(xcenter + halfx+0.5); /* actual dimensions limited to array */
  xmax = min(xmax,arrsize); 
  xmin = (int)(xcenter - halfx+0.5);
  xmin = max(xmin,0); 
  ymax = (int)(ycenter + halfy+0.5);
  ymax = min(ymax,arrsize); 
  ymin = (int)(ycenter - halfy+0.5);
  ymin = max(ymin,0); 
  smalldim = min(width,length) * iscale;
 
			/* before making rect, add up its area */
  minx = xcenter - halfx;
  maxx = xcenter + halfx;
  miny = ycenter - halfy;
  maxy = ycenter + halfy;

  orientrad = MPI / 180 * orient;       /* rotate test point opposite way */
  coso = cos(orientrad);
  sino = sin(orientrad);

  totrect = 0.0;
  for (y=ymin, ysino=(ycenter-y)*sino + xcenter,
	       ycoso=(y-ycenter)*coso + ycenter; y<ymax; y++, ysino-=sino,ycoso+=coso) {
    for (x=xmin,xr=(x-xcenter)*coso+ysino,
		yr=(x-xcenter)*sino+ycoso; x<xmax; x++,xr+=coso,yr+=sino) {
      // if (inrect((double)x, (double)y, xcenter+halfx, xcenter-halfx,
      //                  ycenter+halfy, ycenter-halfy, orient))
       if ((minx <= xr) && (xr < maxx) && (miny <= yr) && (yr < maxy)) totrect++;
    }
  }

/*
  if (totrect==0.0) {
     if (width==0.0 || length==0.0) {
          ncfprintf (stderr,"Stim: Error: rect has zero width\n");
     }
     else {
	static int warned=0;
       if (!warned) {
	 warned = 1;
 	 ncfprintf (stderr,
		"Stim: warning: rect is outside stimulus array\n");
       }
      return;
     }
  }
*/
			/* then weight rect to correct its intensity */

/*
  if (smalldim > 20) weight = 1.0;
  else weight = (xmax-xmin) * (ymax-ymin+1) / totrect;
*/

  weight = 1.0;

  inten *= weight;
  for (y=ymin, ysino=(ycenter-y)*sino + xcenter,
	       ycoso=(y-ycenter)*coso + ycenter; y<ymax; y++, ysino-=sino,ycoso+=coso) {
    for (x=xmin,xr=(x-xcenter)*coso+ysino,
		yr=(x-xcenter)*sino+ycoso; x<xmax; x++,xr+=coso,yr+=sino) {
      // if (inrect((double)x, (double)y, xcenter+halfx, xcenter-halfx,
      //                  ycenter+halfy, ycenter-halfy, orient))
       if ((minx <= xr) && (xr < maxx) && (miny <= yr) && (yr < maxy)) addc (arr, x, y, inten);
    }
  }
  // stimlist() is done in modcode.cc, in xstim();
}

/*------------------------------------*/

void makbar(int arr, int arrsize, double width, double length, 
		double xloc, double yloc, double xcent, double ycent,
		double scale, double orient, double inten, 
		double start, double dur, double wavel, double mask, int stimchan)

/* Make a rectangle in the stimulus array.  */
/*   xcent,ycent represents the offset of the stimulus in the array space. */

{
   int x,y;
   double halfx,halfy,xcenter,ycenter;
   double xir,totrect,weight,smalldim;
   int xmax,xmin,ymax,ymin,iscale;
   double minx,maxx,miny,maxy;
   double coso,sino,ycoso,ysino;
   double orientrad,xr,yr;
 
  if (scale==0.0) scale = 1.0;
  iscale = int(1.0/scale+0.5);
  swavel[arr] = wavel;
  xcenter =  (xloc-xcent) * iscale + arrsize/2;
  ycenter =  (yloc-ycent) * iscale + arrsize/2;
  halfx = width / 2 * iscale;
  halfy = length / 2 * iscale;

  xmax = (int)(xcenter + halfx+0.5); /* actual dimensions limited to array */
  xmax = min(xmax,arrsize); 
  xmin = (int)(xcenter - halfx+0.5);
  xmin = max(xmin,0); 
  ymax = (int)(ycenter + halfy+0.5);
  ymax = min(ymax,arrsize); 
  ymin = (int)(ycenter - halfy+0.5);
  ymin = max(ymin,0); 
  smalldim = min(width,length) * iscale;
 
			/* before making rect, add up its area */
  minx = xcenter - halfx;
  maxx = xcenter + halfx;
  miny = ycenter - halfy;
  maxy = ycenter + halfy;

  orientrad = MPI / 180 * orient;       /* rotate test point opposite way */
  coso = cos(orientrad);
  sino = sin(orientrad);

  totrect = 0.0;
  for (y=ymin, ysino=(ycenter-y)*sino + xcenter,
	       ycoso=(y-ycenter)*coso + ycenter; y<ymax; y++, ysino-=sino,ycoso+=coso) {
    for (x=xmin,xr=(x-xcenter)*coso+ysino,
		yr=(x-xcenter)*sino+ycoso; x<xmax; x++,xr+=coso,yr+=sino) {
      // if (inrect((double)x, (double)y, xcenter+halfx, xcenter-halfx,
      //                  ycenter+halfy, ycenter-halfy, orient))
       if ((minx <= xr) && (xr < maxx) && (miny <= yr) && (yr < maxy)) totrect++;
    }
  }
/*
  if (totrect==0.0) {
     if (width==0.0 || length==0.0) {
          ncfprintf (stderr,"Stim: Error: rect has zero width\n");
     }
     else {
	static int warned=0;
       if (!warned) {
	 warned = 1;
 	 ncfprintf (stderr,
		"Stim: warning: rect is outside stimulus array\n");
       }
      return;
     }
  }
*/
			/* then weight rect to correct its intensity */

/*
  if (smalldim > 20) weight = 1.0;
  else weight = (xmax-xmin) * (ymax-ymin+1) / totrect;
*/

  weight = 1.0;

  inten *= weight;
  for (y=ymin, ysino=(ycenter-y)*sino + xcenter,
	       ycoso=(y-ycenter)*coso + ycenter; y<ymax; y++, ysino-=sino,ycoso+=coso) {
    for (x=xmin,xr=(x-xcenter)*coso+ysino,
		yr=(x-xcenter)*sino+ycoso; x<xmax; x++,xr+=coso,yr+=sino) {
      // if (inrect((double)x, (double)y, xcenter+halfx, xcenter-halfx,
      //                  ycenter+halfy, ycenter-halfy, orient))
       if ((minx <= xr) && (xr < maxx) && (miny <= yr) && (yr < maxy)) addc (arr, x, y, inten);
    }
  }

  if (int(stimdia/scale+0.5) & 1) xir = 0; else xir = 0.5;

  recpat(arr, arrsize, xcent+xir*scale,ycent,scale); /* set photrec inten */
  stimlist(1.0, start, mask,stimchan);               /* make an action list */
  stimlist(-1.0, start+dur, mask,stimchan);
  recback(arr,0.0,1.0);   /* zero photrec stim inten */
}

/*------------------------------------*/

void makcheckerboard(int arr, int arrsize, double xsize, double ysize, int xn, int yn,
		double xloc, double yloc, double xcent, double ycent, 
		double scale, double orient, double tfreq, 
		double inten, double contrast, double start,double dur, 
		double **stim_rndarr, int *stim_nfr, int ckrseed)

/* xloc,yloc represent the position of the checkerboard in neural space.
   xcent,ycent represent the offset of the stimulus in the array space. */

{
  int i,j,k,indx,nframes,framesiz,kindx;
  int iscale,stimchan;
  int xmax,xmin,ymax,ymin;
  double x,y,t;
  double pixwidth, pixheight;
  double xcenter, ycenter;
  double hpixsize, vpixsize;
  double tperiod, tincr, stimend;
  double orientrad, coso, sino, ycoso, ysino;
  double minx,maxx,miny,maxy;
  double rand_inten, stime, sinten, wavel, mask;
  double xoff, yoff, halfx, halfy, xr, yr;
  double *rand_arr;
  char *chpnt;
  unsigned int chseed;

  if (tfreq==0.0) tfreq = 1.0;
  tincr = 1.0 / tfreq;               /* temporal period */
  tincr = abs(tincr);                /* allow negative temporal frequency */
  if (tincr < stiminc) tincr = stiminc; /* time incr = photrec timinc */
  nframes = dur / tincr + 0.5;
  stimend = start + tincr*nframes + RNDUP;

  if (xn < 1) xn = 1;		/* make at least 1 check in checkerboard */
  if (yn < 1) xn = 1;
  pixwidth  = xsize / xn;	/* size of each check */
  pixheight = ysize / yn;
  framesiz = xn * yn;		/* total number of checks in checkerboard */

  /* make a new array to contain the spatio-temporal stimulus */
		    
  if ((rand_arr = (double*)emalloc(xn * yn * nframes * sizeof(double)))==NULL) {
      ncfprintf (stderr,"nc: checkerboard, can't init array\n");
      return;
  }
  if (stim_rndarr != NULL && stim_nfr != NULL) {
     *stim_rndarr = rand_arr; 	/* remember the stimulus array */
     *stim_nfr = nframes;
  }

  /* make a new random number generator and initialize it */

  chseed = ckrseed^7924517;
  chpnt = makrand(RNDSIZ,"checkerboard",1);
  initstate (chseed,chpnt,RNDSIZ);
  drand_setstate(chpnt);

  /* load the array with the random stimulus */

  for (k=indx=0; k<nframes; k++) {
    for (j=0; j<yn; j++) {
      for (i=0; i<xn; i++) {
        rand_arr[indx++] = contrast * drand();
      }
    }
  }
  restorstate();                /* restore orig random sequence */

  orientrad = MPI / 180.0 * orient;
  coso = cos(orientrad);
  sino = sin(orientrad);
  wavel = 0;
  mask = 0;
  stimchan = 0;
  swavel[arr] = wavel;

  if (scale==0.0) scale = 1.0;
  halfx = pixwidth  * 0.5 / scale;
  halfy = pixheight * 0.5 / scale;

  for (k=indx=0; k<nframes; k++) {	/* for all frames */
    for (j=0; j<yn; j++) {  		/* for all pixels */
      for (i=0; i<xn; i++) {

        sinten = inten + rand_arr[indx++];

        x = i*pixwidth  + pixwidth*0.5;	/* relative loc of unrotated rectangle */
        y = j*pixheight + pixheight*0.5;
        xoff = coso*x - sino*y;		/* rel loc of rotated rectangle */
        yoff = sino*x + coso*y;

        xcenter =  (xoff+xloc) / scale + arrsize/2; /* location in stim array */
        ycenter =  (yoff+yloc) / scale + arrsize/2;

        xmax = (int)(xcenter + halfx+0.5); /* actual dimensions limited to array */
        xmax = min(xmax,arrsize); 
        xmin = (int)(xcenter - halfx+0.5);
        xmin = max(xmin,0); 
        ymax = (int)(ycenter + halfy+0.5);
        ymax = min(ymax,arrsize); 
        ymin = (int)(ycenter - halfy+0.5);
        ymin = max(ymin,0); 
 
        minx = xcenter - halfx;
        maxx = xcenter + halfx;
        miny = ycenter - halfy;
        maxy = ycenter + halfy;


        for (y=ymin, ysino=(ycenter-y)*sino + xcenter,
	         ycoso=(y-ycenter)*coso + ycenter; y<ymax; y++, ysino-=sino,ycoso+=coso) {
          for (x=xmin,xr=(x-xcenter)*coso+ysino,
		  yr=(x-xcenter)*sino+ycoso; x<xmax; x++,xr+=coso,yr+=sino) {
             if ((minx <= xr) && (xr < maxx) && (miny <= yr) && (yr < maxy)) 
		     	addc (arr, x, y, sinten);
          }
        }
      } /* for (i;;) */
    }  /* for (j;;) */

    stime = start + k * tincr;

    recpat(arr,arrsize,xcent,ycent,scale);  /* set recept inten */
    stimlist(1.0, stime,mask,stimchan);		/* make a delta action list */
    stimlist(-1.0, stime+tincr,mask,stimchan);
    recback(arr,0.0,1.0);     		/* zero recept stim inten */
  }   /* for (stime= ; ; ) */

}     /* makcheckerboard() */


/*---------------------------------------------------------------------*/

double gaussf1s(double x, double xcent, double rad)

        /* Return a gaussian function of x,
         *    centered in an array of given size. */

{
         double distsq,tempx;

   tempx = (double)(x - xcent);
   distsq = tempx*tempx;
   return (exp(-(distsq) / (rad*rad)));
}

/*---------------------------------------------------------------------*/

double gaussf1s(double x, double xcent, double rad, double rlim)

        /* Return a gaussian function of x,
         *    centered in an array of given size. */

{
         double distsq,tempx;

   tempx = (double)(x - xcent);
   distsq = tempx*tempx;
   if (x > rlim) return 0;
   else
   return (exp(-(distsq) / (rad*rad)));
}

/*---------------------------------------------------------------------*/

double sinsq2 (double theta, int sq)

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

/*------------------------------------*/

#define STIMINC dsintinc	/* min timinc for sine wave (0.002 sec) */
#define STIMRES dsintres	/* time res for sine waves (0.002) */

void makgabor(int arr, int arrsize, double contrast, double inten_add, double inten_mult,
		double xloc, double yloc,
		double xcent, double ycent,
		double xenv, double yenv, 
		int makenv, int sq, double tfreq, double drift, 
		double speriod, double sphase, 
		double scale, double orient, double wavel,
		double mask, int stimchan, double start, double dur)

/* Make a Gabor filter (sine wave inside Gaussian envelope) 
   in the stimulus array.
*/

{
    int x,y,tmod;
    double xcenter,ycenter;
    double stoff,sinttheta,sinten;
    double orientrad,sphasenorm,coso,sino,ysino,ycoso;
    double xedgemax, xedgemin, yedgemax, yedgemin;
    double xlim, ylim, gcontrast;
    double xr,yr,xradsq,yradsq, stheta,sinstheta;
    double  toffset, tperiod, tincr, stime, stimend;
    int xmax,xmin,ymax,ymin;
    int amid;

  if (scale==0.0) scale = 1.0;
  swavel[arr] = wavel;

  amid = arrsize/2;
  xcenter =  (xloc-xcent) / scale + amid;
  ycenter =  (yloc-ycent) / scale + amid;
  xenv /= scale;
  yenv /= scale;

  if (makenv) {                         /* use envelope function */
    xlim = xenv * 5;
    ylim = yenv * 5;
  } 
  else {                                /* sine, make sharp border */
    xlim = xenv; 		/* limit for computing this stimulus */
    ylim = yenv; 
  }
  // xmax = int(xcenter + xlim + 0.5); /* actual dims limited to array */
  // xmin = int(xcenter - xlim + 0.5);
  // ymax = int(ycenter + ylim + 0.5);
  // ymin = int(ycenter - ylim + 0.5);

  // xmax = min(xmax,arrsize); 
  // xmin = max(xmin,0); 
  // ymax = min(ymax,arrsize); 
  // ymin = max(ymin,0); 

  xmax = arrsize; 
  xmin = 0;
  ymax = arrsize;
  ymin = 0; 

  if (speriod <= 0.0) speriod = 1;      /* spatial period from user */
  speriod /= scale;

  orientrad =   MPI / 180.0 * orient;
  sphasenorm =  sphase/360.0;
  coso = cos(orientrad);
  sino = sin(orientrad);
  if (drift || tfreq) {		/* drifting or counterphase grating */
     tmod = 1;			/* temporal modulation */
     if (tfreq==0.0) tfreq = 1.0;
     tperiod = 1.0 / tfreq;	/* temporal period */
     tincr = tperiod * STIMRES;
     tincr = abs(tincr);	/* allow negative temporal frequency */
     if (tincr < STIMINC) tincr = STIMINC; /* time incr = rec.timinc*/
     stimend = start+dur+RNDUP;
  }
  else {			/* static grating, no temporal modulation */
      tmod = 0;			/* no temporal modulation */
      stoff = 0.0;		/* no spatio-temporal offset */
      sinttheta = 1.0;		/* no counterphase modulation */
      stimend = start+dur-RNDUP; /* stop after one step */
      tincr = dur;		/* dummy tincr */
  }

  for (stime=start; stime<stimend; stime += tincr) {

    if (stime+tincr >= stimend) {
      tincr = stimend-stime;            /* fix last step so it doesn't overflow */
    }
    if (tmod) {
      toffset = (start-stime)/tperiod;
      if (drift) {		     /* drifting grating */
        stoff = drift * toffset;     /* add temporal shift into spatial phase */
        sinttheta = 1.0;
      }
      else	{		     /* stationary counterphasing grating */
        stoff = 0.0;
        sinttheta = sin(2*MPI * toffset); /* temporal shift, indep sine func */
      }
    }
    xradsq=xenv*xenv;
    yradsq=yenv*yenv;
    gcontrast = inten_mult * contrast * sinttheta;
    for (y=ymin,ysino=(ycenter-y)*sino,
		ycoso=(y-ycenter)*coso; y<ymax; y++,ysino-=sino,ycoso+=coso) {
       for (x=xmin,xr=(x-xcenter)*coso+ysino,
		   yr=(x-xcenter)*sino+ycoso; x<xmax; x++,xr+=coso,yr+=sino) {
        //
	// save time by not calling gabormod() with all its parameters
        // sinten = inten_add + gabormod(contrast, inten_mult, x, y, xcenter, ycenter, 
	//		  xenv, yenv, sino, coso, sphasenorm, speriod, 
	//		  sinttheta, stoff);
	 stheta = (xr-xcent)/speriod + stoff + sphasenorm;
	 sinstheta = sinsq2(2*MPI * stheta,sq) * exp(-xr*xr/xradsq) * exp(-yr*yr/yradsq);
	 if (xr > xlim && yr > ylim) sinstheta = 0;
         sinten = inten_add + gcontrast * sinstheta; 
         addc (arr, x, y, sinten);
       }
    }
    recpat(arr,arrsize,xcent,ycent,scale);  /* set recept inten */
    stimlist(1.0, stime,mask,stimchan);     /* make a delta action list */
    stimlist(-1.0, stime+tincr,mask,stimchan);
    recback(arr,0.0,1.0);     /* zero recept stim inten */

  }   /* for (stime= ; ; ) */

#ifdef DEBUG
  if (debug & NCSTIM && debugz & 32)
    printa(arr); 		/* print slice through center for debugging */
#endif

}  

/*------------------------------------*/

void makgabor(int arr, int arrsize, double contrast, double inten_mult,
		double xloc, double yloc,
		double xcent, double ycent,
		double xenv, double yenv, 
		int makenv, int sq, double tfreq, double drift, 
		double speriod, double sphase, 
		double scale, double orient, double wavel,
		double mask, double start, double dur)
{
     double inten_add;
     int stimchan;

   makgabor(arr, arrsize, contrast, inten_add=0, inten_mult, xloc, yloc, xcent, ycent, 
		     xenv, yenv, makenv, sq, tfreq, drift, speriod, sphase, 
		     scale, orient, wavel, mask, start, dur, stimchan=0);
}

/*------------------------------------*/

void maksineann(int arr, int arrsize, double contrast, double inten_add, double inten_mult,
		double xloc, double yloc, double xcent, double ycent, 
		double renv, int makenv, int sq, double tfreq,
		double drift, double speriod, double sphase, 
		double scale, double wavel,
		double mask, int stimchan, double start, double dur)

/* Make a sine wave annulus (radial sine wave inside Gaussian envelope) 
   in the stimulus array.
*/

{
    int x,y,tmod;
    double xcenter,ycenter;
    double stoff,sinttheta,stheta,sinten;
    double gcontrast,sphasenorm;
    double rlim,r;
    double xtemp,ytemp,ytemp2;
    double toffset, tperiod, tincr, stime, stimend;
    int xmax,xmin,ymax,ymin;
    int amid;

  if (scale==0.0) scale = 1.0;
  swavel[arr] = wavel;

  amid = arrsize/2;  
  xcenter =  (xloc-xcent) / scale + amid;
  ycenter =  (yloc-ycent) / scale + amid;
  renv /= scale;

 if (makenv) {				/* use envelope function */
    rlim = renv * 5;
 } 
 else {					/* sine, make sharp border */
    rlim = renv;
  }
  // xmax = int(xcenter + rlim + 0.5); /* actual dims limited to array */
  // xmin = int(xcenter - rlim + 0.5);
  // ymax = int(ycenter + rlim + 0.5);
  // ymin = int(ycenter - rlim + 0.5);

  // xmax = min(xmax,arrsize); 
  // xmin = max(xmin,0); 
  // ymax = min(ymax,arrsize); 
  // ymin = max(ymin,0); 

  xmax = arrsize; 
  xmin = 0;
  ymax = arrsize;
  ymin = 0; 

  if (speriod <= 0.0) speriod = 1;      /* spatial period from user */
  speriod /= scale;

  sphasenorm =  sphase/360.0;
  if (drift || tfreq) {		/* drifting or counterphase grating */
     tmod = 1;			/* temporal modulation */
     if (tfreq==0.0) tfreq = 1.0;
     tperiod = 1.0 / tfreq;	/* temporal period */
     tincr = tperiod * STIMRES;
     tincr = abs(tincr);	/* allow negative temporal frequency */
     if (tincr < STIMINC) tincr = STIMINC; /* time incr = rec.timinc*/
     stimend = start+dur+RNDUP;
  }
  else {			/* static grating, no temporal modulation */
      tmod = 0;			/* no temporal offset */
      stoff = 0.0;		/* no spatio-temporal offset */
      sinttheta = 1.0;		/* no temporal modulation */
      stimend = start+dur-RNDUP;/* stop after one step */
      tincr = dur;		/* dummy tincr */
  }

  for (stime=start; stime<stimend; stime += tincr) {

    if (stime+tincr >= stimend) {
      tincr = stimend-stime;            /* fix last step so it doesn't overflow */
    }
    if (tmod) {
      toffset = (start-stime)/tperiod;
      if (drift) {			/* drifting grating */
        stoff = drift * toffset;	/* add temporal shift to spatial phase */
        sinttheta = 1.0;
      }
      else	{			/* stationary counterphasing grating */
        stoff = 0.0;
        sinttheta = sin(2*MPI * toffset); /* temporal shift is indep sine func */
      }
    }

    // sphasenorm += stoff;
    gcontrast = inten_mult * contrast * sinttheta;
    for (y=ymin,ytemp=y-ycenter; y<ymax; y++,ytemp++) {
      ytemp2 = ytemp*ytemp;
      for (x=xmin,xtemp=x-xcenter; x<xmax; x++,xtemp++) {
        //
	// save time by not calling sineannmod() with all its parameters
        // sinten = inten_add + sineannmod(contrast, inten_mult, x, y, xcenter, ycenter, renv,
        //                             sphasenorm, speriod, sinttheta, stoff);
	r = sqrt(xtemp*xtemp + ytemp2);
	stheta = r/speriod + stoff + sphasenorm;
        sinten = inten_add + sinsq2(MPI*2*stheta,sq) * gcontrast * gaussf1s(r,0,renv,rlim);
        addc (arr, x, y, sinten);
      }
    }
    recpat(arr,arrsize,xcent,ycent,scale);  /* set recept inten */
    stimlist(1.0, stime,mask,stimchan);     /* make a delta action list */
    stimlist(-1.0, stime+tincr,mask,stimchan);
    recback(arr,0.0,1.0);     /* zero recept stim inten */

  }   /* for (stime= ; ; ) */

  if (debug & NCSTIM && debugz & 32)
    printa(arr); 		/* print slice through center for debugging */
}  

/*------------------------------------*/

void maksineann(int arr, int arrsize, double contrast, double inten_mult,
		double xloc, double yloc, double xcent, double ycent, 
		double renv, int makenv, int sq, double tfreq,
		double drift, double speriod, double sphase, 
		double scale, double wavel,
		double mask, double start, double dur)

{
     double inten_add;
     int stimchan;

  maksineann(arr, arrsize, contrast, inten_add=0, inten_mult, xloc, yloc, xcent, ycent, 
		  renv, makenv, sq, tfreq, drift, speriod, sphase, 
		  scale, wavel, mask, stimchan=0, start, dur);
}

/*------------------------------------*/

void makwindmill(int arr, int arrsize, double contrast, double inten_add, double inten_mult,
		double xloc, double yloc, double xcent, double ycent, 
		double renv, int makenv, int sq, double tfreq,
		double drift, double speriod, double sphase, 
		double scale, double start, double dur, double wavel, double mask, int stimchan)

/* Make a windmill (sine wave cross inside Gaussian envelope) 
   in the stimulus array.
*/

{
    int x,y,tmod;
    double xcenter,ycenter;
    double stoff,sinttheta,stheta,sinten;
    double gcontrast,sphasenorm;
    double xtemp,ytemp;
    double rlim,r,t;
    double toffset, tperiod, tincr, stime, stimend;
    int xmax,xmin,ymax,ymin;
    int amid;

  if (scale==0.0) scale = 1.0;
  swavel[arr] = wavel;

  amid = arrsize/2;  
  xcenter =  (xloc-xcent) / scale + amid;
  ycenter =  (yloc-ycent) / scale + amid;
  renv /= scale;

 if (makenv) {				/* use envelope function */
    rlim = renv * 5;
 } 
 else {					/* sine, make sharp border */
    rlim = renv; 
  }
  // xmax = int(xcenter + rlim + 0.5); /* actual dims limited to array */
  // xmin = int(xcenter - rlim + 0.5);
  // ymax = int(ycenter + rlim + 0.5);
  // ymin = int(ycenter - rlim + 0.5);

  // xmax = min(xmax,arrsize); 
  // xmin = max(xmin,0); 
  // ymax = min(ymax,arrsize); 
  // ymin = max(ymin,0); 

  xmax = arrsize; 
  xmin = 0;
  ymax = arrsize;
  ymin = 0; 

  if (speriod <= 0.0) speriod = 1;      /* spatial period from user */
  speriod /= scale;
  sphasenorm =  sphase/360.0;
  if (drift || tfreq) {		/* drifting or counterphase grating */
     tmod = 1;			/* temporal modulation */
     if (tfreq==0.0) tfreq = 1.0;
     tperiod = 1.0 / tfreq;	/* temporal period */
     tincr = tperiod * STIMRES;
     tincr = abs(tincr);	/* allow negative temporal frequency */
     if (tincr < STIMINC) tincr = STIMINC; /* time incr = rec.timinc*/
     stimend = start+dur+RNDUP;
  }
  else {			/* static grating, no temporal modulation */
      tmod = 0;			/* no temporal offset */
      stoff = 0.0;		/* no spatio-temporal offset */
      sinttheta = 1.0;		/* no temporal modulation */
      stimend = start+RNDUP;	/* stop after one step */
      tincr = dur;		/* dummy tincr */
  }

  for (stime=start; stime<stimend; stime += tincr) {

    if (stime+tincr >= stimend) {
      tincr = stimend-stime;            /* fix last step so it doesn't overflow */
    }
    if (tmod) {
      toffset = (start-stime)/tperiod;
      if (drift) {			/* drifting grating */
        stoff = -drift * toffset;	/* add temporal shift to spatial phase */
        sinttheta = 1.0;
      }
      else	{			/* stationary counterphasing grating */
        stoff = 0.0;
        sinttheta = sin(2*MPI * toffset); /* temporal shift is indep sine func */
      }
    }

    sphasenorm += stoff;
    gcontrast = inten_mult*contrast*sinttheta;
    for (y=ymin,ytemp=y-ycenter; y<ymax; y++,ytemp++) {
      for (x=xmin,xtemp=x-xcenter; x<xmax; x++,xtemp++) {
        //sinten = inten_add + windmillmod(contrast, inten_mult, x, y, xcenter, ycenter, renv,
        //                             sphasenorm, speriod, sinttheta, stoff);
	t = axy(xtemp,ytemp);
	r = sqrt(xtemp*xtemp + ytemp*ytemp);
	stheta = t*speriod/(MPI*2) + sphasenorm;
        sinten = inten_add + sinsq2(MPI*2*stheta,sq) * gcontrast * gaussf1s(r,0,renv,rlim);
        addc (arr, x, y, sinten);
      }
    }
    recpat(arr,arrsize,xcent,ycent,scale);  /* set recept inten */
    stimlist(1.0, stime,mask,stimchan);     /* make a delta action list */
    stimlist(-1.0, stime+tincr,mask,stimchan);
    recback(arr,0.0,1.0);     /* zero recept stim inten */

  }   /* for (stime= ; ; ) */

  if (debug & NCSTIM && debugz & 32)
    printa(arr); 		/* print slice through center for debugging */
}  

/*------------------------------------*/

void makwindmill(int arr, int arrsize, double contrast, double inten_mult,
		double xloc, double yloc, double xcent, double ycent, 
		double renv, int makenv, int sq, double tfreq,
		double drift, double speriod, double sphase, 
		double scale, double start, double dur, double wavel, double mask)
{
    double inten_add;
    int stimchan;

  makwindmill(arr, arrsize, contrast, inten_add=0, inten_mult, xloc, yloc, xcent, ycent, 
			renv, makenv, sq, tfreq, drift, speriod, sphase, 
			scale, start, dur, wavel, mask, stimchan=0);
}

/*------------------------------------*/

void makimage(int arr, int arrsize, double width, double length, 
		double xoff, double yoff, double orient, double scale, 
		double inten, int type, const char *imgfil)

/* Make an image in the stimulus array.  */

{
   int x,y,val;
   double halfx,halfy,xcenter,ycenter;
   double minx,maxx,miny,maxy;
   int xmax,xmin,ymax,ymin,iscale;
   FILE *fimg;


  if ((fimg=fopen(imgfil,"r"))==NULL) { /* open file */
	ncfprintf (stderr,"%s: can't open image file '%s'\n", "stim",imgfil);
	return;
  }

  if (scale==0.0) scale = 1.0;
  iscale = int(1.0/scale+0.5);
  xcenter =  xoff * iscale + arrsize/2;
  ycenter =  yoff * iscale + arrsize/2;
  halfx = width / 2 * iscale;
  halfy = length / 2 * iscale;

  xmax = (int)(xcenter + halfx+0.5); /* actual dimensions limited to array */
  xmax = min(xmax,arrsize); 
  xmin = (int)(xcenter - halfx+0.5);
  xmin = max(xmin,0); 
  ymax = (int)(ycenter + halfy+0.5);
  ymax = min(ymax,arrsize); 
  ymin = (int)(ycenter - halfy+0.5);
  ymin = max(ymin,0); 

  minx = xcenter - halfx;
  maxx = xcenter + halfx;
  miny = ycenter - halfy;
  maxy = ycenter + halfy;

  val = 0; 
  for (y=ymin; y<=ymax; y++) {
    for (x=xmin; x<xmax; x++) {
       fread(&val,1,1,fimg);
        // if (inrect((double)x, (double)y, xcenter+halfx, xcenter-halfx,
        //                ycenter+halfy, ycenter-halfy, orient))
        //  addc (arr, x, y, inten+val);
        if ((minx <= x) && (x < maxx) && (miny <= y) && (y < maxy)) addc (arr, x, y, inten+val);
    }
  }
}

/*------------------------------------*/

void makspot(int arr, int arrsize, double dia, 
		double xloc, double yloc, double xcent, double ycent, 
		double scale, double inten, double start, double dur, 
		double wavel, double mask, int stimchan, int invert) 

/* Make a spot in the stimulus array.
   If it is a small spot, there may be some sampling
   error in the array.  If blur is larger than the spot 
   size (usually true for a small spot), this error is 
   corrected by finding the ratio of the wanted to actual 
   spot area and weighting the intensity by this factor.
*/

{
   int x,y,iscale,xtemp,ytemp,ytemp2;
   double xcenter,ycenter,radius,totspot,weight;
   double xmax,xmin,ymax,ymin;

  if (scale==0.0) scale = 1.0;
  iscale = int(1.0/scale+0.5);
  swavel[arr] = wavel;
  radius = dia * iscale / 2.0;
  xcenter =  (xloc-xcent) * iscale + arrsize/2;
  ycenter =  (yloc-ycent) * iscale + arrsize/2;
  xmax = xcenter + radius + 1;
  xmax = min(xmax,arrsize); 
  xmin = xcenter - radius - 1;
  xmin = max(xmin,0); 
  ymax = ycenter + radius + 1;
  ymax = min(ymax,arrsize); 
  ymin = ycenter - radius - 1;
  ymin = max(ymin,0); 

			/* before making spot, add up its area */

  for (y=(int)ymin,ytemp=y-ycenter; y<(int)ymax; y++,ytemp++) {
    ytemp2 = ytemp*ytemp;
    for (x=(int)xmin,xtemp=x-xcenter; x<(int)xmax; x++,xtemp++) {
      //if (incirc((double)x, (double)y, xcenter, ycenter, radius))
      if ((sqrt(xtemp*xtemp+ytemp2)<=radius) ^ invert) totspot++;
    }
  }
  if (totspot==0.0) {
     if (dia==0.0) {
          ncfprintf (stderr,"Stim: Error: spot has zero diameter\n");
     }
     else {
	static int warned=0;
       if (!warned) {
	 warned = 1;
 	 ncfprintf (stderr,"Stim: warning: spot is outside stimulus array\n");
       }
       return;
     }
  }

  if (dia > 20) weight = 1.0;   /* weight spot to correct its intensity */
  else          weight = MPI * radius * radius / totspot;
  inten *= weight;

  for (y=(int)ymin,ytemp=y-ycenter; y<(int)ymax; y++,ytemp++) {
    ytemp2 = ytemp*ytemp;
    for (x=(int)xmin,xtemp=x-xcenter; x<(int)xmax; x++,xtemp++) {
      //if (incirc((double)x, (double)y, xcenter, ycenter, radius))
      if ((sqrt(xtemp*xtemp+ytemp2)<=radius) ^ invert ) addc (arr, x, y, inten);
    }
  }
  recpat(arr,arrsize,xcent,ycent,scale); /* set photrec inten */
  stimlist(1.0, dur,mask,stimchan);
  recback(arr,0.0,1.0);     /* zero photrec stim inten */
}

/*------------------------------------*/

recnod *maksnod(void)

/* make a new stimulus and link it to the list. */

{
    recnod *rpnt;

  if ((rpnt=(recnod *)emalloc(sizeof(recnod))) == NULL) {
     ncfprintf (stderr,"no space left for recnod %d\n", reccum+1);
     return ((recnod*)NULL);  
  }
  rpnt->next = (recnod*)NULL;
  rpnt->last = reclend;
  if (!reclist) reclist = rpnt;         /* save head if first synap */
  if (reclend)
    reclend->next = rpnt;
  reclend = rpnt;

  reccum++;                     /* increment total */
  return (rpnt); 
}

/*------------------------------------*/

void delsnod(recnod *rpnt)

/* delete a stimulus and remove it from the list. */

{
  if (!rpnt) return;
  if (rpnt->next) rpnt->next->last = rpnt->last;
  if (rpnt->last) rpnt->last->next = rpnt->next;
  if (reclist==rpnt) reclist = rpnt->last;         /* save head if first synap */
  if (reclend==rpnt) reclend = rpnt->next;
  efree(rpnt);
  reccum--;                     /* increment total */
  return; 
}

/*------------------------------------*/

void erasrecnod(elem *epnt)

/* erase a recnod if the element is a photoreceptor */

{
   recnod *rpnt, *nrpnt;
   int nod1,nod2,nod3,nod4;

  switch (epnt->ctype) {
    case ROD:
    case CONE:
    case VTRANSDUCER:
    case ITRANSDUCER:
     nod1 = epnt->node1a;
     nod2 = epnt->node1b;
     nod3 = epnt->node1c;
     nod4 = epnt->node1d;
     for (rpnt=reclist; rpnt; rpnt=nrpnt) {
       nrpnt = rpnt->next;
       if (rpnt->recnm1==nod1 && rpnt->recnm2==nod2 &&
           rpnt->recnm3==nod3 && rpnt->recnm4==nod4) {
         delsnod(rpnt);
       }
     }
     break;
     default: break;
  }
}

/*------------------------------------*/

