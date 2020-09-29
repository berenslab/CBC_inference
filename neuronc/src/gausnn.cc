/* Module gausnn in program gausnn */

/* Makes an array with gaussian distribution of
   nearest-neighbor distance */

/* This module can be used to build a standalone program,
    or it may be included in "nc" as a subroutine. */

 
/*  R.G.Smith */

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif

#ifdef __cplusplus
}
#endif

#include "ncio.h"

#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)>(b) ? (a) : (b))

#define MAXPTS 20000

#define FRAMEX  200  /* Frame size in microns; X-dimension */
#define FRAMEY  200  /* Frame size in microns; Y-dimension */

#define XCENT  0     /* location of X center in microns */

#define NUMCELLS 400
#define NUMN 10

char outfilnam[] = {"gausn.dat"};

FILE *outfil=(FILE*)NULL;
static FILE *fp;
static FILE *textout=stdout;

#ifdef __cplusplus
extern "C" {
#endif

 double sqrt(double);

#ifdef __cplusplus
}
#endif

#ifndef NONC
#include "ncsub.h"
#include "ncomp.h"
int noduniq(nodeint nodea, nodeint nodeb, 
			nodeint nodec, nodeint noded, int seed);
char *makrand (int rndsiz, const char *mesg, int num);
void initstate (unsigned long int seed, char *state, int n);
void drand_setstate (char *arg_state);
void restorstate(void);

#ifdef __cplusplus
extern "C" {
#endif
 int getpid(void);
 int time(int*);
#ifdef __cplusplus
}
#endif
#endif

double drand();
double nndist(int, double*, double*, int);
double packing(double);
double elastprob(int,double*,double*,int,double,double);
double filter(int);
char *emalloc(unsigned int n);
void sendxy(int n, double x, double y, int filout, int printfl);
void nndistm(int i, double *xv, double *yv, int n, double *nnd, int *nn);

/*--------------------------------------------------------------*/

int gausnn(double mean, double stdev, double density, double ms, int grseed, 
	double framex, double framey, double xcent, double ycent, 
	int numcells, double **xarr, double **yarr, int first_center, int filout, int textfl,
	int info)
                                                            
/* Creates array of points in xarr[] and yarr[] with
   distribution that is a Gaussian deviate with mean and
   stdev calculated from input parameters.  Returns address
   of xarr[] and yarr[] and number of (x,y) points generated.
*/

{
   int i,j,numpts,meanfl,initcells;
   int tries,tottries,itercrit,maxiter;
   int cn, iter;
   double p,okdist,mdist,mstd,fiter,initden,xcenter,ycenter;
   double cms,cmean,cstdev,n2mean;
   static double nnd[NUMN],*xval,*yval;
   static int nn[NUMN],varn[4];
   double sum[4],sumsq[4],val,area;
   double tvar=0.0, tstdev=0.0,tmean=0.0;
   double tstdevx=0.0, tmeanstd=0.0;
#ifndef NONC
   static char* gstate;
#endif

#define ROUNDUP 1e-12

   for (i=0; i<NUMN; i++) {
     nn[i] = 0;
     nnd[i] = 0;
   }
   for (i=0; i<4; i++) {
     varn[i] = 0;
   }
  
/* ncfprintf
(stderr,"m %g s %g d %g ms %g fx %g fy %g xc %g yc %g n %d f %d p %d\n",
      mean,stdev,density,ms,framex,framey,xcenter,ycenter,
		numcells,filout,info);			/* */

  if (!framex) framex = FRAMEX;
  if (!framey) framey = framex;
  area = framex * framey;
  if (xcent > 1e10) xcent = XCENT;
  if (ycent > 1e10) ycent = xcent;
  xcenter = xcent - framex * 0.5;		/* center of array */
  ycenter = ycent - framey * 0.5;

  /* Set the density and number of cells: */

  meanfl = (mean>0);
  if (numcells) {		/* Number specified overrides everything. */
    density = numcells / area; 
    n2mean = sqrt (1.0/density);/* approx avg neighbor mean */
    if (!stdev) {		/* either stdev or ms must be specified */
       double tmean;

       if (!ms) ms = 5;		/* default mean/stdev = 5 */
       stdev = n2mean / ms; 
       tmean = n2mean - stdev;
       stdev = tmean / ms; 
    }
    if (!meanfl)
       mean = n2mean - stdev;	/* Stdev must be specified. */ 
    ms = mean / stdev;		/* find approx mean/stdev */
    n2mean *= packing(ms);	/* correct avg neighbor mean */
    if (!meanfl)
	mean = n2mean - stdev;	/* final target nnd mean */
  }
  else if (density) {		/* density specified, not mean */
    numcells = int(((density * area) + 0.5 + ROUNDUP));
  }
  else if (meanfl) {		/* do nothing here because it's done below */
  }
  else {
    mean = 20;			/* default nnd == 20 um */
    meanfl = 1;
    if (!ms) ms = 5;		/* default mean/stdev = 5 */
  }

  /* Set mean and std. deviation: */

  if (density && ms) {		/* "mean/stdev" ratio specified */
    n2mean = sqrt(1.0/density)*packing(ms); /* avg neighbor mean */
    stdev = n2mean / (ms + 1.0);	/* find approx stdev */
    if (!meanfl) mean = n2mean - stdev;  /* find approx nnd mean */
  }
  else if (meanfl) {		/* only mean, stdev specified */
     if (!stdev) {
        if (ms) stdev = mean / ms;
        else {
	   ncfprintf (stderr,"gausnn: stdev must be specified with mean\n");
           return (0);
        }
     }
     if (!ms) ms = mean / stdev;
     n2mean = mean + stdev;	/* average neighbor dist = nnd + stdev */
     density = packing(ms)*packing(ms)/(n2mean*n2mean);
     if (!numcells) numcells = (int)((density * area) + 0.5 + ROUNDUP);
  }

/* Two adjustments are necessary in order
   to make mean and stdev correct. 

     1) Since this algorithm uses two neighbors to calculate
         probability of acceptance, stdev is less than it should be
         (if only one neighbor was used) by a factor of 2.  Therefore,
         to compensate, we must multiply stdev by a factor of 2.

     2) The nearest-neighbor mean is less than the "average-neighbor"
         mean (approximated here by the 2nd nearest neighbor) by a
         factor equal to the stdev.  Therefore to compensate we 
         increment the mean by the stdev.

   With these modifications, mean and stdev approach correct values
   for regular arrays (mean/stdev > 10).  But of course with less
   regular arrays we need to compensate both mean and stdev
   in other ways...
*/

  ms = mean / stdev;
  cmean = mean;
  cstdev = stdev;
  cms = cmean / cstdev;
/*  if (cms<10) cstdev *= pow((10/cms),0.0);	   /* fine tune stdev */
  cstdev *=  2.0; 			/* stdev created by this alg is low  */
  if (cms<10) cmean += cstdev; /* *pow((10/ms),0.0);	/* */
  else        cmean += cstdev; 		/* nnd is always less than 2nnd */


  if (numcells > MAXPTS){
	ncfprintf (stderr, "Gausnn: too many cells: %d\n",numcells);
	numcells = MAXPTS;
  }
if (filout) ncfprintf (textout,"Gausnn,  Sept 1992:\n\n");

textfl = 1;
if (textfl && (info >= 3)) {
ncfprintf (textout,"# Array size:     X=%g,  Y=%g, area=%-8.4g\n",
				framex,framey,area);
ncfprintf (textout,"# Target density: %-6.4g\n",density);
ncfprintf (textout,"# Target: n %4d mean %-6.4g stdev %-6.4g m/s ratio %-6.4g\n",
		numcells,mean,stdev,ms);
/*
ncfprintf(textout,"# Set:    n %4d mean %-6.4g stdev %-6.4g m/s ratio %-6.4g\n",
		numcells,cmean,cstdev,cmean/cstdev); /* */
}

if (filout) {
  if (outfil==NULL) {
     if ((fp=fopen (outfilnam,"w"))==NULL) {
        ncfprintf (stderr,"Can't open '%s'\n",outfilnam);
        outfil = stdout;
     }
     else outfil = fp;
  }
}

  tottries = 0;
  numpts = 0;

  if (numcells <= 0) {
      if (textfl && (info >= 3)) {
        ncfprintf (textout,"# Gausnn: can't allocate array with zero cells\n");
        ncfprintf (textout,"#  allocating 1 cell instead\n");
      }
      numcells = 1;
  }
  if ((xval=(double *)emalloc((size_t)(numcells*sizeof(double)))) == NULL) {
	ncfprintf (textout,"Gausnn: can't allocate xval size %d\n",numcells);
        return (0);
  }
  if ((yval=(double *)emalloc((size_t)(numcells*sizeof(double)))) == NULL) {
	ncfprintf (textout,"Gausnn: can't allocate yval size %d\n",numcells);
        return (0);
  }
 
  for (i=0; i<numcells; i++) {		/* zero the x,y locations */
        xval[i] = 0;
        yval[i] = 0;
  }

  if (first_center > 0) {			/* make first cell be in the center */
      i=0;
      xval[i] = framex*0.5 + xcenter; 		/* make one point in center */
      yval[i] = framey*0.5 + ycenter; 
      sendxy (i,xval[i],yval[i],filout,info>=3);
      numpts++;
  }
					/*  make some initial seed points */
					/*  but not too close */
  initcells = (int)((numcells/(ms*6) + 0.5 + ROUNDUP));  	/* */
  if (ms>15) initcells = 1; 		/* only one start if very regular */
  else if (initcells==0) initcells = 1;

#ifndef NONC	/* no need for separate random number gen if not in nc */

       /* set up random number generator with unique seed */

#define IMASK 0x7fffffff
  if (grseed >=0) {		/* use default random number gen if neg */
    if (grseed==0) {
      grseed = noduniq(numcells&IMASK, int(cmean*1e8)&IMASK, 
		int(cstdev*1e8)&IMASK,0,grseed^7272);
    }
  }
  else {
     if (grseed<0) grseed = (getpid()*65536+time((int*)NULL)) & IMASK;
  }
  gstate = (char *)makrand (RNDSIZ,"gausnn",numcells);
  if(gstate) initstate(((grseed&IMASK)+3297657)&IMASK,gstate,RNDSIZ);
#undef IMASK
  if (gstate) drand_setstate(gstate);

#endif

  initden = initcells / area;
  okdist = 0.82*sqrt(1.0 / initden);            
  for (i=numpts; i<initcells; numpts++,i++) { 
     for (mdist=0,tries=0; mdist<okdist && tries<1000; tries++) { 
        tottries++;
        xval[i] = drand() * framex + xcenter; 
        yval[i] = drand() * framey + ycenter; 
        mdist = nndist(i,xval,yval,numpts);
     }
     if (mdist<okdist) {i--; continue;}
     sendxy(i,xval[i],yval[i],filout,info>=4);
  }

		/* make maxiter proportional to number of cells */
		/*   and square of regularity */

/* MAXITER controls how many incorrect points per cell are
   not allowable.  ITERCRIT controls how many incorrect points per
   cell is allowable, as long as the next cell is found "quickly".
*/

#define MAXITER .5
#define ITERCRIT .10
#define FITERSIZ 7

  mstd = min (ms,50);
  maxiter  = (int)((MAXITER  * numcells * ms * mstd) + 0.5 + ROUNDUP);
  itercrit = (int)((ITERCRIT * numcells * ms * mstd) + 0.5 + ROUNDUP);

  for (i=numpts; i<numcells; numpts++,i++) {
     for (iter=1; iter<maxiter; iter++) {
        tottries++;
        xval[i] = drand() * framex + xcenter; 
        yval[i] = drand() * framey + ycenter; 
        p = elastprob(i,xval,yval,numpts,cmean,cstdev);
	/* ncfprintf (stderr,"n %d p %g\n",numpts,p); /* */
        if (p > drand()) break; 
     }
     if (iter >= maxiter) {
        if (info>=4) 
	 ncfprintf(textout,"Number of tries (%d) too many, stopping...\n",iter);
	   break;
     }
     sendxy(i,xval[i],yval[i],filout,info>=4);
     if (info>=4) ncfprintf (textout,"iter %d\n",iter);
     fiter = filter(iter);
     if (info>=4) 
        ncfprintf(textout,"threshold %d filter %g\n",itercrit,fiter);
     if (fiter >= itercrit ) {
        if (info>=4) 
           ncfprintf(textout,"Slowing down on last %d tries, stopping.\n",
				FITERSIZ);
        numpts++;
	break;
     }

  }  /* for (i;;) */

  for (j=0; j<2; j++) {			/* clear arrays for stdev calc. */
     sum[j] = 0;
     sumsq[j] = 0;
  }

  for (i=0; i<numpts; i++) { 
    nndistm(i,xval,yval,numpts,nnd,nn);		/* find nearest 4 neighbs */
    for (j=0; j<2; j++) {
      val = nnd[j];
      sum[j] += val;
      sumsq[j] += val*val;
    }
  } 
  varn[0] = numpts;
  varn[1] = numpts;
/*  varn[1] += varn[0];			/* */
/*  sum[1]  += sum[0];			/* add nearest to next neighbor dist */
/*  sumsq[1] += sumsq[0];			/* */

  for (j=0; j<2; j++) {
    if (varn[j] <= 1) varn[j] = 2;
    tmean = sum[j] / varn[j];
    tvar = (sumsq[j] - (sum[j] * tmean)) / (varn[j] - 1);
    tstdev = sqrt (tvar);
    tstdevx = tstdev;
    if (tstdevx == 0.0) tstdevx = 1.0;
    tmeanstd = tmean / tstdevx;
    if (textfl && (info>=3)) {
      ncfprintf (textout,"# nn #%d:  ",j+1);
      ncfprintf (textout,"n %4d mean %-6.4g stdev %-6.4g m/s ratio %-6.4g\n",
		varn[j],tmean,tstdev,tmeanstd);
    }
  }
 if (textfl && (info>=3)) {
   ncfprintf (textout,"# Final density   %-6.4g\n",numpts/area);
   ncfprintf (textout,"# Total tries     %d\n",tottries);
 }
#ifndef NONC
 if (gstate) restorstate();
#endif
 *xarr = xval;
 *yarr = yval;
 return (numpts);
#undef ROUNDUP
}

/* -------------------------------------------------------------- */

double filter(int ival)
            

/* Low-pass filter for trycount.
   Implemented with circular buffer.

   Uses individual weights for beginning and ending points,
   to reduce influence of single values (high frequencies).
*/

#define WT0 0.7
#define WT1 0.9
#define WT2 0.5
#define WT3 0.1

{
   static double tbuf[FITERSIZ] = {0};
   static double *inp=tbuf, *outp=tbuf;
   double outval;
   int i;

  for (i=0; i<FITERSIZ; i++) tbuf[i] = 0;
  *inp = ival;
  if (++inp >= tbuf+FITERSIZ) inp = tbuf;
  for (outval=0.0,outp=inp,i=0; i<FITERSIZ; i++) {
    switch (i) { 
      case 0:
          outval += *outp * WT0;
	  break;
      case (FITERSIZ-3):
          outval += *outp * WT1;
	  break;
      case (FITERSIZ-2):
          outval += *outp * WT2;
	  break;
      case (FITERSIZ-1):
          outval += *outp * WT3;
	  break;
      default:
    	  outval += *outp;
	  break;
    }  /* switch */
    if (--outp < tbuf) outp = tbuf+FITERSIZ-1;
  }   
  return (outval/(FITERSIZ-(4-(WT0+WT1+WT2+WT3))));   
}

/* -------------------------------------------------------------- */

double packing (double reg)

/* Find approximate packing radius.
   Weight sqrt ( 2 / sqrt(3)) == 1.07457 for triangular packing
   when regularity is high.  When packing is more
   square or random, reduce packing radius.
   The parameter "reg" is equal to mean / stdev of the
   nearest neighbor distance (not avg. neighbor distance). 
   This function is used by multiplying with the the average
   neighbor distance.  To get nearest neighbor distance from 
   avg. neighbor distance, subtract the nnd stdev.

   To get more cells or higher packing density, 
   increase values in "packdat[]".
*/

{
#define PACKSIZ 11
  double a,a10f,r,rfactor;
  static double packdat[PACKSIZ]= {1.04, 0.985, 0.975, 0.965, 0.955,
				   0.95, 0.94,  0.94,  0.94,  0.93, .92};
  int a10i;

 reg = max(reg,1);
 a = 2/reg;
 if (a<0) a = 0;
 else if (a>1) a=.99;
 a10f = a * 10;
 a10i = (int)(a10f);
 r = a10f - a10i;
 rfactor = (1-r) * packdat[a10i] + r * (packdat[a10i+1]);

/* ncfprintf (stderr,"a %g  a10i %d rfactor %g\n",a,a10i,rfactor); /* */
 return (rfactor);    
}

/* -------------------------------------------------------------- */

void sendxy(int n, double x, double y, int filout, int printfl)

/* print a cell's position */

{
  if (filout) {
    ncfprintf(outfil,"%g  %g\n",x,y);
    fflush (outfil);
  }
  if (printfl) ncfprintf (textout,"N %d %g %g\n",n+1,x,y);
}

/* -------------------------------------------------------------- */

double fgauss(double x, double mu, double sigma)
{
double value,r;
        r = (x-mu) / sigma;
        value = exp(-r*r);
        return(value);
}

/* ---------------------------------------------------------- */

void nndistm(int i, double *xv, double *yv, int n, double *nnd, int *nn)

/* Find nearest 4 neighbors, and their distances from point i */

{
int j;
double dist,xt,yt,md1,md2,md3,md4;
int n1=0,n2=0,n3=0,n4=0;

        for (md4=md3=md2=md1=1e10,j=0; j<n; j++){
                if (i!=j) {
                   xt = xv[i]-xv[j];
                   yt = yv[i]-yv[j];
                   dist = xt*xt + yt*yt;        
                   if (dist < md1) {
                         md4 = md3;
                         md3 = md2;
                         md2 = md1;
                         md1 = dist;
                         n4  = n3;
                         n3  = n2;
                         n2  = n1;
                         n1  = j;
                   }
                   else if (dist < md2) {
                         md4 = md3;
                         md3 = md2;
                         md2 = dist;
                         n4  = n3;
                         n3  = n2;
                         n2  = j;
                   }
                   else if (dist < md3) {
                         md4 = md3;
                         md3 = dist;
                         n4  = n3;
                         n3  = j;
                   }
                   else if (dist < md4) {
                         md4 = dist;
                         n4  = j;
                   }

                 }
        }
   nnd[3] = sqrt(md4);
   nnd[2] = sqrt(md3);
   nnd[1] = sqrt(md2);
   nnd[0] = sqrt(md1);
   nn[3]  = n4;
   nn[2]  = n3;
   nn[1]  = n2;
   nn[0]  = n1;
}

/* ---------------------------------------------------------- */

double nndist(int i, double *xv, double *yv, int n)

/* find nearest neighbor distance */

{
int j;
double dist,xt,yt,min_dist;

        for (min_dist=1e10,j=0; j<n; j++){
                if (i!=j){
                   xt = xv[i]-xv[j];
                   yt = yv[i]-yv[j];
                   dist = xt*xt + yt*yt;        
                   if (dist < min_dist) {
                         min_dist = dist;
                   }
                }
        }
return sqrt(min_dist);
}

/* ---------------------------------------------------------- */

double elastprob(int i, double *xv, double *yv, int n, 
		double emean, double estdev)

/* Find probability of acceptance. */

{
  int j;
  double prob;
  double nnd[NUMN];
  int nn[NUMN];

  nndistm(i,xv,yv,n,nnd,nn); /* find nearest 4 neighbs */
  for (prob=1.0,j=0; j<1; j++) {
         if (j>0 && nnd[j] > emean*2+estdev*2) continue; /* */
        /* if (j>0 && nnd[j] > 1000) continue; /* */
        prob *= fgauss(nnd[j],emean,estdev);
   /* ncfprintf (stderr,"n %d j %d nndist %g prob %9.4g %6.4g %6.4g\n",n,j,nnd[j],prob,emean,estdev);  /*  */
  }

/* for (prob=1.0,j=0; j<n; j++){
        if (i!=j) {
           xt = xv[i]-xv[j];
           yt = yv[i]-yv[j];
           dist = sqrt(xt*xt + yt*yt);  
           if (dist < emean) {
                prob *= fgauss(dist,emean,estdev);
 /* ncfprintf (stderr,"prob %g\n",prob);  /*  */
/*         }
        }
  }
*/

  return (prob);
}

/* ---------------------------------------------------------- */

