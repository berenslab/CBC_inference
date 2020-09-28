/* Module synfilt in program nc */

/* Synapatic filters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "ncsub.h"
#include "ncomp.h"
#include "control.h"
#include "ncio.h"

extern int cumlpfilt;

char *emalloc(unsigned int n);

void efree (void *ptr);

double ncabs(double val);

#ifdef __cplusplus
extern "C" {
#endif

#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif

#ifdef __cplusplus
}
#endif

double calctau(double tau, double timestep);
double makfilter(unsigned int options, int order, double raw_alpha1, double raw_alpha2, 
                        double *xcoeffs, double *ycoeffs);

/*------------------------------------*/

double *makfiltarr (int newsiz, int oldsiz, double *oarr, double val)

/* make larger filter array, copy old values to it */

{
    int i;
    double *narr;

  if (newsiz <= oldsiz) return((double *)NULL); /* do nothing if not larger */

  if ((narr=(double *)emalloc((newsiz)*sizeof(double))) == (double *)NULL) {
    ncfprintf (stderr,"makfiltarr: no space left for filter array\n");
    return ((double *)NULL);  
  }
  for (i=0; i<oldsiz; i++) {
    narr[i] = oarr[i]; 			/* copy the old values */
  }
  for (; i<newsiz; i++) {
    narr[i] = val; 			/* set the new spaces to val */
  }
  if (oarr) efree (oarr);		/* delete old filter array */
  return narr;
}

/*------------------------------------*/

double *copyfiltarr (int siz, double *oarr)

/* copy filter array to "save heap". */ 

{
    int i;
    double *narr;

  if ((narr=(double *)emalloc((siz)*sizeof(double))) == (double *)NULL) {
    ncfprintf (stderr,"copyfiltarr: no space left for filter array\n");
    return ((double *)NULL);  
  }
  for (i=0; i<siz; i++) {
    narr[i] = oarr[i]; 			/* copy the old values */
  }
  if (oarr) efree (oarr);		/* delete old filter array */
  return narr;
}

/*------------------------------------*/

// /* calctau now defined in "digfilt.cc" */
//
// double calctau(double tau, double timestep)
// 
// /* tau input in msec, output is multiplier for digital filter */
// /* M.v.Rossum, 1996 */
// 
// {
//    double k,taugran;
//    int sign;
// 
// #define TIMERES  (.1)
// #define MINTAU  (TIMERES)
// #define TAUGRAN (.632120558)  /* (1 - 1/e) */
// 
//  sign = (tau >= 0);
//  tau = ncabs(tau);
//  if (tau < timestep) tau =  timestep;
//  if (tau==0.0) tau = 1.0;
// 
// /*  taugran = timestep * pow(TAUGRAN,timestep/tau);  /* the old way */
// /*  k = taugran / tau; /* */
// 
//  k = 1.0 - exp(-timestep/tau);		/* fix error for small tau*/
//  if (k > 1.0) k = 1.0;			/* limit k to reasonable values */
//  else if (k<=0.0) k = .0001;
//  return ( sign ? k : -k);
// }

/*------------------------------------*/

lpfilt *maksynfilt (int stype, int nfilt, double offset, double tfall, double val, 
		double *timec, double timestep)

/* make a temporal filter */

{
     int i;
     lpfilt *fpnt;
     double *ftau,*lfilt;

  if (nfilt <= 0) return ((lpfilt*)NULL);
  if (nfilt > NUMFILT) nfilt = NUMFILT;
  nfilt++;				  /* make one more */
  if ((fpnt=(lpfilt *)emalloc(sizeof(lpfilt))) == (lpfilt *)NULL) {
      ncfprintf (stderr,"no space left for filter in synapse\n");
      return ((lpfilt *)NULL);  
  }
  lfilt = makfiltarr(nfilt, 0, (double*)NULL, 0.0);
  ftau  = makfiltarr(nfilt, 0, (double*)NULL, 0.0);
  fpnt->lfilt = lfilt;
  fpnt->ftau  = ftau;
  fpnt->offset  = offset;		/* offset for high-pass filter */
  fpnt->ctype = LPFILT;
  fpnt->stype = stype;			 /* lp or hp */
  fpnt->nfilt = nfilt;			 /* number of filters */
  fpnt->tfall = tfall;			 /* lp fall tau or hp gain */
  fpnt->xv = NULL;			 
  fpnt->yv = NULL;			 
  fpnt->gain = 0;			 
  for (i=0; i<nfilt-1; i++) {
 	 double timecn;
      if (timec) timecn = timec[i];	 /* check for missing timec */
      else       timecn = 1.0;
      fpnt->lfilt[i] = val;		 /* init low/high-pass filt array */
      fpnt->ftau[i] = calctau(timecn,timestep); /* time const, time step (ms)*/
  }
  fpnt->lfilt[i] = val;		 	/* init low/high-pass filt array */
  fpnt->ftau[i] = 0;			/* zero last one */
  cumlpfilt++;
  return (fpnt);
}

/*------------------------------------*/

lpfilt *makbessfilt (unsigned int options, int order, double alpha)

/* make a Bessel filter, given type, order, cutoff freq / sample freq = alpha*/
/* from: http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html */
{
     int i, nfilt;
     lpfilt *fpnt;
     double *xv,*yv;
     double *xcoeffs = NULL, *ycoeffs = NULL; 

#define opt_z  0x08000  /* use matched z-transform, defined in makfilter.cc  */

  nfilt = order + 1;				/* one more than the order (4th order bessel) */
  if ((fpnt=(lpfilt *)emalloc(sizeof(lpfilt))) == (lpfilt *)NULL) {
      ncfprintf (stderr,"no space left for filter in synapse\n");
      return ((lpfilt *)NULL);  
  }
  if (!(options & opt_z)) {		/* if not z-transform, make xv, xccoeffs */
      fpnt->xv = makfiltarr(nfilt, 0, (double*)NULL, 0.0);
  }
  else fpnt->xv = NULL;
  xcoeffs  = makfiltarr(nfilt, 0, (double*)NULL, 0.0);
  ycoeffs  = makfiltarr(nfilt, 0, (double*)NULL, 0.0);
  fpnt->yv = makfiltarr(nfilt, 0, (double*)NULL, 0.0);
  fpnt->gain = makfilter (options, order=4, alpha, alpha, xcoeffs, ycoeffs);
   // fprintf (stderr,"gain %g options %x\n",fpnt->gain, options);
   // for (i=0; i<nfilt; i++) { fprintf (stderr,"xcoeffs %d %g\n",i,xcoeffs[i]);}
   // for (i=0; i<nfilt; i++) { fprintf (stderr,"ycoeffs %d %g\n",i,ycoeffs[i]);}
  fpnt->lfilt = xcoeffs;
  fpnt->ftau  = ycoeffs;
  fpnt->offset  = 0;			/* offset for high-pass filter */
  fpnt->ctype = BESSFILT;
  fpnt->stype = 0;			 /* lp or hp */
  fpnt->nfilt = nfilt;			 /* number of filters */
  fpnt->tfall = 0;			 /* lp fall tau or hp gain */
  cumlpfilt++;
  return (fpnt);
}

/*------------------------------------*/

/* Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
 *    Command line: /www/usr/fisher/helpers/mkfilter -Be -Lp -o 4 -a 4.0000000000e-02 0.0000000000e+00 -l */

/* designed for 50K/sec sampling rate (20e-6), 2000 Hz cutoff */
/* see: http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html */
/* from: http://www-users.cs.york.ac.uk/~fisher/software/mkfilter/current */

#define NZEROS 4
#define NPOLES 4
#define GAIN 1.330668083e+03

// static double xv[NZEROS+1] = {0};
// static double yv[NPOLES+1] = {0};

double bessfilt (lpfilt *lpnt, double val)
{
  int i, order;
  double *xv, *yv;
  double *xcoeffs, *ycoeffs;
  double gain;

  order   = lpnt->nfilt - 1; 
  xv      = lpnt->xv; 
  xcoeffs = lpnt->lfilt;
  yv      = lpnt->yv; 
  ycoeffs = lpnt->ftau;
  gain    = lpnt->gain;

  if (xv!=NULL) {   /* bilinear transform method */

  //  xv[0] = xv[1]; xv[1] = xv[2]; xv[2] = xv[3]; xv[3] = xv[4];
  //  xv[4] = val / GAIN;
  //  yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3]; yv[3] = yv[4];
  //  yv[4] =   (xv[0] + xv[4]) + 4 * (xv[1] + xv[3]) + 6 * xv[2]
  //              + ( -0.3041592568 * yv[0]) + (  1.5960375869 * yv[1])
  //              + ( -3.1910200543 * yv[2]) + (  2.8871176889 * yv[3]);
  // return yv[4];

    for (i=0; i<order; i++) xv[i] = xv[i+1];
    xv[order] = val / gain;
    for (i=0; i<order; i++) yv[i] = yv[i+1];
    yv[order] = 0;
    for (i=0; i<lpnt->nfilt; i++) yv[order] += xv[i] * xcoeffs[i];
    for (i=0; i<order; i++)       yv[order] += yv[i] * ycoeffs[i];

  } else {

   /* matched z-transform method */

  // xv[0] = next input value / GAIN;
  // yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3]; yv[3] = yv[4]; 
  // yv[4] =   xv[0]
  //            + ( -0.3045507035 * yv[0]) + (  1.5983910456 * yv[1])
  //            + ( -3.1956640387 * yv[2]) + (  2.8901773622 * yv[3]);
  // next output value = yv[4];

    for (i=0; i<order; i++) yv[i] = yv[i+1];
    yv[order] = val / gain;
    for (i=0; i<order; i++) yv[order] += yv[i] * ycoeffs[i];
  }
  return yv[order];
}

#undef NZEROES
#undef NPOLES
#undef GAIN

/*------------------------------------*/

double synfilt(lpfilt *lpnt, double transrate)

/* Make spnt->nfilt first-order low pass filters.
   Each filter has a time constant of spnt->ftau[i]/10 * stiminc (0.1 msec) */
{
   int i,prev,nfilt;
   double ctr,k,tr,oldtr;

  if (!lpnt) return transrate;
  switch (lpnt->ctype) {

  case LPFILT:
   tr=transrate;
   if (!(nfilt=lpnt->nfilt-1))		/* number of low pass stages   */
	return tr;			/*  leave one for output       */   
   oldtr = lpnt->lfilt[0];
   lpnt->lfilt[0] = tr;			/* zero filter holds input to filter */
   for (i=1; i<nfilt; i++) {		/* do 1 less filter than actual num */
     prev = i-1;
     k=lpnt->ftau[prev];
     if (synaptau > 0.0) k /= synaptau;	/* possibly modify time constant */
     if (k >= 1.0) k = .9999;		/* limit tau to reasonable value */
     ctr = lpnt->lfilt[i];
     if (lpnt->stype==LP) {
 	if (k >= 0) ctr += (tr - ctr) * k;   	/* low pass filter */
	else { 					/* neg k, high pass filter */
	  ctr += tr - oldtr;
          ctr *= (1+k);
          oldtr = lpnt->lfilt[i];
        }
     }
     else {				/* high pass */
        ctr += tr - oldtr;
        ctr *= (1-k);
        oldtr = lpnt->lfilt[i];
     }
     lpnt->lfilt[i] = tr = ctr; 
   }
		/* do the last one differently: separate falling phase */

   if (lpnt->stype==LP) {
     prev = i-1;
     k=lpnt->ftau[prev];
     if (synaptau > 0.0) k /= synaptau;	/* possibly modify time constant */
     if (k >= 1.0) k = .9999;		/* limit tau to reasonable value */
     ctr = lpnt->lfilt[i];
     if (lpnt->tfall>0.0) {		/* separate falling phase */
	 double trel;

       trel = (tr - ctr) * k; 	 	/* transmitter released */
       if (trel > 0) ctr += trel; 	/* rising phase, trel is non-neg */
       ctr *= lpnt->tfall;		/* falling phase */
       tr = ctr;
    }
    else {
 	if (k >= 0) ctr += (tr - ctr) * k;	/* neg k, low pass filter */
	else { 					/* high pass filter */
	  ctr += tr - oldtr;
          ctr *= (1+k);
        }
    }
   tr = ctr; 
   lpnt->lfilt[i] = tr;
   return (ctr);				/* LP: return last filter val */
  }			/* if LP */

  else {		/* if HP */
    prev = i-1;
    k=lpnt->ftau[prev];
    if (synaptau > 0.0) k /= synaptau;	/* possibly modify time constant */
    if (k >= 1.0) k = .9999;		/* limit tau to reasonable value */
    ctr = lpnt->lfilt[i];
    ctr += tr - oldtr;
    ctr *= (1-k);
    lpnt->lfilt[i] = ctr;
    return (ctr);		/* HP: tfall is gain */
  }
  break;

  case BESSFILT: return (bessfilt(lpnt,transrate)); break;
   
  case MFILT: return (0);
  break; 
  }
 return 0;
}

/*------------------------------------*/
