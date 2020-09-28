/* gaussfit */
/* 2D curve fitting to a Gaussian */

/* R.G. Smith   June, 2009 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ncfuncs.h"
#include "ncinit.h"

#include "lmmin.h"
#include "lm_eval.h"

double kc, kc_max, kc_min;
double rc, rc_max, rc_min;
double ks, ks_max, ks_min;
double rs, rs_max, rs_min;
double xo, xo_max, xo_min;
double yo, yo_max, yo_min;

/*-----------------------------------------------------------------------------*/

double gauss (double x, double r)

{
  double val;

  val = -(x*x/(r*r));
  if (val < -745.0) val = 0;
  else val = exp(val);
  return val;
}

/*-----------------------------------------------------------------------------*/

double gauss_2d (double x, double y, double *coeff, int modelnum) 

/* 2D Gaussian for curve fitting. 
  Coefficients are: 
	coeff[0] = KC (amplitude)
	coeff[1] = RC (radius)
	coeff[2] = KS (amplitude)
	coeff[3] = RS (radius)
	coeff[4] = X offset
	coeff[5] = Y offset
*/

#define _KC 0
#define _RC 1
#define _KS 2
#define _RS 3
#define _X0 4
#define _Y0 5
#define NCOEFF 6

{
    double r, rc, rs;
    double xr,yr, gaussval;

  xo = coeff[_X0];
  yo = coeff[_Y0];
  xr = x - xo;
  yr = y - yo;
  r = sqrt (xr*xr + yr*yr);

  rc = coeff[_RC];
  if (rc == 0) rc = 1e-30;

  rs = coeff[_RS];
  if (rs == 0) rs = 1e-30;
  kc = coeff[_KC];
  ks = coeff[_KS];

  gaussval = kc * gauss(r, rc) - ks * gauss(r, rs);
  return gaussval;
}

/*-----------------------------------------------------------------------------*/

#define DATASIZE  50  // size of 2D data array

double coeff[NCOEFF] = {0};
double *coeffc = 0;
double xydata[2*DATASIZE*DATASIZE];

main (int argc, char **argv) 
{
    int i, datasize, maxdata;
    double dx;
    const char *filenam;
    const char *progname;
    FILE *fp;

  progname = argv[0];
  if (argc<2) {
     fprintf (stderr,"Usage: gaussfit <parameters>\n");
     fprintf (stderr,"	     --kc n       center amplitude (1.0)\n");
     fprintf (stderr,"	     --rc n       center radius (5.0)\n");
     fprintf (stderr,"       --rc_max n   center rad max \n");
     fprintf (stderr,"       --rc_min n   center rad min \n");
     fprintf (stderr,"	     --ks n       surround amplitude\n");
     fprintf (stderr,"	     --rs n       surround radius \n");
     fprintf (stderr,"	     --x0 n       X offset \n");
     fprintf (stderr,"	     --y0 n       y offset \n");
     fprintf (stderr,"	     --dx n       spatial incr (1.0)\n");
     fprintf (stderr,"	     --file xxx   file name (file.dat)\n");
     fprintf (stderr,"	     --datasize n size of square array (15)\n");
     }

  setptr("kc",      &kc);
  setptr("kc_max",  &kc_max);
  setptr("kc_min",  &kc_min);
  setptr("rc",      &rc);
  setptr("rc_max",  &rc_max);
  setptr("rc_min",  &rc_min);
  setptr("ks",      &ks);
  setptr("ks_max",  &ks_max);
  setptr("ks_min",  &ks_min);
  setptr("rs",      &rs);
  setptr("rs_max",  &rs_max);
  setptr("rs_min",  &rs_min);
  setptr("xo",      &xo);
  setptr("xo_max",  &xo_max);
  setptr("xo_min",  &xo_min);
  setptr("yo",      &yo);
  setptr("yo_max",  &yo_max);
  setptr("yo_min",  &yo_min);

  setptr("dx",   &dx);
  setptr("info", &info);

  setptrn("file", &filenam);
  setptrn("datasize", &datasize);

  kc = 1.0;
  rc = 5;
  ks = 0.04;
  rs = 25;
  xo = 0;
  yo = 0;

  dx = 1;			// x increment for checkerboard
  info = 0;
  filenam = "file.dat";
  datasize = 15;

  setvars(argc,argv);		// set variables from command line

  if (xo==0) xo = (datasize-1)/2 * dx;
  if (yo==0) yo = (datasize-1)/2 * dx;

  for (i=0; i<datasize; i++) {
     xydata[i] = dx*i;		// xydata[0][0..datasize-1] -> X values
     xydata[datasize+i] = dx*i; // xydata[0][xsize..2*datasize-1] -> Y values
  }
  if (!(fp=fopen (filenam, "r"))) {
      fprintf (stderr,"%s: can't open file %s\n",progname,filenam);
      exit(0);	  
  }
  maxdata = datasize*datasize;
  for (i=0; i<maxdata; i++) {
       int n;
       double val;
     if ((n=fscanf(fp,"%lg ",&val))<1) {	// read in square array of data values from file
	fprintf (stderr,"%s: can't read data array\n",progname);
	break;
     }
     xydata[maxdata+i] = val;	// xydata[1][0..n] -> Z values to fit with least-squares 
  }

  // Starting coefficients:
  
  coeff[_KC] = kc;
  coeff[_RC] = rc;
  coeff[_KS] = ks;
  coeff[_RS] = rs;
  coeff[_X0] = xo;
  coeff[_Y0] = yo;

#define LC 0
#define UC 1

  coeffc = (double *)emalloc(NCOEFF*2*sizeof(double));

  coeffc[_KC*2+LC] = kc_min;
  coeffc[_KC*2+UC] = kc_max;
  coeffc[_RC*2+LC] = rc_min;
  coeffc[_RC*2+UC] = rc_max;
  coeffc[_KS*2+LC] = ks_min;
  coeffc[_KS*2+UC] = ks_max;
  coeffc[_RS*2+LC] = rs_min;
  coeffc[_RS*2+UC] = rs_max;
  coeffc[_X0*2+LC] = xo_min;
  coeffc[_X0*2+UC] = xo_max;
  coeffc[_Y0*2+LC] = yo_min;
  coeffc[_Y0*2+UC] = yo_max;
  
 lmfit2d (gauss_2d, maxdata, xydata, NCOEFF, coeff, coeffc);
}
