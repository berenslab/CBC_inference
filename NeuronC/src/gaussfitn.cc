/* gaussfit */
/* 2D curve fitting to a Gaussian */

/* R.G. Smith   June, 2009 */


#include <stdio.h>
#include <math.h>
#include <unistd.h>

//#include "ncfuncs.h"

#include "ncio.h"
#include "lmmin.h"
#include "lm_eval.h"
#include "lm_funcs.h"

typedef struct STYPE {
       const char *name;
       int type;
       union {
            double val;
            int    *iptr;
            char   *cptr;
            char   **sptr;
            float  *fptr;
            double *dptr;
      };
} stype;

#define VVSIZE 100
stype varval[VVSIZE];
int varset=0;
int info=0;

double kc, kc_max, kc_min;
double rc, rc_max, rc_min;
double ks, ks_max, ks_min;
double rs, rs_max, rs_min;
double xo, xo_max, xo_min;
double yo, yo_max, yo_min;

void setptrn (const char *name, int    *iptr);
void setptrn (const char *name, float  *fptr);
void setptrn (const char *name, double *dptr);
void setptrn (const char *name, const char **sptr);
int setvars(int argc, char **argv);

#ifdef __cplusplus
extern "C" {
#endif
double atof(char *ptr);
void exit(int n);
#ifdef __cplusplus
}
#endif


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

double gauss_2d (double x, double y, double *coeff) 

/* 2D Gaussian for curve fitting. 
  Coefficients are: 
	coeff[0] = KC (amplitude)
	coeff[1] = RC (radius)
	coeff[2] = KS (amplitude)
	coeff[3] = RS (radius)
	coeff[4] = X offset
	coeff[5] = Y offset
*/

#define KC 0
#define RC 1
#define KS 2
#define RS 3
#define X0 4
#define Y0 5
#define NCOEFF 6

{
    double r, rc, rs;
    double xr,yr, gaussval;


  /* check against max,min range */

  xo = coeff[X0];
  if (xo > xo_max) xo = xo_max;
  if (xo < xo_min) xo = xo_min;
  yo = coeff[Y0];
  if (yo > yo_max) yo = yo_max;
  if (yo < yo_min) yo = yo_min;
  xr = x - xo;
  yr = y - yo;
  r = sqrt (xr*xr + yr*yr);

  rc = coeff[RC];
  if (rc > rc_max) rc = rc_max;
  if (rc < rc_min) rc = rc_min;
  if (rc == 0) rc = 1e-30;
  rs = coeff[RS];
  if (rs > rs_max) rs = rs_max;
  if (rs < rs_min) rs = rs_min;
  if (rs == 0) rs = 1e-30;

  kc = coeff[KC];
  if (kc > kc_max) kc = kc_max;
  if (kc < kc_min) kc = kc_min;
  ks = coeff[KS];
  if (ks > ks_max) ks = ks_max;
  if (ks < ks_min) ks = ks_min;

  gaussval = kc * gauss(r, rc) - ks * gauss(r, rs);
  return gaussval;
}

/*-----------------------------------------------------------------------------*/

#define DATASIZE  50  // size of 2D data array

double coeff[NCOEFF] = {0};
double xydata[2*DATASIZE*DATASIZE];

main (int argc, char **argv) 
{
    int i, j, k, datasize, maxdata;
    double dx;
    const char *filenam = "file.dat";
    const char *progname;
    FILE *fp;
    char *cptr;

  progname = argv[0];
  if (argc<2) {
     fprintf (stderr,"Usage: %s <parameters>\n",progname);
     fprintf (stderr,"	     --kc n       center amplitude (1.0)\n");
     fprintf (stderr,"	     --rc n       center radius (5.0)\n");
     fprintf (stderr,"	     --rc_max n   center rad max \n");
     fprintf (stderr,"	     --rc_min n   center rad min \n");
     fprintf (stderr,"	     --ks n       surround amplitude\n");
     fprintf (stderr,"	     --rs n       surround radius \n");
     fprintf (stderr,"	     --xo n       X offset \n");
     fprintf (stderr,"	     --yo n       y offset \n");
     fprintf (stderr,"	     --dx n       spatial incr (1.0)\n");
     fprintf (stderr,"	     --file xxx   file name (file.dat)\n");
     fprintf (stderr,"	     --datasize n size of square array (15)\n");
     exit(0);
   }

  setptrn("kc",      &kc);
  setptrn("kc_max",  &kc_max);
  setptrn("kc_min",  &kc_min);
  setptrn("rc",      &rc);
  setptrn("rc_max",  &rc_max);
  setptrn("rc_min",  &rc_min);
  setptrn("ks",      &ks);
  setptrn("ks_max",  &ks_max);
  setptrn("ks_min",  &ks_min);
  setptrn("rs",      &rs);
  setptrn("rs_max",  &rs_max);
  setptrn("rs_min",  &rs_min);
  setptrn("xo",      &xo);
  setptrn("xo_max",  &xo_max);
  setptrn("xo_min",  &xo_min);
  setptrn("yo",      &yo);
  setptrn("yo_max",  &yo_max);
  setptrn("yo_min",  &yo_min);

  setptrn("dx",   &dx);
  setptrn("info", &info);

  setptrn("file", &filenam);
  setptrn("datasize", &datasize);

  kc = 1.0;  kc_max = 1e6; kc_min = -1e6;
  rc = 5;    rc_max = 1e6; rc_min = -1e6;
  ks = 0.04; ks_max = 1e6; ks_min = -1e6;
  rs = 25;   rs_max = 1e6; rs_min = -1e6;
  xo = 0;    xo_max = 1e6; xo_min = -1e6;
  yo = 0;    yo_max = 1e6; yo_min = -1e6;

  dx = 1;			// x increment for checkerboard
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
  /*
  for (k=0; k<2; k++) {
    for (j=0; j<datasize; j++) {
      for (i=0; i<datasize; i++) {
        fprintf (stderr,"%6.3g ",xydata[k*maxdata+j*datasize+i]);
      } 
      fprintf (stderr,"\n");
    } 
  } 
  */

  // Starting coefficients:
  
  coeff[KC] = kc;
  coeff[RC] = rc;
  coeff[KS] = ks;
  coeff[RS] = rs;
  coeff[X0] = xo;
  coeff[Y0] = yo;
  
 lmfit2d (gauss_2d, maxdata, xydata, NCOEFF, coeff);
}
