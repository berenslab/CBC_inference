
/* makgauss */
/* output 2D Gaussian */

/* R.G. Smith   July, 2009 */


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

  xr = x - coeff[X0];
  yr = y - coeff[Y0];
  r = sqrt (xr*xr + yr*yr);
  rc = coeff[RC];
  rs = coeff[RS];
  gaussval = coeff[KC] * gauss(r, rc) - coeff[KS] * gauss(r, rs);
  return gaussval;
}

/*-----------------------------------------------------------------------------*/

#define DATASIZE  50  // size of 2D data array

double coeff[NCOEFF] = {0};

main (int argc, char **argv) 
{
    int i, j, k, datasize, maxdata;
    double dx,gaussval;
    double kc, rc;
    double ks, rs;
    double x0, y0;
    const char *filenam = "file.dat";
    const char *progname;
    FILE *fp;
    char *cptr;

  progname = argv[0];
  if (argc<2) {
     fprintf (stderr,"Usage: %s <parameters>\n",progname);
     fprintf (stderr,"	     --kc n       center amplitude (1.0)\n");
     fprintf (stderr,"	     --rc n       center radius (5.0)\n");
     fprintf (stderr,"	     --ks n       surround amplitude (0.04)\n");
     fprintf (stderr,"	     --rs n       surround radius (25.0)\n");
     fprintf (stderr,"	     --x0 n       X offset \n");
     fprintf (stderr,"	     --y0 n       y offset \n");
     fprintf (stderr,"	     --dx n       spatial incr (1)\n");
     fprintf (stderr,"	     --file xxx   file name (file.dat)\n");
     fprintf (stderr,"	     --datasize n size of square array (50)\n");
     exit(0);
   }

  setptrn("kc",   &kc);
  setptrn("rc",   &rc);
  setptrn("ks",   &ks);
  setptrn("rs",   &rs);
  setptrn("x0",   &x0);
  setptrn("y0",   &y0);
  setptrn("dx",   &dx);

  setptrn("file", &filenam);
  setptrn("datasize", &datasize);

  kc = 1.0;
  rc = 5;
  ks = 0.04;
  rs = 25;
  dx = 1;			// x increment for checkerboard
  x0 = 0;
  y0 = 0;
  filenam = "file.dat";		// not used
  datasize = 15;

  setvars(argc,argv);		// set variables from command line

  if (x0==0) x0 = (datasize-1)/2 * dx;
  if (y0==0) y0 = (datasize-1)/2 * dx;

  coeff[KC] = kc;
  coeff[RC] = rc;
  coeff[KS] = ks;
  coeff[RS] = rs;
  coeff[X0] = x0;
  coeff[Y0] = y0;

  maxdata = datasize*datasize;
  
    for (j=0; j<datasize; j++) {
      for (i=0; i<datasize; i++) {
	gaussval = gauss_2d(i,j,coeff);
        printf ("%8.6g ",gaussval);
      } 
      printf ("\n");
    } 

}
