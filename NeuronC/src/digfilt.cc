/* module digfilt in program nc */
/* digital filters for general use */

/* RGS 8/2009 */

#include <stdio.h>
#include <math.h>

/*------------------------------------*/

double calctau(double tau, double timestep)

 /* tau input in msec, output is multiplier for digital filter */
 /* M.v.Rossum, 1996 */
 /* Called from "synfilt.cc" so this links with retsim. */

{
   double k,taugran;
   int sign;

#define TIMERES  (.1)
#define MINTAU  (.01)
#define TAUGRAN (.632120558)  /* (1 - 1/e) */

   sign = (tau >= 0);
   if (tau<0) tau = -tau;
   if (tau < (timestep*MINTAU)) tau =  timestep*MINTAU;
   if (tau==0.0) tau = 1.0;

   /*  taugran = timestep * pow(TAUGRAN,timestep/tau);  /* the old  way */
   /*  k = taugran / tau; /* */
   
   k = 1.0 - exp(-timestep/tau);          /* fix error for small tau*/
   if (k > 1.0) k = 1.0;                  /* limit k to reasonable values */
   else if (k<=0.0) k = 1e-8;
   // fprintf (stderr,"k %g\n",k);
   return ( sign ? k : -k);
}

/*------------------------------------*/

static double oldval1=0;
static double oldval2=0;
static double digfilt_k1 = 1;
static double digfilt_k2 = 1;


void init_digfilt (double init_val, double tau, double timestep)
{
   digfilt_k1 = calctau(tau, timestep);
   oldval1 = init_val;
}

/*------------------------------------*/

void init_digfilt2 (double init_val, double tau, double timestep)
{
   digfilt_k2 = calctau(tau, timestep);
   oldval2 = init_val;
}

/*------------------------------------*/

double digfilt (double val)

{
   oldval1 += (val - oldval1) * digfilt_k1;          /* low pass filter */
   return oldval1;
}

/*------------------------------------*/

double digfilt2 (double val)

{
   oldval2 += (val - oldval2) * digfilt_k2;          /* low pass filter */
   return oldval2;
}

/*------------------------------------*/

/* Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Be -Lp -o 4 -a 4.0000000000e-02 0.0000000000e+00 -l */

/* designed for 50K/sec sampling rate (20e-6), 2000 Hz cutoff */

#define NZEROS 4
#define NPOLES 4
#define GAIN 1.330668083e+03

static double xv[NZEROS+1] = {0};
static double yv[NPOLES+1] = {0};

double bessfilt (double val)
{
  xv[0] = xv[1]; xv[1] = xv[2]; xv[2] = xv[3]; xv[3] = xv[4]; 
  xv[4] = val / GAIN;
  yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3]; yv[3] = yv[4]; 
  yv[4] =   (xv[0] + xv[4]) + 4 * (xv[1] + xv[3]) + 6 * xv[2]
               + ( -0.3041592568 * yv[0]) + (  1.5960375869 * yv[1])
               + ( -3.1910200543 * yv[2]) + (  2.8871176889 * yv[3]);
  return yv[4];
}

#undef NZEROES
#undef NPOLES
#undef GAIN

/*------------------------------------*/

/* 
main (int argc, char **argv)

{
  double t, exptdur, scontrast, meanval;

  exptdur = 1;

  meanval = 1;
  tau = 0.01;
  tstep = 0.001;
  scontrast = 0.3;

  init_digfilt  (meanval, tau, tstep);
  init_digfilt2 (meanval, tau, tstep);

  for (t=0; t<exptdur; t+= tstep){
        double start, dur, inten;

    inten = digfilt2(digfilt((gasdev()*scontrast + meanval)));
    if (inten < 0) inten = 0;
    printf ("%g\n",inten);
  }
}
*/

