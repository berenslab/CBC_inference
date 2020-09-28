
/* 
 Curve fit of 2-D line spread function based on line-spread 
 function of the cat eye, according to Robson and Enroth-Cugell (1978).
 
        I (r) = 1 / (1 + (r/k)^(5/2))
 
   where r = 4.4 um matches a Gaussian of radius 6.

*/

#include <stdio.h>
#include <math.h>

#define minutes_degree 60
#define microns_degree 300

#define scale 1
#define xwidth (100*scale)
#define xstop (2*xwidth)
#define xarrsiz xstop
#define yarrsiz (xstop*4)
#define xstart 0		/* microns */

double sarr[xarrsiz][yarrsiz];
double g2arr[xarrsiz][yarrsiz];

main (int argc, char **argv) {

int x,y,yo;
int incr, xmid;
double a1, b1, s1;
double a2, b2, s2;
double g1, g2, g3, g4, gf;
double sp;

double blur_rad; 
double ppow;
double ampl_overall;
double total, peak;
double xpos,ypos;

double ampl_scatter;

xmid = (xstop - xstart) / 2;
incr = 1;

s1 = 11.0;			/* 22 um dia Gaussian for cat eye */
sp = s1 * .82;

blur_rad = sp * 2.5;
ppow = 2.5;

ampl_scatter = .15;

ampl_overall = 1.0;

/* initialize the arrays */

 for (y=xstart; y<yarrsiz; y+= incr) {
  for (x=xstart; x<xstop; x+= incr) {
     sarr[x][y] = 0;
    g2arr[x][y] = 0;
  };
 };

  /* Fill up source array with single blur function */

  /* The idea is, here we assume the smaller (dia) of the 2 Gaussians */
  /*  is the "blur" function which always starts out with amplitude of 1. */
  /*  The larger (dia) Gaussian is the "scatter" function which must be */
  /*  set to the correct relative amplitude */

 for (y=xstart; y<xstop; y+= incr) {
  for (x=xstart; x<xstop; x+= incr) {
     double r1, r2, r3, rg, rp, gg, rad;
 
   xpos = (double)(x - xmid) / scale;
   ypos = (double)(y - xmid) / scale;
   rad = sqrt (xpos*xpos + ypos*ypos); 
   rp = rad/blur_rad;
   rg = rad/sp;
   gg = exp (-(rg*rg));
   gf= ampl_scatter * 1 / (1 + pow(rp,ppow));
   gg += gf;

   sarr[x][y] = gg;

  };
 };

yo = 0;
for (yo=xstart; yo<xstop*3; yo+= incr) {
 for (y=xstart; y<xstop; y+= incr) {
  for (x=xstart; x<xstop; x+= incr) {
 
   g2arr[x][y+yo] += sarr[x][y];

  };
 };
};

/* Normalize line spread peak to 1 */

  peak = g2arr[xmid][xstop];
  for (y=xstart; y<yarrsiz; y+= incr) {
   for (x=xstart; x<xstop; x+= incr) {
    g2arr[x][y] /= peak;
   };
  };

y = xstop;
ypos = 0;
for (x=xstart; x<xarrsiz; x+= incr) {
     double rg, r1, r2, r3, rp, rad;
     double gg, gp;

   xpos = (double)(x - xmid) / scale;
   rad = sqrt (xpos*xpos + ypos*ypos); 
   rg = rad/s1;
   r1 = rad/blur_rad;
   rp = rad/sp;

   gg= 1.0 *exp(-(rg*rg)); 		/* Gaussian line spread function */
   g1= 1.0 *exp(-(r1*r1)); 		/* Gaussian line spread function */
   gp= 1 / (1 + pow(rp,1.5));

   gf = g2arr[x][y];

  printf ("%g %g %g %g\n",((double)x)/scale, gg, gp, gf);
  //printf ("%g %g %g\n", ((double)x)/scale, g1, gp - gf);
  //printf ("%g %g\n", ((double)x)/scale, gp - gf);
};

};
