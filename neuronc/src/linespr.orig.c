
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

double blur_rad, scatter_rad1, scatter_rad2, scatter_rad3;
double ampl_scatter1, ampl_scatter2, ampl_scatter3;
double ampl_overall;
double gfact, total, peak;
double xpos,ypos;


xmid = (xstop - xstart) / 2;
incr = 1;

//blur_rad = 2.215;
//scatter_rad = 10.175;

s1 = 11.0;			/* 22 um dia Gaussian for cat eye */
sp = 9;

blur_rad = 9;
scatter_rad1 = 16;
scatter_rad2 = 36;
scatter_rad3 = 72;


gfact = 21;		/*  with    linespread normalization */
gfact = 1.0;		/*  without linespread normalization */
ampl_scatter1 = 0.20;
ampl_scatter2 = 0.025;
ampl_scatter3 = 0.002;

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
     double r1, r2, r3, r4, rad;
 
   xpos = (double)(x - xmid) / scale;
   ypos = (double)(y - xmid) / scale;
   rad = sqrt (xpos*xpos + ypos*ypos); 
   r1 = rad/blur_rad;
   r2 = rad/scatter_rad1;
   r3 = rad/scatter_rad2;
   r4 = rad/scatter_rad3;
   gf = ampl_overall * ( 1.0 * exp (-(r1*r1)) +
                ampl_scatter1 * exp (-(r2*r2))
	       	+ ampl_scatter2 * exp (-(r3*r3)) 
	       	+ ampl_scatter3 * exp (-(r4*r4)) 
		);

   sarr[x][y] = gf;

  };
 };

yo = 0;
for (yo=xstart; yo<xstop*3; yo+= incr) {
 for (y=xstart; y<xstop; y+= incr) {
  for (x=xstart; x<xstop; x+= incr) {
     double r1, r2, rad;
 
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
   r2 = rad/scatter_rad1;
   r3 = rad/scatter_rad2;
   rp = rad/sp;

   gg= 1.0 *exp(-(rg*rg)); 		/* Gaussian line spread function */
   g1= 1.0 *exp(-(r1*r1)); 		/* Gaussian line spread function */
   g2= ampl_scatter1 *exp(-(r2*r2)); 	/* Gaussian line spread function */
   g3= ampl_scatter2 *exp(-(r3*r3)); 	/* Gaussian line spread function */
   gp= 1 / (1 + pow(rp,2.5));

   gf = g2arr[x][y];

  printf ("%g %g %g %g %g %g\n",((double)x)/scale, gg, g2, g3, gp, gf*gfact);
  //printf ("%g %g %g %g %g\n", ((double)x)/scale, g1, g2, g3, gp - gf * gfact);
  //printf ("%g %g\n", ((double)x)/scale, gp - gf * gfact);
};

};
