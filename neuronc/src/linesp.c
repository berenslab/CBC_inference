
/* 
 Curve fit of 1-D line spread function based on fit
 of Campbell and Gubisch (1966) by Geisler (1984)
 to an equivalent 2-D point spread function.
 Original walues are in terms of minutes of visual angle.

 Courtesy of A. Hsu, 1996

 This is implemented in "pointspfov()" inside "stimsub.cc".

 The lesson learned here is that the radii of the Gaussian functions
  are not changed by the transformation to a point spread function,
  but their weighting is.
 
*/

#include <stdio.h>
#include <math.h>

#define minutes_degree 60
#define microns_degree 300

#define scale 5
#define xwidth (30*scale)
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
double g1, g2;

double blur_rad, scatter_rad, scatter_radr;
double ampl_scatter, ampl_overall;
double gfact, total, peak;
double xpos,ypos;

a1 = 0.684;
b1 = 0.443;
s1 = 0.443 / minutes_degree * microns_degree;  /* = 2.215 um radius */
a2 = 0.587;
b2 = 2.035;
s2 = 2.035 / minutes_degree * microns_degree;  /* = 10.175 um radius */


xmid = (xstop - xstart) / 2;
incr = 1;


//blur_rad = 2.215;
//scatter_rad = 10.175;

blur_rad = 2.215;
scatter_rad = 10.175;

ampl_scatter = (a2/b2)  / (a1/b1);
ampl_overall = (a1/(2*b1) + a2/(2*b2)) / (1+ampl_scatter);

//print ampl_scatter, ampl_overall;

//ampl_scatter = 0.18681912;	/* for line spread function */
//ampl_overall = 0.77200903;
ampl_overall = 1.0;

//gfact = .146364;		/* without linespread normalization */
gfact = 7.783265;		/*  with   linespread normalization */

scatter_radr = 10.157/2.215;	/* ratio of radii */
ampl_scatter = 1 / (scatter_radr * scatter_radr * 1.1077382); /* 0.04278007; */

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
     double r1, r2, rad;
 
   xpos = (double)(x - xmid) / scale;
   ypos = (double)(y - xmid) / scale;
   rad = sqrt (xpos*xpos + ypos*ypos); 
   r1 = rad/blur_rad;
   r2 = rad/scatter_rad;
   g2 =     ampl_overall * ( 1.0 * exp (-0.5*(r1*r1)) +
                 ampl_scatter * exp (-0.5*(r2*r2)));

   sarr[x][y] = g2;

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

/* Normalize line spread volume to 1 */
/*  Do it here to make it similar to Geisler's fit */
/*  But, no need to make absolute amplitude correct because */
/*  blur volume is always normalized before use. */

/* 
 total = 0;
  y = xstop;
  for (x=xstart; x<xarrsiz; x+= incr) {
    total += g2arr[x][y];
  };
  total /= scale;
  for (y=xstart; y<yarrsiz; y+= incr) {
   for (x=xstart; x<xstop; x+= incr) {
    g2arr[x][y] /= total;
   };
  };
*/

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
     double r1, r2, rad;

   xpos = (double)(x - xmid) / scale;
   rad = sqrt (xpos*xpos + ypos*ypos); 
   r1 = rad/s1;
   r2 = rad/s2;
   g1=(a1/(2*b1))*exp(-0.5*(r1*r1)) +	/* original line spread function */
      (a2/(2*b2))*exp(-0.5*(r2*r2));	/* from Geisler, (1984) */

   g2 = g2arr[x][y];

  printf ("%g %g\n", ((double)(x-xmid))/scale/microns_degree, g2);
  //printf ("%g %g %g\n", ((double)x)/scale, g1, g2 * gfact);
  //printf ("%g %g\n", ((double)x)/scale, g1 - g2 * gfact);
};

};
