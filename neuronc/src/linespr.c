
/* 
 Curve fit of 2-D point spread function based on line-spread 
 function.

 Fit to line-spread function for white light taken from either
  
  1) Robson & Enroth-Cugell, (1978) 4 mm pupil, cat.

  2) Campbell & Gubisch, 1966,  5.7 mm pupil, human.

  3) Guirao, et al. (2001), 5.8 mm pupil, "best refracted" average human eye.

  4) Guirao, et al (2001), 5.8 mm pupil, "uncorrected" average human eye.


 Form of function taken from Robson and Enroth-Cugell (1978).
 
        I (r) = 1 / (1 + (r/k)^(5/2))
 
*/

#include <stdio.h>
#include <math.h>

#define minutes_degree 60
#define microns_degree 300

#define scale 3		/* need at least 3-5 for accurate line-spread func */
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
double g1, g2, g3, g4, gf;
double blur_rad;
double scatter_rad; 
double scatter_ampl;
double ppow;
double overall_ampl;

double total, peak;
double xpos,ypos;
double oscale;


xmid = (xstop - xstart) / 2;
incr = 1;

s1 = 2.215;

/* from Robson and Enroth-Cugell (1978) */
/*  for a 4 mm pupil, cat */

//scatter_rad  = 22.55;
//ppow         = 2.5;
//scatter_ampl = .15;
//blur_rad     = 11;

/* From Campbell and Gubisch (1966) */
/*   For a 5.7 mm pupil.  See "cg58.txt". */

//blur_rad    = 1.0;      /* no Gaussian blur here */
//scatter_rad = 2.7;
//ppow        = 1.95;

scatter_rad  = 2.7;
ppow         = 1.95;
scatter_ampl = 10.0;
blur_rad     = 1.5; 

/* For "best refracted" average human eye, from Guirao et al, (2001). */
/*   For a 5.7 mm pupil.  See "refr.txt". */

// scatter_rad  = 1.3;
// ppow         = 1.85;
// scatter_ampl = 1000.0;
// blur_rad     = 1.5; 

/* for "uncorrected" average human eye, from Guirao et al, (2001). */
/*    for a 5.7 mm pupil.  See "uncorr.txt". */

// scatter_rad  = 1.5;
// ppow         = 1.0;
// scatter_ampl = 1000.0;
// blur_rad     = 1.5; 

oscale = 1.0;
overall_ampl = 1.0;

/* initialize the arrays */

 for (y=xstart; y<yarrsiz; y+= incr) {
  for (x=xstart; x<xstop; x+= incr) {
     sarr[x][y] = 0;
    g2arr[x][y] = 0;
  };
 };

  /* Fill up source array with single blur function */

  /* For a 5.8 mm pupil the Gaussian "diffraction image" is much smaller */
  /*  the "scatter image" so the curve is best fit with a blur function  */
  /*  like the Enroth-Cugell & Robson scatter function. */

 for (y=xstart; y<xstop; y+= incr) {
  for (x=xstart; x<xstop; x+= incr) {
     double r1, r2, r3, rg, rs, gg, rad;
 
   xpos = (double)(x - xmid) / scale;
   ypos = (double)(y - xmid) / scale;
   rad = sqrt (xpos*xpos + ypos*ypos); 
   rs = rad/scatter_rad;
   rg = rad/blur_rad;
   gg = exp (-(rg*rg));
   gf= scatter_ampl * 1 / (1 + pow(rs,ppow));
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

/* Normalize line spread volume to 1 */
/*  Do it here to make it similar to Geisler's (1984) fit */
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

y = xstop;	/* take section at least one blur's worth into conv arr */
ypos = 0;
for (x=xstart; x<xarrsiz; x+= incr) {
     double rg, r1, r2, r3, rs, rad;
     double gg, gp;

   xpos = (double)(x - xmid) / scale;
   rad = sqrt (xpos*xpos + ypos*ypos); 
   rg = rad/s1;
   r1 = rad/blur_rad;
   rs = rad/scatter_rad;

   gg= 1.0 *exp(-(rg*rg)); 		/* Gaussian line spread function */
   g1= 1.0 *exp(-(r1*r1)); 		/* Gaussian line spread function */
   gp= 1 / (1 + pow(rs,1.5));

   gf = g2arr[x][y];

  //printf ("%g %g\n",((double)x)/scale, gf);
  //printf ("%g %g\n",((double)x)/scale/microns_degree, gf);
  printf ("%g %g\n",((double)(x-xmid))/scale/microns_degree, gf*oscale);
  //printf ("%g %g %g %g\n",((double)x)/scale, gg, gp, gf);
  //printf ("%g %g %g\n", ((double)x)/scale, g1, gp - gf);
  //printf ("%g %g\n", ((double)x)/scale, gp - gf);
};

};
