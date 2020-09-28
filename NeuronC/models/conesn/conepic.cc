
/* conepic */

#include <stdio.h>
#include <math.h>

double xsize;
double ysize;
double ncones = 100;
double scale;
double maxgray = 255;
double midgray = 127;
double microns_degree = 300;
double cone_spac = 2.5;
double scontrast;
double speriod;
double spac_freq;
double envrad;

main (int argc, char **argv)

{
   int x,y;
   double midx, rad;
   double gaborfunc, sinefunc;

  scale = 5;
  spac_freq = 8;
  speriod = microns_degree/spac_freq; 
  envrad = speriod*3;
  scontrast = 1;
  ysize = .15 * microns_degree * scale;
  xsize = 1 * microns_degree * scale;
  printf ("P5\n%d %d\n%d\n",int(xsize),int(ysize),int(maxgray));

  midx = xsize / 2.0;
  for (y=0; y<ysize; y++) {
    for (x=0; x<xsize; x++) {
        sinefunc = scontrast * sin(2*M_PI * ((x-midx)/(speriod*scale)+.25));
        rad = (x-midx)/(envrad*scale);
        gaborfunc = exp(-rad*rad);
        printf ("%c", int((sinefunc*gaborfunc+1)*midgray));
     }
   }
}

